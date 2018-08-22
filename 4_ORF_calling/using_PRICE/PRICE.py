#!/usr/bin/env python
import traceback
import sys
import os
import sqlite3 as sqlite
import time
import argparse
from collections import defaultdict
import re
from multiprocessing import Pool


__author__ = 'Steven Verbruggen'


def main():

    starttime = time.time()

    print
    print "##########################"
    print "# ORF Calling with PRICE #"
    print "##########################"
    print
    print "Info about PRICE: https://github.com/erhard-lab/gedi/wiki/Price"
    print

    #Parse arguments from command line
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS,
                                     description="This tool is part of the PROTEOFORMER pipeline. It calls translated "
                                                 "ORFs with the use of the PRICE tool.\n"
                                                 "!! PRICE only works in combination with STAR mapping. Generate PRICE "
                                                 "specific alignment files during mapping.",
                                     add_help=False)

    man_args = parser.add_argument_group("Mandatory parameters")
    man_args.add_argument("--result_db", "-r", action="store", required=True, nargs="?", metavar="DB", type=str,
                          help="The SQLite results database")

    opt_args = parser.add_argument_group("Optional parameters")
    opt_args.add_argument("--help", "-h", action="help", help="Show help message and exit")
    opt_args.add_argument("--workdir", "-w", action="store", required=False, nargs="?", metavar="FOLDER", default="",
                        type=str, help="Working directory (default: CWD)")
    opt_args.add_argument("--tmp", "-t", action="store", required=False, nargs="?", metavar="FOLDER", default="",
                           type=str, help="Temporary folder (default: workdir/tmp)")
    opt_args.add_argument("--fdr", "-f", action="store", required=False, nargs="?", metavar="float", default=0.1,
                          type=float, help="FDR value for running PRICE (default: 0.1)")
    opt_args.add_argument("--max_ram", "-m", action="store", required=False, nargs="?", metavar="integer", default=0,
                          type=int, help="Maximum amount of gigabytes available for PRICE (default: Price default)")

    args = parser.parse_args()

    #default workdir is CWD, default tmp, default fdr
    if args.workdir == "":
        args.workdir = os.getcwd()
    if args.tmp == "":
        args.tmp = args.workdir+"/tmp/"
    if not os.path.isdir(args.tmp):
        os.mkdir(args.tmp)
    price_tmp = args.tmp+"PRICE/"
    if not os.path.isdir(price_tmp):
        os.mkdir(price_tmp)
    if args.fdr == "":
        args.fdr = 0.1

    #Get parameters from arguments table
    table_args = argparse.Namespace()
    table_args.igenomes_root, table_args.species, table_args.ens_v, table_args.ensdb, table_args.readtype, \
        table_args.bam_untr, table_args.bam_tr, table_args.cores, table_args.tr_calling_method, table_args.run_name, \
        table_args.mapper, table_args.uniq = get_arguments(args.result_db)

    # Conversion of species terminology
    speciesLatin = "Mus_musculus" if table_args.species == "mouse" else \
        "Homo_sapiens" if table_args.species == "human" else \
        "Arabidopsis_thaliana" if table_args.species == "arabidopsis" else \
        "Drosophila_melanogaster" if table_args.species == "fruitfly" else ""
    speciesShort = "mmu" if table_args.species == "mouse" else \
        "hsa" if table_args.species == "human" else \
        "ath" if table_args.species == "arabidopsis" else \
        "dme" if table_args.species == "fruitfly" else ""
    # Assembly
    assembly = "GRCm38" if table_args.species == "mouse" and table_args.ens_v >= 70 else \
        "NCBIM37" if table_args.species == "mouse" and table_args.ens_v < 70 else \
        "GRCh38" if table_args.species == "human" and table_args.ens_v > 75 else \
        "GRCh37" if table_args.species == "human" and table_args.ens_v <= 75 else \
        "TAIR10" if table_args.species == "arabidopsis" else \
        "BDGP5" if table_args.species == "fruitfly" else ""

    #List parameters
    print "Parameters:"
    for arg in vars(args):
        print '    %-15s\t%s' % (arg, getattr(args, arg))
    print "Parameters from arguments table:"
    for arg in vars(table_args):
        print '    %-15s\t%s' % (arg, getattr(table_args, arg))
    print
    sys.stdout.flush()

    #Check if bam files are indexed
    index_bams(table_args.bam_untr, table_args.bam_tr)

    #Define gtf file path
    gtf_file = table_args.igenomes_root+"/"+speciesLatin+"/Ensembl/"+assembly+"/Annotation/Genes/genes_"+str(table_args.ens_v)+".gtf"

    #Get chromosomes
    chrs={}
    print "Get chromosomes"
    print
    chromosomeSizesFile = table_args.igenomes_root+"/"+speciesLatin+"/Ensembl/"+assembly+"/Annotation/Genes/ChromInfo.txt"
    if os.path.isfile(chromosomeSizesFile):
        chrs = get_chrs(chromosomeSizesFile, table_args.species)
    else:
        chromosomeSizesFile = table_args.igenomes_root + "/" + speciesLatin + "/Ensembl/" + assembly + "/Sequence/WholeGenomeFasta/GenomeSize.xml"
        if os.path.isfile(chromosomeSizesFile):
            chrs = get_chrsXml(chromosomeSizesFile, table_args.species)
        else:
            print "ERROR: chromosome sizes file could not be found."
            sys.exit()

    #Create binary chromosomes if they don't exist
    BIN_chrom_dir = table_args.igenomes_root+"/"+speciesLatin+"/Ensembl/"+assembly+"/Sequence/Chromosomes_BIN"
    if os.path.isdir(BIN_chrom_dir):
        print "Binary chromosomes already exist"
        print
    else:
        print "Create binary chromosomes"
        print
        create_BIN_chrom(table_args.cores,BIN_chrom_dir,chrs,table_args.igenomes_root,speciesLatin,assembly)
    sys.stdout.flush()


    ########
    # MAIN #
    ########


    #Init
    price_prefix = price_tmp+"output"
    all_orfs_tsv = price_prefix+".orfs.tsv"
    cit_file = price_prefix+".orfs.cit"
    selected_orfs_loc = price_prefix+".orfs.loc"

    #Check if file not already exist
    if not (os.path.isfile(all_orfs_tsv) and os.path.isfile(selected_orfs_loc)):

        #Install PRICE
        print
        print "[%s]: Install PRICE" % (convert_time(time.time()-starttime))
        sys.stdout.flush()
        price_operatable, price_install_folder = install_price(price_tmp)
        print "[%s]: COMPLETE: PRICE installation" % (convert_time(time.time() - starttime))
        sys.stdout.flush()

        #Prepare genome reference
        print
        print "[%s]: Prepare PRICE genome reference" % (convert_time(time.time() - starttime))
        sys.stdout.flush()
        fasta = table_args.igenomes_root+"/"+speciesLatin+"/Ensembl/"+assembly+"/Sequence/WholeGenomeFasta/genome.fa"
        gtf_file = table_args.igenomes_root + "/" + speciesLatin + "/Ensembl/" + assembly + "/Annotation/Genes/genes_" + str(
            table_args.ens_v) + ".gtf"
        gref = prepare_genome_price(price_operatable, fasta, gtf_file, speciesLatin.lower(), str(table_args.ens_v))
        print "[%s]: COMPLETE: PRICE genome reference" % (convert_time(time.time() - starttime))
        sys.stdout.flush()

        #Run PRICE
        print
        print "[%s]: Run PRICE" % (convert_time(time.time() - starttime))
        sys.stdout.flush()
        run_price(price_operatable, price_prefix, gref, table_args.readtype, table_args.bam_untr, table_args.bam_tr,
                  table_args.cores, args.fdr, price_tmp, args.max_ram)
        print "[%s]: COMPLETE: PRICE run" % (convert_time(time.time() - starttime))
        sys.stdout.flush()

        #Run ViewCIT
        print
        print "[%s]: Convert CIT file" % (convert_time(time.time() - starttime))
        sys.stdout.flush()
        convert_cit(cit_file, selected_orfs_loc, price_operatable)
        print "[%s]: COMPLETE: CIT file conversion" % (convert_time(time.time() - starttime))
        sys.stdout.flush()

        #Delete PRICE installation
        #os.system("rm -rf "+price_install_folder)

    else:
        print "PRICE output files already found"
        sys.stdout.flush()

    #Parse the PRICE output files
    print
    print "[%s]: Parse PRICE output" % (convert_time(time.time() - starttime))
    sys.stdout.flush()
    orfs = parse_price_output(all_orfs_tsv, selected_orfs_loc)
    print "[%s]: COMPLETE: PRICE output parsing" % (convert_time(time.time() - starttime))
    sys.stdout.flush()

    ##Get all data needed for PROTEOFORMER transcripts table
    print
    print "[%s]: Convert to assembly table features" % (convert_time(time.time() - starttime))
    sys.stdout.flush()
    tr_features = construct_tr_features(orfs, table_args.ensdb, gtf_file, BIN_chrom_dir, args.result_db, price_tmp,
                                        table_args.run_name, table_args.mapper, table_args.uniq, table_args.ens_v)
    print "[%s]: COMPLETE: assembly table features" % (convert_time(time.time() - starttime))
    sys.stdout.flush()

    #Drop possible existing table and initiate new table, Save data in TIS_overview table
    print
    print "[%s]: Construct SQLite table" % (convert_time(time.time() - starttime))
    sys.stdout.flush()
    tis_id = get_id(args.result_db, table_args.tr_calling_method, args.fdr)
    table_name = "TIS_"+str(tis_id)+"_transcripts"
    create_new_table(args.result_db, table_name)

    #Construct csv
    csv_name = price_tmp+"/assembly_table.csv"
    filtered_out_csv = price_tmp+"/filtered_out.csv"
    count_orfs = construct_csv(csv_name, filtered_out_csv, tr_features)
    print "Amount of ORFs found by PRICE and saved to assembly table: "+str(count_orfs)
    sys.stdout.flush()
    #Dump
    dump_into_sqlite(args.result_db, table_name, csv_name)
    print "[%s]: COMPLETE: SQLite table" % (convert_time(time.time() - starttime))
    sys.stdout.flush()

    #Remove tmp files
    #os.system("rm -rf "+price_tmp)


    #End of program message
    print
    print "-----------------------"
    print "[%s]: PROGRAM COMPLETE" % (convert_time(time.time()-starttime))
    print "-----------------------"
    print
    sys.stdout.flush()


    return

## Dump into sqlite
def dump_into_sqlite(results_db, table_name, csv_name):

    #Dump
    try:
        os.system("sqlite3 -separator , "+results_db+" \".import "+csv_name+" "+table_name+"\"")
    except:
        print "ERROR: dumping of file "+csv_name+" into database "+results_db+" failed!"
        sys.exit()

    return

## Construct csv dump file
def construct_csv(csv_name, filtered_out_csv, tr_features):

    counter=0

    with open(csv_name, 'w') as FW:
        FiltFW = open(filtered_out_csv, 'w')
        for price_id in tr_features:
            if tr_features[price_id]['start_codon']!="Unmatching_start_codon" and tr_features[price_id]['aa_seq'].split(',')[0]!="early_stop" and \
                    tr_features[price_id]['aa_seq'] != "no_stop":
                line = tr_features[price_id]['tr_stable_id']+","+\
                    tr_features[price_id]['chr']+","+\
                    str(tr_features[price_id]['strand'])+","+\
                    str(tr_features[price_id]['start'])+","+\
                    tr_features[price_id]['start_codon']+","+\
                    str(tr_features[price_id]['stop'])+","+\
                    tr_features[price_id]['starts_list']+","+\
                    tr_features[price_id]['ends_list']+","+\
                    str(tr_features[price_id]['dist_to_transcript_start'])+","+\
                    str(tr_features[price_id]['dist_to_aTIS'])+","+\
                    tr_features[price_id]['annotation']+","+\
                    tr_features[price_id]['biotype']+","+\
                    tr_features[price_id]['aTIS_call']+","+\
                    str(tr_features[price_id]['peak_shift'])+","+\
                    str(tr_features[price_id]['count'])+","+\
                    str(tr_features[price_id]['Rltm_min_Rchx'])+","+\
                    str(tr_features[price_id]['coverage'])+","+\
                    str(tr_features[price_id]['FPKM'])+","+\
                    tr_features[price_id]['SNP']+","+\
                    tr_features[price_id]['INDEL']+","+\
                    tr_features[price_id]['secs']+","+\
                    tr_features[price_id]['tr_seq']+","+\
                    tr_features[price_id]['aa_seq']+","+\
                    str(price_id)+","+\
                    tr_features[price_id]['price_annotation']+","+\
                    tr_features[price_id]['price_pval']+"\n"
                counter+=1
                FW.write(line)
            else:
                line = tr_features[price_id]['tr_stable_id']+","+ \
                    tr_features[price_id]['chr'] + "," + \
                    str(tr_features[price_id]['strand']) + "," + \
                    str(tr_features[price_id]['start']) + "," + \
                    tr_features[price_id]['start_codon'] + "," + \
                    str(tr_features[price_id]['stop']) + "," + \
                    tr_features[price_id]['starts_list'] + "," + \
                    tr_features[price_id]['ends_list'] + "," + \
                    tr_features[price_id]['annotation'] + "," + \
                    str(tr_features[price_id]['dist_to_transcript_start']) + "," + \
                    str(tr_features[price_id]['dist_to_aTIS']) + "," + \
                    tr_features[price_id]['tr_seq'] + "," + \
                    tr_features[price_id]['aa_seq'] + "," + \
                    tr_features[price_id]['secs'] + "," + \
                    str(price_id)+"\n"
                FiltFW.write(line)
        FiltFW.close()

    return counter

## Create new ORF table
def create_new_table(result_db, table_name):

    #Connect to sqlite db
    con = sqlite.connect(result_db)

    with con:
        cur = con.cursor()

        #Delete existing table
        drop_query = "DROP TABLE IF EXISTS '"+table_name+"';"
        cur.execute(drop_query)

        #Create new table
        create_query = "CREATE TABLE IF NOT EXISTS '"+table_name+"' (" \
                        "'tr_stable_id' varchar(128) NOT NULL default ''," \
                        "'chr' char(50) NOT NULL default ''," \
                        "'strand' int(2) NOT NULL default ''," \
                        "'start' int(10) NOT NULL default ''," \
                        "'start_codon' varchar(128) NOT NULL default ''," \
                        "'stop' int(10) NOT NULL default ''," \
                        "'starts_list' varchar(512) NOT NULL default ''," \
                        "'ends_list' varchar(512) NOT NULL default ''," \
                        "'dist_to_transcript_start' int(10) NOT NULL default ''," \
                        "'dist_to_aTIS' int(10) NOT NULL default ''," \
                        "'annotation' varchar(128) NOT NULL default ''," \
                        "'biotype' varchar(128) NOT NULL default ''," \
                        "'aTIS_call' varchar(128) NOT NULL default ''," \
                        "'peak_shift' int(2) NOT NULL default ''," \
                        "'count' float default NULL," \
                        "'Rltm_min_Rchx' decimal(11,8) NOT NULL default '0'," \
                        "'coverage' decimal(11,8) NOT NULL default '0'," \
                        "'FPKM' decimal(11,8) NOT NULL default '0'," \
                        "'SNP' varchar(256) NOT NULL default '0'," \
                        "'INDEL' varchar(256) NOT NULL default '0'," \
                        "'secs' varchar(256) NOT NULL default ''," \
                        "'tr_seq' TEXT NOT NULL default ''," \
                        "'aa_seq' TEXT NOT NULL default '',"\
                        "'price_id' varchar(256) NOT NULL default '',"\
                        "'price_annotation' varchar(256) NOT NULL default ''," \
                        "'price_pval' varchar(256) NOT NULL default '');"
        cur.execute(create_query)

    return

## Get TIS ID
def get_id(result_db, tr_calling_method,fdr):

    #Connect to db
    try:
        conn = sqlite.connect(result_db)
    except:
        print "ERROR: could not connect to "+result_db
        sys.exit()

    with conn:
        cur = conn.cursor()

        #Create table if not yet exists
        query_create = "CREATE TABLE IF NOT EXISTS 'TIS_overview' (" \
                        "'ID' INTEGER primary key, " \
                        "'local_max' int(10) NOT NULL default '', " \
                        "'min_count_aTIS' int(10) NOT NULL default '', " \
                        "'R_aTis' decimal(11,8) NOT NULL default '', " \
                        "'min_count_5UTR' int(10) NOT NULL default '', " \
                        "'R_5UTR' decimal(11,8) NOT NULL default '', " \
                        "'min_count_CDS' int(10) NOT NULL default '', " \
                        "'R_CDS' decimal(11,8) NOT NULL default '', " \
                        "'min_count_3UTR' int(10) NOT NULL default '', " \
                        "'R_3UTR' decimal(11,8) NOT NULL default '', " \
                        "'min_count_ntr' int(10) NOT NULL default '', " \
                        "'R_ntr' decimal(11,8) NOT NULL default '', " \
                        "'PRICE_FDR' decimal(11,8) NOT NULL default ''," \
                        "'SPECTRE_FDR' decimal(11,8) NOT NULL default ''," \
                        "'SNP' varchar(20) default '', " \
                        "'indel' varchar(20) default '', " \
                        "'filter' varchar(20) default '', " \
                        "'tr_calling' varchar(20) default '', " \
                        "'TIS_calling' varchar(20) default '');"
        cur.execute(query_create)

        #Add parameters and get ID
        query_add = "INSERT INTO 'TIS_overview'(filter,tr_calling,TIS_calling,SNP,indel, PRICE_FDR) "\
                    "VALUES('none','"+tr_calling_method+"','PRICE','NO','NO','"+str(fdr)+"');"
        cur.execute(query_add)
        TIS_id = int(cur.lastrowid)

    return TIS_id

## Construct all necessary translation features for orf table
def construct_tr_features(orfs, ens_db, gtf_file, BIN_chrom_dir, result_db, price_tmp, run_name, mapper, uniq, ensversion):

    #Init
    tr_features = defaultdict(lambda: defaultdict())

    #Get underlying info from reads table, gtf and ensembl
    secs = get_secs(gtf_file)
    store_secs_in_results_db(result_db, secs, price_tmp)
    biotypes = get_biotypes(ens_db)
    reads_for, reads_rev = get_reads(result_db)
    sample_name = run_name+"_"+mapper+"_"+uniq+"_"+str(ensversion)+".fastq1"
    seq_depth = get_seq_depth(result_db, sample_name)
    seq_depthM = seq_depth/1000000

    for price_id in orfs:
        #Get some easy to obtain values
        tr_features[price_id]['tr_stable_id'] = orfs[price_id]['tr_stable_id']
        tr_features[price_id]['chr'] = orfs[price_id]['chr']
        tr_features[price_id]['strand'] = orfs[price_id]['strand']
        tr_features[price_id]['start_codon'] = orfs[price_id]['start_codon']
        #Get start and stop coordinates
        (tr_features[price_id]['start'], tr_features[price_id]['stop'], tr_features[price_id]['starts_list'],
            tr_features[price_id]['ends_list'], coding_parts) = get_start_and_stop_coordinates(orfs, price_id, tr_features[price_id]['strand'])
        #Get underlying info from Ensembl based on stable id: Ensembl transcript id and exons
        (transcript_id, exons) = get_transcript_id_and_exons(tr_features[price_id]['tr_stable_id'], ens_db)
        #Get canonical start and stop information
        (aTIS, stop) = get_aTIS(transcript_id, tr_features[price_id]['strand'], ens_db)
        #Calculate distances between the predicted TIS, the aTIS and the transcript start
        (tr_features[price_id]['dist_to_transcript_start'], tr_features[price_id]['dist_to_aTIS']) = \
            calc_dists(tr_features[price_id]['start'], aTIS, exons, tr_features[price_id]['strand'])
        tr_features[price_id]['annotation'] = get_annotation(tr_features[price_id]['start'], aTIS, stop,
                                                             tr_features[price_id]['strand'])
        #Construct sequence features
        (tr_features[price_id]['start_codon'], tr_features[price_id]['tr_seq'], tr_features[price_id]['aa_seq'],
            tr_features[price_id]['secs']) = construct_seq(coding_parts, tr_features[price_id]['strand'], BIN_chrom_dir,
            tr_features[price_id]['chr'], tr_features[price_id]['tr_stable_id'], price_id, secs,
            tr_features[price_id]['start_codon'])
        # Check if unmatching start codon, no stop or early stop
        if tr_features[price_id]['start_codon'] == "Unmatching_start_codon" or tr_features[price_id]["aa_seq"].split(',')[
            0] == "early_stop" or tr_features[price_id]["aa_seq"] == "no_stop":
            continue
        #Fetch biotype from biotypes ensembl info
        tr_features[price_id]['biotype'] = biotypes[transcript_id]
        # aTIS call and peak shift
        tr_features[price_id]['aTIS_call'] = "NA"
        tr_features[price_id]['peak_shift'] = "NA"
        #Calc coverage based on read table info
        if(tr_features[price_id]['strand']=="1"):
            (tr_features[price_id]['coverage'], tr_features[price_id]['count'], tr_features[price_id]['FPKM']) = \
                calc_coverage_count_fpkm(reads_for, coding_parts, tr_features[price_id]['chr'],
                                         len(tr_features[price_id]['tr_seq'])*1.0/1000, seq_depthM)
        elif(tr_features[price_id]['strand']=="-1"):
            (tr_features[price_id]['coverage'], tr_features[price_id]['count'], tr_features[price_id]['FPKM']) = \
                calc_coverage_count_fpkm(reads_rev, coding_parts, tr_features[price_id]['chr'],
                                         len(tr_features[price_id]['tr_seq'])*1.0/1000, seq_depthM)
        #We do not need the Rltm-Rchx value, as the filtering is done by PRICE itself
        tr_features[price_id]['Rltm_min_Rchx'] = "NA"
        #SNP and indel empty for the moment
        tr_features[price_id]['SNP'] = ""
        tr_features[price_id]['INDEL'] = ""
        #PRICE specific features
        tr_features[price_id]['price_pval'] = orfs[price_id]['pval']
        tr_features[price_id]['price_annotation'] = orfs[price_id]['price_annotation']

    return tr_features


## Determine coverage
def calc_coverage_count_fpkm(reads, coding_parts, chr, length_orfK, seq_depthM):

    #Init
    covered = 0
    not_covered = 0
    coverage = 0
    count = 0

    #Go over exon positions and find covered positions
    for exon_part in range(1, len(coding_parts)+1):
        for pos in range(coding_parts[exon_part]['start'], coding_parts[exon_part]['end']+1):
            if pos in reads[chr].keys():
                covered+=1
                count+=reads[chr][pos]
            else:
                not_covered+=1

    #Calc coverage
    if covered!=0 or not_covered!=0:
        coverage = float(covered) / (float(covered) + float(not_covered))

    #Calc fpkm
    fpm = count*1.0/seq_depthM
    fpkm = fpm*1.0/length_orfK

    return (coverage, count, fpkm)

## Get untr reads
def get_reads(result_db):

    #Init
    reads_for = defaultdict(lambda: defaultdict())
    reads_rev = defaultdict(lambda: defaultdict())

    #Connect
    try:
        con = sqlite.connect(result_db)
    except:
        print "ERROR: could not connect to "+result_db
        sys.exit()

    with con:
        cur = con.cursor()
        query = "SELECT start, chr, count FROM count_fastq1 WHERE strand = '1';"
        cur.execute(query)
        output = cur.fetchall()
        for i in range(0, len(output)):
            reads_for[output[i][1]][output[i][0]] = int(output[i][2])

        query = "SELECT start, chr, count FROM count_fastq1 WHERE strand = '-1';"
        cur.execute(query)
        output = cur.fetchall()
        for i in range(0, len(output)):
            reads_rev[output[i][1]][output[i][0]] = int(output[i][2])

    return reads_for, reads_rev

## Get biotypes of all transcripts in Ensembl
def get_biotypes(ens_db):

    #Init
    biotypes = defaultdict()

    try:
        conn = sqlite.connect(ens_db)
    except:
        print "Could not connect to Ensembl DB"
        sys.exit()

    with conn:
        cur = conn.cursor()
        query = "SELECT transcript_id, biotype FROM transcript;"
        cur.execute(query)
        output = cur.fetchall()
        for row in output:
            biotypes[row[0]] = row[1]

    return biotypes

## Get sequece of start codon, ORF cDNA
def construct_seq(coding_parts, strand, BIN_chrom_dir, chr, tr_stable_id, price_id, secs, start_codon):

    #Init
    tr_seq = ""
    aa_seq = ""
    sec_str = ""
    secs_in_seq = []
    #Go over all coding parts, situated over all exons
    for exon_part in coding_parts:
        ##If sec in between start and stop of coding part, method to know position of sec in translated coordinates.
        # Based on the length of the tmp tr_seq you already have, to know the coordinates of the exon parts already built in the translated sequence
        if tr_stable_id in secs.keys():
            for sec_id in secs[tr_stable_id].keys():
                #Check if sec in that coding part
                if secs[tr_stable_id][sec_id]['start']>coding_parts[exon_part]['start'] and secs[tr_stable_id][sec_id]['end']<coding_parts[exon_part]['end']:
                    pos_sec = len(tr_seq) + secs[tr_stable_id][sec_id]['start'] - coding_parts[exon_part]['start'] + 1
                    secs_in_seq.append(pos_sec)
        #Add sequence of coding part to tmp tr seq
        tr_seq += get_seq(coding_parts[exon_part]['start'], coding_parts[exon_part]['end'], chr, BIN_chrom_dir)
    #For reverse strand: take reverse complement in the end
    if strand == "-1":
        tr_seq = reverse_complement(tr_seq)
        converted_pos_secs = []
        for pos_sec in secs_in_seq:
            converted_pos_secs.append(len(tr_seq)-pos_sec+1-2) # -2 because the position of the start (T) is 2 positions earlier than the end of the TGA codon (A)
        secs_in_seq = converted_pos_secs
    # Check start codon with start codon from PRICE
    if (start_codon!=tr_seq[:3]):
        start_codon = "Unmatching_start_codon"
    else:
        #Get amino acid sequence
        aa_seq, sec_str = translate_seq(tr_seq, tr_stable_id, price_id, secs_in_seq)
    return start_codon, tr_seq, aa_seq, sec_str

## Translate sequence to amino acids
def translate_seq(tr_seq, tr_stable_id, price_id, secs_in_seq):

    #Init
    aa_seq = ""
    sec_str = ""

    # The codon table represents the corresponding amino acid for every codon.
    codonTable = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
                  'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
                  'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                  'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                  'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                  'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                  'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                  'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
                  'AAN': 'X', 'ATN': 'X', 'ACN': 'T', 'AGN': 'X', 'ANN': 'X', 'ANA': 'X', 'ANG': 'X', 'ANC': 'X',
                  'ANT': 'X', 'TAN': 'X', 'TTN': 'X', 'TCN': 'S', 'TGN': 'X', 'TNN': 'X', 'TNA': 'X', 'TNG': 'X',
                  'TNC': 'X', 'TNT': 'X', 'CAN': 'X', 'CTN': 'L', 'CCN': 'P', 'CGN': 'R', 'CNN': 'X', 'CNA': 'X',
                  'CNG': 'X', 'CNC': 'X', 'CNT': 'X', 'GAN': 'X', 'GTN': 'V', 'GCN': 'A', 'GGN': 'G', 'GNN': 'X',
                  'GNA': 'X', 'GNG': 'X', 'GNC': 'X', 'GNT': 'X', 'NAN': 'X', 'NAA': 'X', 'NAG': 'X', 'NAC': 'X',
                  'NAT': 'X', 'NTN': 'X', 'NTA': 'X', 'NTG': 'X', 'NTC': 'X', 'NTT': 'X', 'NGN': 'X', 'NGA': 'X',
                  'NGG': 'X', 'NGC': 'X', 'NGT': 'X', 'NCN': 'X', 'NCA': 'X', 'NCG': 'X', 'NCC': 'X', 'NCT': 'X',
                  'NNN': 'X', 'NNA': 'X', 'NNG': 'X', 'NNC': 'X', 'NNT': 'X'}

    #Go over sequence per triplet
    codon=""
    for pos in range(0, len(tr_seq)-2, 3):
        codon = tr_seq[pos:pos+3]
        #Check for selenocysteine (Do not build in a sec if we are at the end of the tr_seq)
        if((pos+1 in secs_in_seq) and (codon=="TGA") and (pos!=len(tr_seq)-3)):
            aa_seq+="U"
            sec_str = sec_str+str(secs_in_seq[pos+1])+"_"
        else:
            aa_seq+=codonTable[codon]
            #Control if stop codon is only at the end
            if codonTable[codon]=='*':
                if pos!=(len(tr_seq)-3):
                    #print "Error: stop earlier than expected (transcript ID: "+tr_stable_id+", PRICE ID: "+str(price_id)+")"
                    aa_seq = "early_stop,"+aa_seq
                    #print "tr seq: "+tr_seq
                    #print "AA seq: "+aa_seq
                    return aa_seq
    #Check if end codon is a STOP codon
    if codon!='TAA' and codon!='TAG' and codon!='TGA':
        #print "Error: no stop codon (TAA, TAG or TGA) found at the end (transcript ID: "+tr_stable_id+", PRICE ID: "+str(price_id)+")"
        #print "tr seq: " + tr_seq
        #print "AA seq: " + aa_seq
        aa_seq="no_stop"
        return aa_seq
    #Check for near-cognate starts and replace the first amino acid for methionine then
    if re.search('[ACTG]TG|A[ACTG]G|AT[ACTG]', tr_seq[:3]):
        aa_seq = "M" + aa_seq[1:]
    #print "tr seq: "+tr_seq
    #print "aa seq: "+aa_seq

    return aa_seq, sec_str[:-1]

## Get sequence based on start and end coordinate
def get_seq(start, stop, chromosome, BIN_chrom_dir):

    #Descirbe fasta file path
    BIN_chrom_fasta = BIN_chrom_dir+"/"+chromosome+".fa"

    #Convert to coordinates used by file reader
    offset = start - 1
    length = stop - start + 1

    with open(BIN_chrom_fasta, 'r') as FR:
        FR.seek(offset)
        seq = FR.read(length)

    return seq

#Give reverse complement of sequence
def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c':'g', 'g':'c', 't':'a'}
    return ''.join([complement[base] for base in dna[::-1]])

## Load all selenocysteines
def get_secs(gtf_file):

    #Init
    secs=defaultdict(lambda: defaultdict(lambda: defaultdict()))
    selenocysteine_id = 1

    with open(gtf_file) as FR:
        for line in FR:
            features = line.rsplit("\t")
            if features[2] == "Selenocysteine":
                chr = features[0]
                start = int(features[3])
                end = int(features[4])
                m = re.search('transcript_id "(\S+?)";', features[8])
                if m:
                    transcript_id = m.group(1)
                else:
                    print "Could not find transcript ID for selenocysteine "+chr+"_"+start
                strand=""
                if(features[6]=="+"):
                    strand="1"
                else:
                    strand="-1"
                secs[transcript_id][selenocysteine_id]["start"] = start
                secs[transcript_id][selenocysteine_id]["end"] = end
                secs[transcript_id][selenocysteine_id]["chr"] = chr
                secs[transcript_id][selenocysteine_id]["strand"] = strand

                selenocysteine_id += 1

    return secs


## Store selenocysteines in results DB
def store_secs_in_results_db(results_db, secs, price_tmp):
    table = "secs"

    # Connect to results db
    try:
        con = sqlite.connect(results_db)
    except:
        print "Could not connect to results DB!"
        sys.exit()

    with(con):
        cur = con.cursor()

        # Drop existing sec table
        query_drop = "DROP TABLE IF EXISTS '" + table + "';"
        cur.execute(query_drop)

        # Create table
        query_create = "CREATE TABLE IF NOT EXISTS '" + table + "' (" \
                        "'id' int(10) NOT NULL default '', " \
                        "'tr_stable_id' varchar(128) NOT NULL default '', " \
                        "'chr' char(5) NOT NULL default '', " \
                        "'strand' int(2) NOT NULL default '', " \
                        "'start' int(10) NOT NULL default '', " \
                        "'stop' int(10) NOT NULL default '');"
        cur.execute(query_create)

    # Create csv
    tmp_sec_file = price_tmp + "/secs_db.csv"
    FW = open(tmp_sec_file, 'w')
    for transcript_id in secs.keys():
        for sec_id in secs[transcript_id].keys():
            FW.write(str(sec_id) + "," + transcript_id + "," + str(secs[transcript_id][sec_id]["chr"]) + "," +
                     secs[transcript_id][sec_id]["strand"] + "," + str(
                secs[transcript_id][sec_id]["start"]) + "," + str(secs[transcript_id][sec_id]["end"]) + "\n")
    FW.close()

    # Dump
    try:
        os.system("sqlite3 -separator , " + results_db + " \".import " + tmp_sec_file + " " + table + "\"")
    except:
        print "Could not dump csv file " + tmp_sec_file + " into table " + table + " of DB " + results_db
        sys.exit()
    os.system("rm -rf " + tmp_sec_file)

    return

## determine annotation
def get_annotation(start, aTIS, stop, strand):

    annotation=""

    if strand == '1':
        if aTIS=="NA":
            annotation = "ntr"
        elif start == aTIS:
            annotation = "aTIS"
        elif start < aTIS:
            annotation = "5UTR"
        elif start > aTIS and start<=stop:
            annotation = "CDS"
        else:
            annotation = "3UTR"
    elif strand == '-1':
        if aTIS=="NA":
            annotation = "ntr"
        elif start == aTIS:
            annotation = "aTIS"
        elif start > aTIS:
            annotation = "5UTR"
        elif start < aTIS and start>=stop:
            annotation = "CDS"
        else:
            annotation = "3UTR"

    return annotation

## Calc distances to aTIS and transcript start in cDNA positions
def calc_dists(start_site_pos, aTIS, exons, strand):

    #Init
    dist_to_transcript_start = 0
    dist_to_aTIS = 0

    #Key is to determine the cDNA position of the transcript start (is of course position 1), of the aTIS and of the start sit
    tr_length=0
    cdnapos_aTIS=0
    cdnapos_start_site_pos=0
    atis_found = False
    start_site_found = False
    if aTIS=="NA":
        dist_to_aTIS="NA"
        atis_found = True
    if strand == "1":
        for rank in sorted(exons.iterkeys()):
            if atis_found==False:
                if aTIS<exons[rank]['end'] and aTIS>exons[rank]['start']:
                    cdnapos_aTIS = aTIS - exons[rank]['start'] + 1 + tr_length
                    atis_found = True
            if start_site_found==False:
                if start_site_pos<exons[rank]['end'] and start_site_pos>exons[rank]['start']:
                    cdnapos_start_site_pos = start_site_pos - exons[rank]['start'] + 1 + tr_length
                    start_site_found = True
            if atis_found==True and start_site_found==True:
                break
            tr_length = tr_length + exons[rank]['end'] - exons[rank]['start'] + 1
        dist_to_transcript_start = cdnapos_start_site_pos
        if dist_to_aTIS!="NA":
            dist_to_aTIS = cdnapos_start_site_pos - cdnapos_aTIS
    elif strand == "-1":
        for rank in sorted(exons.iterkeys()):
            if atis_found==False:
                if aTIS < exons[rank]['end'] and aTIS > exons[rank]['start']:
                    cdnapos_aTIS = exons[rank]['end'] - aTIS + 1 + tr_length
                    atis_found = True
            if start_site_found==False:
                if start_site_pos < exons[rank]['end'] and start_site_pos > exons[rank]['start']:
                    cdnapos_start_site_pos = exons[rank]['end'] - start_site_pos + 1 + tr_length
                    start_site_found = True
            if atis_found == True and start_site_found == True:
                break
            tr_length = tr_length + exons[rank]['end'] - exons[rank]['start'] + 1
        dist_to_transcript_start = cdnapos_start_site_pos
        if dist_to_aTIS!="NA":
            dist_to_aTIS = cdnapos_start_site_pos - cdnapos_aTIS

    return dist_to_transcript_start, dist_to_aTIS

## Get aTIS and transcript start position
def get_aTIS(transcript_id, strand, ens_db):

    #Init
    aTIS=0
    stop=0

    #Connection to ens_db
    try:
        con = sqlite.connect(ens_db)
    except:
        print "ERROR: could not connect to database "+ens_db
        sys.exit()

    with con:
        cur = con.cursor()

        if strand == "1":
            query = "SELECT translation.seq_start, exon.seq_region_start FROM translation JOIN exon ON " \
                    "translation.start_exon_id = exon.exon_id JOIN transcript ON translation.translation_id=" \
                    "transcript.canonical_translation_id WHERE transcript.transcript_id = '" + str(transcript_id) + "';"
            if cur.execute(query):
                out = cur.fetchone()
                if out is not None:
                    start_exon_seq_region_start = int(out[1])
                    seq_start = int(out[0])
                    aTIS = start_exon_seq_region_start + seq_start - 1
            query = "SELECT translation.seq_end, exon.seq_region_start FROM translation JOIN exon ON " \
                    "translation.end_exon_id = exon.exon_id JOIN transcript ON translation.translation_id = " \
                    "transcript.canonical_translation_id WHERE transcript.transcript_id = '" + str(transcript_id) + "';"
            if cur.execute(query):
                out = cur.fetchone()
                if out is not None:
                    end_exon_seq_region_start = int(out[1])
                    seq_end = int(out[0])
                    stop = end_exon_seq_region_start + seq_end - 1
        elif strand == "-1":
            query = "SELECT translation.seq_start, exon.seq_region_end FROM translation JOIN exon ON " \
                    "translation.start_exon_id = exon.exon_id JOIN transcript ON translation.translation_id=" \
                    "transcript.canonical_translation_id WHERE transcript.transcript_id = '" + str(transcript_id) + "';"
            if cur.execute(query):
                out = cur.fetchone()
                if out is not None:
                    start_exon_seq_region_end = int(out[1])
                    seq_start = int(out[0])
                    aTIS = start_exon_seq_region_end - seq_start + 1
            query = "SELECT translation.seq_end, exon.seq_region_end FROM translation JOIN exon ON " \
                    "translation.end_exon_id = exon.exon_id JOIN transcript ON translation.translation_id = " \
                    "transcript.canonical_translation_id WHERE transcript.transcript_id = '" +str(transcript_id)+ "';"
            if cur.execute(query):
                out = cur.fetchone()
                if out is not None:
                    end_exon_seq_region_end = int(out[1])
                    seq_end = int(out[0])
                    stop = end_exon_seq_region_end - seq_end + 1

    #If no translation info
    if aTIS == 0:
        aTIS = "NA"
    if stop == 0:
        stop = "NA"

    return aTIS, stop

## Get start and stop coordinates, also from the exons
def get_start_and_stop_coordinates(orfs, price_id, strand):

    #Init
    start=0
    stop=0
    starts_list=""
    ends_list=""
    coding_parts = defaultdict(lambda: defaultdict())
    exon_number=1

    if strand == "1":
        for exon in sorted(orfs[price_id]['exons'].keys()):
            if exon==1:
                start = orfs[price_id]['exons'][exon]['start']
            if exon==len(orfs[price_id]['exons'].keys()):
                stop = orfs[price_id]['exons'][exon]['end']
            starts_list = starts_list+str(orfs[price_id]['exons'][exon]['start'])+"_"
            ends_list = ends_list+str(orfs[price_id]['exons'][exon]['end'])+"_"
            coding_parts[exon_number]['start'] = orfs[price_id]['exons'][exon]['start']
            coding_parts[exon_number]['end'] = orfs[price_id]['exons'][exon]['end']
            exon_number+=1
    elif strand=="-1":
        for exon in sorted(orfs[price_id]['exons'].keys(), reverse=True):
            if exon==1:
                start = orfs[price_id]['exons'][exon]['start']
            if exon==len(orfs[price_id]['exons'].keys()):
                stop = orfs[price_id]['exons'][exon]['end']
            starts_list = starts_list+str(orfs[price_id]['exons'][exon]['end'])+"_"
            ends_list = ends_list+str(orfs[price_id]['exons'][exon]['start'])+"_"
            coding_parts[exon_number]['start'] = orfs[price_id]['exons'][exon]['end'] #Make strand undependent
            coding_parts[exon_number]['end'] = orfs[price_id]['exons'][exon]['start']
            exon_number+=1

    #Cut of trailing underscores in lists
    starts_list = starts_list.rstrip('_')
    ends_list = ends_list.rstrip('_')

    return (start, stop, starts_list, ends_list, coding_parts)

# Get ensembl transcript ID based on stable id
def get_transcript_id_and_exons(stable_id, ens_db):

    #Init
    transcript_id = 0
    exons = defaultdict(lambda: defaultdict())

    #Connection
    try:
        con = sqlite.connect(ens_db)
    except:
        print "ERROR: could not connect to "+ens_db
        sys.exit()

    with con:
        cur = con.cursor()
        #Transcript ID
        query = "SELECT transcript_id FROM transcript WHERE stable_id='"+stable_id+"';"
        cur.execute(query)
        transcript_id = cur.fetchone()[0]
        #exons
        query = "SELECT exon_transcript.rank, exon.seq_region_start, exon.seq_region_end FROM exon_transcript JOIN exon" \
                " ON exon_transcript.exon_id=exon.exon_id WHERE exon_transcript.transcript_id='"+str(transcript_id)+"' ORDER " \
                "BY exon_transcript.rank;"
        cur.execute(query)
        output = cur.fetchall()
        for i in range(0, len(output)):
            rank = output[i][0]
            exons[rank]['start'] = output[i][1]
            exons[rank]['end'] = output[i][2]

    return transcript_id, exons

## Get orfs from PRICE output
def parse_price_output(tsv, loc):

    #Read in tsv and loc files
    all_orfs = read_tsv(tsv)
    filtered_orfs = read_loc(loc)

    #Combine data of both files to obtain filtered orfs from price
    price_orfs = combine_data(all_orfs, filtered_orfs)

    return price_orfs

## Combine orfs to get filtered orfs from PRICE
def combine_data(all, filtered):

    #Init
    orfs = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict())))

    #Only filtered orfs to keep
    for orf in filtered:
        id = all[orf]['id']
        orfs[id]['start_codon'] = all[orf]['start_codon']
        orfs[id]['price_annotation'] = all[orf]['price_annotation']
        orfs[id]['pval'] = filtered[orf]['pval']
        orfs[id]['tr_stable_id'] = id.split('_')[0]
        #Parse location string
        m = re.search('^(.+)([\+\-]):(.+)',all[orf]['location'])
        if m:
            orfs[id]['chr'] = str(m.group(1))
            if m.group(2)=='+':
                orfs[id]['strand'] = '1'
                #Parse start and end of all exons
                exon_str = m.group(3)
                exons = exon_str.split('|')
                exon_number = 1
                for exon in exons:
                    coor = exon.split('-')
                    orfs[id]['exons'][exon_number]['start'] = int(coor[0])+1 # Switch to 1-based coordinates (Ensembl)
                    orfs[id]['exons'][exon_number]['end'] = int(coor[1])
                    exon_number += 1
            else:
                orfs[id]['strand'] = '-1'
                # Parse start and end of all exons
                exon_str = m.group(3)
                exons = exon_str.split('|')
                exon_number = len(exons)
                for exon in exons:
                    coor = exon.split('-')
                    orfs[id]['exons'][exon_number]['start'] = int(coor[1])
                    orfs[id]['exons'][exon_number]['end'] = int(coor[0]) +1 #Switch to 1 based coordinates
                    exon_number = exon_number-1

    return orfs

## Read filtered orfs from loc file
def read_loc(loc):

    #Init
    filt_orfs = defaultdict(lambda: defaultdict())

    with open(loc, 'r') as FR:
        for line in FR:
            #Chomp new line
            line.rstrip('\n')
            #Split features
            prior_location, other_features = line.split('\t')
            other_features = other_features.split(' ')
            #Parse orfs
            m = re.search('^p=(.+)$', other_features[3])
            if m:
                filt_orfs[prior_location]['pval'] = m.group(1)

    return filt_orfs

## Read all orf:q
# s from tsv file
def read_tsv(tsv):

    all_orfs = defaultdict(lambda: defaultdict())

    with open(tsv, 'r') as FR:
        #Skip header
        next(FR)
        for line in FR:
            #Chomp newline and split by tab
            line.rstrip('\n')
            orf_features = line.split('\t')
            #Parse, store by prior_location (location before start codon assignment) as this is the key we will find in the location file (.loc)
            prior_location = orf_features[3]
            all_orfs[prior_location]['location'] = orf_features[2]
            all_orfs[prior_location]['id'] = orf_features[1]
            all_orfs[prior_location]['start_codon'] = orf_features[4]
            all_orfs[prior_location]['price_annotation'] = orf_features[5]

    return all_orfs

## Convert CIT file to LOC file
def convert_cit(cit, loc, price_operatable):

    command = price_operatable+" -e ViewCIT -o "+loc+" -m Location "+cit
    print command
    os.system(command)

    return

## Run PRICE
def run_price(price_operatable, price_prefix, gref, readtype, bam_untr, bam_tr, cores, fdr, price_tmp, max_ram):

    #make bamlist file if treated condition is available
    reads = ""
    if readtype=="ribo":
        reads = price_tmp+"/reads.bamlist"
        with open(reads, 'w') as FW:
            FW.write(bam_untr+"\n")
            FW.write(bam_tr+"\n")
    else:
        reads = bam_untr

    #If max RAM is not default, change gedi bash wrapper
    if max_ram!=0:
        min_ram = int(max_ram/4)
        operatable_copy = "copy_operatable_price"
        with open(price_operatable, 'r') as FR:
            with open(operatable_copy, 'w') as FW:
                for line in FR:
                    m = re.search('^(java.+-Xmx)\d+?m -Xms\d+?m(.+)$', line)
                    if m:
                        output_line = m.group(1)+str(max_ram)+"g -Xms"+str(min_ram)+"g"+m.group(2)
                        FW.write(output_line)
                    else:
                        FW.write(line)
        os.system("mv "+operatable_copy+" "+price_operatable)
        os.system("chmod 775 "+price_operatable)

    #Execute
    command = price_operatable+" -Dsmile.threads=1 -e Price -reads "+reads+" -genomic "+gref+" -prefix "+price_prefix+\
                " -nthreads "+str(cores)+" -fdr "+str(fdr)
    if readtype=="ribo":
        command = command+" -percond"
    print command
    print
    os.system(command)

    #Delete bamlist
    if readtype=="ribo":
        os.system("rm -rf "+reads)

    return

## Prepare price genome reference
def prepare_genome_price(price_operatable, fasta, gtf, species, ens_v):

    #Execute
    gref_name = species+"_"+ens_v
    command = price_operatable+" -e IndexGenome -s "+fasta+" -a "+gtf+" -n "+gref_name+" -nobowtie -nostar -nokallisto"
    print command
    os.system(command)

    return gref_name

## Install PRICE
def install_price(price_tmp):

    #Download
    try:
        os.system("wget --no-check-certificate --quiet https://github.com/erhard-lab/gedi/releases/download/Price_1.0.2/Price_1.0.2.tar.gz")
    except:
        print "Could not download PRICE 1.0.2"
        sys.exit()

    #Move to tmp folder
    os.system("mv Price_1.0.2.tar.gz "+price_tmp)

    #Unpack
    try:
        os.system("gunzip --quiet "+price_tmp+"/Price_1.0.2.tar.gz")
    except:
        print "Could not unpack with gunzip"
        sys.exit()
    try:
        os.system("tar -xf "+price_tmp+"/Price_1.0.2.tar -C "+price_tmp)
    except:
        print "Could not unpack with tar"
        sys.exit()
    #Remove tar file
    os.system("rm -rf "+price_tmp+"/Price_1.0.2.tar")

    #Link to operatable
    price_operatable = price_tmp+"/Price_1.0.2/gedi"
    price_folder = price_tmp+"/Price_1.0.2"

    return price_operatable, price_folder

## Get chromosomes ##
def get_chrs(chrSizeFile, species):
    #Init
    chrs = {}

    #open file
    try:
        FR = open(chrSizeFile, 'r')
    except:
        print "ERROR with opening chromosome sizes file"
    #parse
    for line in FR:
        parts = re.split('\W+',line)
        if species=="fruitfly" and parts[0]=="M":
            chrs["dmel_mitochondrion_genome"] = parts[1]
        else:
            chrs[parts[0]] = parts[1]

    return chrs

## Get chromosomes in XML ##
def get_chrsXml(chrSizeFile, species):
    #Init
    chrs = {}

    #open file
    try:
        FR = open(chrSizeFile, 'r')
    except:
        print "ERROR with opening chromosome sizes xml file"
    #parse
    for line in FR:
        pattern = re.compile('contigName=\"(\w{1,2})\" totalBases\W\"(\d+)\"')
        m = pattern.search(line)
        if m:
            if species=="fruitfly" and m.group(1)=="M":
                chrs["dmel_mitochondrion_genome"] = m.group(2)
            else:
                chrs[m.group(1)] = m.group(2)

    return chrs

## Create BIN chromosomes
def create_BIN_chrom(cores,BIN_chrom_dir,chrs,igenomes_root,speciesLatin,assembly):

    #create dir
    os.mkdir(BIN_chrom_dir)

    #Split over multiple processes
    pool = Pool(processes=int(cores))
    [pool.apply_async(create_BIN_chrom_chr, args=(chr, BIN_chrom_dir, igenomes_root, speciesLatin, assembly)) for chr in chrs.keys()]
    pool.close()
    pool.join()

    return

def create_BIN_chrom_chr(chr, BIN_chrom_dir, igenomes_root, speciesLatin, assembly):

    chr_file = igenomes_root+"/"+speciesLatin+"/Ensembl/"+assembly+"/Sequence/Chromosomes/"+chr+".fa"
    BIN_chr_file = BIN_chrom_dir+"/"+chr+".fa"

    with open(chr_file, 'r') as FR:
        with open(BIN_chr_file, 'w') as FW:
            for line in FR:
                if not line.startswith('>'):
                    line.rstrip('\n')
                    FW.write(line)

    return

##Get total mapped reads
def get_seq_depth(result_db, sample_name):

    #Init
    seq_depth=""

    #Make database connection
    try:
        con = sqlite.connect(result_db)
    except:
        print "ERROR: could not connect to "+result_db
        sys.exit()

    with con:
        cur = con.cursor()

        query = "SELECT mapped_T FROM statistics WHERE sample='"+sample_name+"' AND type='genomic';";
        if cur.execute(query):
            seq_depth = float(cur.fetchone()[0])

    return seq_depth

## Get arguments from arguments table
def get_arguments(db):

    #Make result db connection
    try:
        con = sqlite.connect(db)
    except:
        print "ERROR: could not connect to "+db
        sys.exit()

    with con:
        cur = con.cursor()

        query = "SELECT value FROM arguments WHERE variable='igenomes_root';"
        if cur.execute(query):
            igenomes_root = str(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the igenomes_root argument in "+db
            print query
            sys.exit()

        query = "SELECT value FROM arguments WHERE variable='species';"
        if cur.execute(query):
            species = str(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the species argument in "+db
            print query
            sys.exit()

        query = "SELECT value FROM arguments WHERE variable='ensembl_version';"
        if cur.execute(query):
            version = int(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the ensembl version argument in "+db
            print query
            sys.exit()

        query = "SELECT value FROM arguments WHERE variable='ens_db';"
        if cur.execute(query):
            ens_db = str(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the ensembl database argument in "+db
            print query
            sys.exit()

        query = "SELECT value FROM arguments WHERE variable='readtype';"
        if cur.execute(query):
            readtype = str(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the readtype argument in " + db
            print query
            sys.exit()

        query = "SELECT value FROM arguments WHERE variable='price_bam_untr';"
        if cur.execute(query):
            bam_untr = str(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the untreated BAM argument in " + db
            print query
            sys.exit()

        bam_tr = ""
        if readtype=="ribo":
            query = "SELECT value FROM arguments WHERE variable='price_bam_tr';"
            if cur.execute(query):
                bam_tr = str(cur.fetchone()[0])
            else:
                print "ERROR: could not fetch the treated BAM argument in " + db
                print query
                sys.exit()
        elif (readtype=="ribo_untr"):
            bam_tr = "Not available"

        query = "SELECT value FROM arguments WHERE variable='nr_of_cores';"
        if cur.execute(query):
            cores = str(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the number of cores argument in " + db
            print query
            sys.exit()

        query = "SELECT value FROM arguments WHERE variable='tr_calling';"
        if cur.execute(query):
            tr_calling_method = str(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the tr_calling argument in "+db
            print query
            sys.exit()

        query = "SELECT value FROM arguments WHERE variable='run_name';"
        if cur.execute(query):
            run_name = str(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the run_name argument in "+db
            print query
            sys.exit()

        query = "SELECT value FROM arguments WHERE variable='mapper';"
        if cur.execute(query):
            mapper = str(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the mapper argument in "+db
            print query
            sys.exit()

        query = "SELECT value FROM arguments WHERE variable='unique';"
        if cur.execute(query):
            uniq = str(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the unique argument in "+db
            print query
            sys.exit()

    return igenomes_root, species, version, ens_db, readtype, bam_untr, bam_tr, cores, tr_calling_method, run_name, mapper, uniq

#Index bams if needed
def index_bams(bam_untr, bam_tr):

    for bam in [bam_untr, bam_tr]:
        #Check if indexed file already exists
        if not os.path.isfile(bam+".bai"):
            #If not: index bam file with samtools
            os.system("samtools index "+bam)

    return

def convert_time(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)

## Data Dumper for recursively printing nested dictionaries and defaultDicts ##
def print_dict(dictionary, indent='', braces=0):
    """
    :param dictionary: dict or defaultdict to be printed
    :param indent: indentation
    :param braces:
    :return: void
    """

    for key, value in dictionary.iteritems():
        if isinstance(value, dict):
            print '%s%s%s%s' % (indent, braces * '[', key, braces * ']')
            print_dict(value, indent + '  ', braces)
        else:
            print indent + '%s = %s' % (key, value)




###### DIRECT TO MAIN ############
if __name__ == "__main__":
    try:
        main()
    except Exception, e:
        traceback.print_exc()
##################################