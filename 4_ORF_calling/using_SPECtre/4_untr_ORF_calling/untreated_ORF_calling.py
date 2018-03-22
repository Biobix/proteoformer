#!/usr/bin/env python
import traceback
import sys
import os
import sqlite3 as sqlite
import time
import argparse
from collections import defaultdict
from multiprocessing import Pool
import math
import re


__author__ = 'Steven Verbruggen'


'''
UNTREATED ORF CALLING

usage: untreated_ORF_calling.py --result_db [DB] [--help] [--workdir [FOLDER]]
                                [--tmp [FOLDER]] [--untr_bam [BAM]]
                                [--offsets [LIST]] [--cores [INTEGER]]
                                [--threads_per_chrom [INTEGER]]

This tool is part of the PROTEOFORMER pipeline. It calls translatedORFs
without the use of LTM/HARR treated samples and TIS calling.

Mandatory parameters:
  --result_db [DB], -r [DB]
                        The SQLite results database

Optional parameters:
  --help, -h            Show help message and exit
  --workdir [FOLDER], -w [FOLDER]
                        Working directory (default: CWD)
  --tmp [FOLDER], -t [FOLDER]
                        Temporary folder (default: workdir/tmp)
  --untr_bam [BAM], -b [BAM]
                        BAM file of untreated sample (default: get BAM file
                        path out of results database)
  --offsets [LIST], -o [LIST]
                        custom list of offsets, but defaultly the program
                        takes the untreated P-site offsets already stored in
                        the database during mapping
  --cores [INTEGER], -c [INTEGER]
                        Defaultly, the program takes the amount of cores
                        suggested during mapping, but you can suggest a
                        different number of cores for this part here
  --threads_per_chrom [INTEGER], -x [INTEGER]
                        The number of threads used per chromosome to run
                        chromosomal SPECTRE runs with
'''

def main():
    starttime = time.time()

    print
    print "#########################"
    print "# Untreated ORF Calling #"
    print "#########################"
    print

    #Parse arguments from command line
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS,
                                     description="This tool is part of the PROTEOFORMER pipeline. It calls translated"
                                                 "ORFs without the use of LTM/HARR treated samples and TIS calling.",
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
    opt_args.add_argument("--untr_bam", "-b", action="store", required=False, nargs="?", metavar="BAM", default="",
                          type=str, help="BAM file of untreated sample (default: get BAM file path out of results "
                                         "database)")
    opt_args.add_argument("--offsets", "-o", action="store", required=False, nargs="?", metavar="LIST", default="",
                          type=str, help="custom list of offsets, but defaultly the program takes the untreated P-site "
                                         "offsets already stored in the database during mapping. Caution: SPECtre is not "
                                         "yet built for handling bigger offsets list like in PROTEOFOMER -> Use most-"
                                         "abundant RPF lengths until further notice.")
    opt_args.add_argument("--cores", "-c", action="store", required=False, nargs="?", metavar="INTEGER", default="",
                          type=int, help="Defaultly, the program takes the amount of cores suggested during mapping,"
                                         " but you can suggest a different number of cores for this part here")
    opt_args.add_argument("--threads_per_chrom", "-x", action="store", required=False, nargs="?", metavar="INTEGER",
                          default=1, type=int, help="The number of threads used per chromosome to run chromosomal "
                                                    "SPECTRE runs with")

    args = parser.parse_args()

    #default workdir is CWD, default tmp, default BAM, default offsets, default cores
    if args.workdir == "":
        args.workdir = os.getcwd()
    if args.tmp == "":
        args.tmp = args.workdir+"/tmp/"
    if not os.path.isdir(args.tmp):
        os.mkdir(args.tmp)
    spectre_tmp = args.tmp+"SPECtre/"
    if not os.path.isdir(spectre_tmp):
        os.mkdir(spectre_tmp)
    if args.untr_bam == "":
        args.untr_bam = get_bam_from_arguments(args.result_db)
    if args.offsets == "":
        args.offsets = get_offsets_as_list(args.result_db)
    if args.cores == "":
        args.cores = get_cores(args.result_db)

    #Get parameters from arguments table
    table_args = argparse.Namespace()
    table_args.tr_calling_method, table_args.igenomes_root, table_args.species, table_args.ens_v, table_args.ensdb = \
        get_arguments(args.result_db)

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

    # GTF file
    gtf_file = table_args.igenomes_root + "/" + speciesLatin + "/Ensembl/" + assembly + "/Annotation/Genes/genes_" + str(
        table_args.ens_v) + ".gtf"

    #List parameters
    print "Parameters:"
    for arg in vars(args):
        print '    %-15s\t%s' % (arg, getattr(args, arg))
    print "Parameters from arguments table:"
    for arg in vars(table_args):
        print '    %-15s\t%s' % (arg, getattr(table_args, arg))
    print
    sys.stdout.flush()

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
        create_BIN_chrom(args.cores,BIN_chrom_dir,chrs,table_args.igenomes_root,speciesLatin,assembly)

    ##########
    #  MAIN  #
    ##########

    #Generate isoform file of translated transcript table
    print "[%s]: Generate isoform file from translated transcript table" % (convert_time(time.time()-starttime))
    sys.stdout.flush()
    isoform_file = generate_iso_file(args.result_db, spectre_tmp, table_args.tr_calling_method)
    print "[%s]: COMPLETE: Generate isoform file from translated transcript table" % (convert_time(time.time()-starttime))
    print
    sys.stdout.flush()


    #Download SPECtre installation
    downloadlink_spectre = "https://github.com/mills-lab/spectre/archive/master.zip"
    installation_dir = download_spectre(spectre_tmp, downloadlink_spectre)

    #Run SPECtre
    print
    print "[%s]: Run SPECtre" % (convert_time(time.time()-starttime))
    sys.stdout.flush()
    merged_results_file = run_spectre(installation_dir, args.untr_bam, spectre_tmp, isoform_file, gtf_file, args.offsets, chrs,
                args.cores, args.threads_per_chrom, starttime)
    print "[%s]: COMPLETE: SPECtre run" % (convert_time(time.time() - starttime))
    print
    sys.stdout.flush()

    #Remove SPECtre installation files
    #os.system("rm -rf "+installation_dir)

    #Parse spectre results into assembly table
    print
    print "[%s]: Parse SPECtre results" % (convert_time(time.time()-starttime))
    sys.stdout.flush()
    parse_results(args.result_db, spectre_tmp, table_args.ensdb, BIN_chrom_dir, merged_results_file, table_args.tr_calling_method)
    print "[%s]: COMPLETE: SPECtre results parsing" % (convert_time(time.time() - starttime))
    print
    sys.stdout.flush()

    #Remove tmp spectre files
    #os.system("rm -rf "+spectre_tmp)

    #End of program message
    print
    print "[%s]: PROGRAM COMPLETE" % (convert_time(time.time()-starttime))
    print
    sys.stdout.flush()

    return

## Parse SPECtre results into assembly table
def parse_results(results_db, spectre_tmp, ens_db, BIN_chrom_dir, merged_results_file, tr_calling_method):

    #Load in log file for statistics
    threshold_translation = load_log_file(spectre_tmp)

    #Load in SPECtre results file
    spectre_results = load_spectre_results(spectre_tmp, merged_results_file)

    #Get all needed features
    tr_features = get_assembly_features(spectre_results, ens_db, results_db, BIN_chrom_dir, threshold_translation)

    print
    print "SPECtre identified "+str(len(spectre_results))+" coding possibilities"
    print str(len(tr_features))+" passed the translation threshold"
    print


    #Drop possible existing table and initiate new table
    tis_id = get_id(results_db, tr_calling_method)
    table_name = "TIS_"+str(tis_id)+"_transcripts"
    create_new_table(results_db, table_name)

    # Construct csv
    csv_name = spectre_tmp+"/assembly_table.csv"
    filtered_out_csv = spectre_tmp+"/filtered_out.csv"
    count_orfs = construct_csv(csv_name, filtered_out_csv, tr_features)
    print str(count_orfs)+" to be saved in the new assembly table"

    # dump into table
    dump_into_sqlite(results_db, table_name, csv_name)

    # Remove csv
    #os.system("rm -rf "+csv_name)

    return

## Load log file
def load_log_file(spectre_tmp):

    #Init
    translation_threshold = -1

    #Def file path
    log_file = spectre_tmp+"/spectre_merge.log"

    #Read and search needed statistics
    with open(log_file, 'r') as FR:
        pre_metric_line = FR.readline()
        while pre_metric_line:
            m = re.search('##########', pre_metric_line)
            if m:
                #Read six lines with metrics
                for i in range(0,6):
                    pre_metric_line = FR.readline()
                    m = re.search('__main__ - INFO - (.+)$', pre_metric_line)
                    if m:
                        metric_line = m.group(1).rstrip('\n')
                        print metric_line
                        m = re.search('Translational Threshold = (\d\.\d+)', metric_line)
                        if m:
                            translation_threshold = float(m.group(1))
                    else:
                        print "Warning: no match for metrics in: "+pre_metric_line
            pre_metric_line = FR.readline()
    print

    return translation_threshold

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
        for spectre_id in tr_features:
            if tr_features[spectre_id]['start_codon']!="Unvalid" and tr_features[spectre_id]['aa_seq'].split(',')[0]!="early_stop" and \
                    tr_features[spectre_id]['aa_seq'] != "no_stop":
                line = tr_features[spectre_id]['tr_stable_id']+","+\
                    tr_features[spectre_id]['chr']+","+\
                    str(tr_features[spectre_id]['strand'])+","+\
                    str(tr_features[spectre_id]['start'])+","+\
                    tr_features[spectre_id]['start_codon']+","+\
                    str(tr_features[spectre_id]['stop'])+","+\
                    tr_features[spectre_id]['starts_list']+","+\
                    tr_features[spectre_id]['ends_list']+","+\
                    str(tr_features[spectre_id]['dist_to_transcript_start'])+","+\
                    str(tr_features[spectre_id]['dist_to_aTIS'])+","+\
                    tr_features[spectre_id]['annotation']+","+\
                    tr_features[spectre_id]['biotype']+","+\
                    tr_features[spectre_id]['aTIS_call']+","+\
                    str(tr_features[spectre_id]['peak_shift'])+","+\
                    str(tr_features[spectre_id]['count'])+","+\
                    str(tr_features[spectre_id]['Rltm_min_Rchx'])+","+\
                    str(tr_features[spectre_id]['coverage'])+","+\
                    str(tr_features[spectre_id]['FPKM'])+","+\
                    tr_features[spectre_id]['SNP']+","+\
                    tr_features[spectre_id]['INDEL']+","+\
                    tr_features[spectre_id]['tr_seq']+","+\
                    tr_features[spectre_id]['aa_seq']+","+\
                    str(spectre_id)+"\n"
                counter+=1
                FW.write(line)
            else:
                line = tr_features[spectre_id]['tr_stable_id']+","+ \
                    tr_features[spectre_id]['chr'] + "," + \
                    str(tr_features[spectre_id]['strand']) + "," + \
                    str(tr_features[spectre_id]['start']) + "," + \
                    tr_features[spectre_id]['start_codon'] + "," + \
                    str(tr_features[spectre_id]['stop']) + "," + \
                    tr_features[spectre_id]['starts_list'] + "," + \
                    tr_features[spectre_id]['ends_list'] + "," + \
                    str(tr_features[spectre_id]['dist_to_transcript_start']) + "," + \
                    str(tr_features[spectre_id]['dist_to_aTIS']) + "," + \
                    tr_features[spectre_id]['tr_seq'] + "," + \
                    tr_features[spectre_id]['aa_seq'] + "," + \
                    str(spectre_id)+"\n"
                FiltFW.write(line)
        FiltFW.close()

    return counter

### Create new assembly table
def create_new_table(results_db, table_name):

    #Connect to sqlite db
    con = sqlite.connect(results_db)

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
                        "'tr_seq' TEXT NOT NULL default ''," \
                        "'aa_seq' TEXT NOT NULL default '',"\
                        "'spectre_id' int(10) NOT NULL default '');"
        cur.execute(create_query)

    return

## Get TIS ID
def get_id(result_db, tr_calling_method):

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
                        "'SNP' varchar(20) default '', " \
                        "'indel' varchar(20) default '', " \
                        "'filter' varchar(20) default '', " \
                        "'tr_calling' varchar(20) default '', " \
                        "'TIS_calling' varchar(20) default '');"
        cur.execute(query_create)

        #Add parameters and get ID
        query_add = "INSERT INTO 'TIS_overview'(filter,tr_calling,TIS_calling,SNP,indel) "\
                    "VALUES('none','"+tr_calling_method+"','SPECtre','NO','NO');"
        cur.execute(query_add)
        TIS_id = int(cur.lastrowid)

    return TIS_id

### Parse results from SPECtre results to assembly table features
def get_assembly_features(spectre_results, ens_db, results_db, BIN_chrom_dir, threshold_translation):

    #Init
    tr_features = defaultdict(lambda: defaultdict())

    reads_for, reads_rev = get_reads(results_db)

    for spectre_id in spectre_results.keys():
    #for spectre_id in ['4172']:
        if spectre_results[spectre_id]['SPEC_score_posterior_CDS']=="NA":
            #SPECtre fails to calculate SPEC score for some ID's. Skip these.
            continue
        if float(spectre_results[spectre_id]['SPEC_score_posterior_CDS'])>threshold_translation:
            #Stable ID
            tr_features[spectre_id]['tr_stable_id'] = spectre_results[spectre_id]['transcript_id']
            #Ensembl transcript ID
            transcript_id, exons = get_transcript_id_and_exons(spectre_results[spectre_id]['transcript_id'], ens_db)

            #Chromosome
            tr_features[spectre_id]['chr'] = spectre_results[spectre_id]['chr']
            #Strand
            if spectre_results[spectre_id]['strand']=="+":
                tr_features[spectre_id]['strand'] = "1"
            elif spectre_results[spectre_id]['strand']=="-":
                tr_features[spectre_id]['strand'] = "-1"
            else:
                print "ERROR: strand unrecognised for SPECtre id "+spectre_id+" (transcript ID "+spectre_results[spectre_id]['transcript_id']+")"
            #Coordinates
            tr_features[spectre_id]['start'], tr_features[spectre_id]['stop'], tr_features[spectre_id]['starts_list'], \
                tr_features[spectre_id]['ends_list'], coding_parts = convert_coordinates(spectre_results, spectre_id)
            #Distance to transcript start and aTIS
            aTIS, stop = get_aTIS(transcript_id, tr_features[spectre_id]['strand'], ens_db)
            tr_features[spectre_id]['dist_to_transcript_start'], tr_features[spectre_id]['dist_to_aTIS'] = \
                calc_dists(tr_features[spectre_id]['start'], aTIS, exons, tr_features[spectre_id]['strand'])
            #Start codon, tr_seq, aa_seq
            tr_features[spectre_id]['start_codon'], tr_features[spectre_id]['tr_seq'], tr_features[spectre_id]['aa_seq'] = \
                construct_seq(coding_parts, tr_features[spectre_id]['strand'], BIN_chrom_dir, tr_features[spectre_id]['chr'], \
                                    tr_features[spectre_id]['tr_stable_id'], spectre_id)

            # print
            # print_dict(tr_features)
            # print
            #
            # print "exons:"
            # print_dict(coding_parts)
            # print "tr seq: "+tr_features[spectre_id]['tr_seq']
            # print "aa seq: "+tr_features[spectre_id]['aa_seq']
            # sys.stdout.flush()
            # sys.exit()


            #Check if unvalid start codon, no stop or early stop
            if tr_features[spectre_id]['start_codon'] == "Unvalid" or tr_features[spectre_id]["aa_seq"].split(',')[0]=="early_stop" \
                    or tr_features[spectre_id]["aa_seq"] == "no_stop":
                continue

            #Annotation based on aTIS and translation stop
            tr_features[spectre_id]['annotation'] = get_annotation(tr_features[spectre_id]['start'], aTIS, stop, tr_features[spectre_id]['strand'])
            tr_features[spectre_id]['biotype'] = spectre_results[spectre_id]['gene_type']
            #aTIS call and peak shift
            tr_features[spectre_id]['aTIS_call'] = "NA"
            tr_features[spectre_id]['peak_shift'] = "NA"
            #Count and coverage
            if tr_features[spectre_id]['strand'] == "1":
                tr_features[spectre_id]['count'], tr_features[spectre_id]['coverage'] = count_and_coverage(reads_for, coding_parts, tr_features[spectre_id]['chr'])
            else:
                tr_features[spectre_id]['count'], tr_features[spectre_id]['coverage'] = count_and_coverage(reads_rev, coding_parts,  tr_features[spectre_id]['chr'])
            #R is impossible to calculate without LTM
            tr_features[spectre_id]['Rltm_min_Rchx'] = "NA"
            #FPKM
            tr_features[spectre_id]['FPKM'] = float(spectre_results[spectre_id]['ribo_fpkm'])
            #SNP and indels: for the moment empty
            tr_features[spectre_id]['SNP'] = ""
            tr_features[spectre_id]['INDEL'] = ""

    return tr_features

## Get sequece of start codon, ORF cDNA
def construct_seq(coding_parts, strand, BIN_chrom_dir, chr, tr_stable_id, spectre_id):

    #Init
    tr_seq = ""
    aa_seq = ""
    #Go over all coding parts, situated over all exons
    for exon_part in coding_parts:
        tr_seq += get_seq(coding_parts[exon_part]['start'], coding_parts[exon_part]['stop'], chr, BIN_chrom_dir)
    #For reverse strand: take reverse complement in the end
    if strand == "-1":
        tr_seq = reverse_complement(tr_seq)

    #Get start codon
    start_codon = tr_seq[:3]
    # Check start codon for (near-)cognate
    if not re.search('[ACTG]TG|A[ACTG]G|AT[ACTG]', start_codon):
        start_codon = "Unvalid"
    else:
        #Get amino acid sequence
        aa_seq = translate_seq(tr_seq, tr_stable_id, spectre_id)

    return start_codon, tr_seq, aa_seq

## Translate sequence to amino acids
def translate_seq(tr_seq, tr_stable_id, spectre_id):

    #Init
    aa_seq = ""

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
        aa_seq+=codonTable[codon]
        #Control if stop codon is only at the end
        if codonTable[codon]=='*':
            if pos!=(len(tr_seq)-3):
                #print "Error: stop earlier than expected (transcript ID: "+tr_stable_id+", SPECtre ID: "+str(spectre_id)+")"
                aa_seq = "early_stop,"+aa_seq
                #print "tr seq: "+tr_seq
                #print "AA seq: "+aa_seq
                return aa_seq
    #Check if end codon is a STOP codon
    if codon!='TAA' and codon!='TAG' and codon!='TGA':
        #print "Error: no stop codon (TAA, TAG or TGA) found at the end (transcript ID: "+tr_stable_id+", SPECtre ID: "+str(spectre_id)+")"
        #print "tr seq: " + tr_seq
        #print "AA seq: " + aa_seq
        aa_seq="no_stop"
        return aa_seq
    #Check for near-cognate starts and replace the first amino acid for methionine then
    if re.search('[ACTG]TG|A[ACTG]G|AT[ACTG]', tr_seq[:3]):
        aa_seq = "M" + aa_seq[1:]
    #print "tr seq: "+tr_seq
    #print "aa seq: "+aa_seq

    return aa_seq

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

## Determine count and coverage
def count_and_coverage(reads, coding_parts, chr):

    #Init
    covered = 0
    not_covered = 0
    coverage = 0
    count = 0

    #Go over exon positions and find counts and coverage
    for exon_part in range(0, len(coding_parts)):
        for pos in range(coding_parts[exon_part]['start'], coding_parts[exon_part]['stop']+1):
            if pos in reads[chr]:
                covered+=1
                count = count + reads[chr][pos]
            else:
                not_covered+=1

    if covered!=0 or not_covered!=0:
        coverage = float(covered) / (float(covered) + float(not_covered))

    return count, coverage

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

## Convert coordinates from SPECtre to what assembly table needs
def convert_coordinates(spectre_results, spectre_id):

    #Init
    start = 0
    stop = 0
    starts_list = ""
    ends_list = ""

    if spectre_results[spectre_id]['strand'] == "+":
        coding_parts = get_exons_from_spectre(spectre_results, spectre_id)
        for i in range(0, len(coding_parts)):
            if i==0 and i==len(coding_parts)-1:
                start = coding_parts[i]['start']
                starts_list = str(coding_parts[i]['start'])
                stop = coding_parts[len(coding_parts) - 1]['stop']
                ends_list = str(coding_parts[i]['stop'])
            elif i == 0:
                start = coding_parts[i]['start']
                starts_list = str(coding_parts[i]['start'])
                ends_list = str(coding_parts[i]['stop'])
            elif i==len(coding_parts)-1:
                stop = coding_parts[len(coding_parts)-1]['stop']
                starts_list = starts_list+"_"+str(coding_parts[i]['start'])
                ends_list = ends_list+"_"+str(coding_parts[i]['stop'])
            else:
                starts_list = starts_list+"_"+str(coding_parts[i]['start'])
                ends_list = ends_list+"_"+str(coding_parts[i]['stop'])

    elif spectre_results[spectre_id]['strand'] == "-":
        coding_parts = get_exons_from_spectre(spectre_results, spectre_id)
        for i in range(0, len(coding_parts)):
            if i==0 and i == len(coding_parts) - 1:
                stop = coding_parts[i]['start']
                starts_list = str(coding_parts[i]['start'])
                start = coding_parts[i]['stop']
                ends_list = str(coding_parts[i]['stop'])
            elif i==0:
                stop = coding_parts[i]['start']
                starts_list = str(coding_parts[i]['start'])
                ends_list = str(coding_parts[i]['stop'])
            elif i == len(coding_parts) - 1:
                start = coding_parts[i]['stop']
                starts_list = starts_list + "_" + str(coding_parts[i]['start'])
                ends_list = ends_list + "_" + str(coding_parts[i]['stop'])
            else:
                starts_list = starts_list+"_"+str(coding_parts[i]['start'])
                ends_list = ends_list+"_"+str(coding_parts[i]['stop'])

    return start, stop, starts_list, ends_list, coding_parts

## Parse SPECtre coding parts of exons, forward strand
def get_exons_from_spectre(spectre_results, spectre_id):

    #Init
    coding_parts = defaultdict(lambda: defaultdict())

    parts = spectre_results[spectre_id]['coordinates_CDS'].split(',')
    exon_number = 0
    for i in range(0, len(parts)):
        start_stop = map(int,parts[i].split('-'))
        if exon_number == 0: #First exon (seen from 5')
            coding_parts[exon_number]['start'] = start_stop[0] + 1 #0based to 1based coordinates
            coding_parts[exon_number]['stop'] = start_stop[1]
            exon_number+=1
        else:
            #Check if spectre exon is right behind the previous one
            #In that case, concat exons as one exon => change stop of previous exon for stop of this exon
            if start_stop[0] == coding_parts[exon_number-1]['stop'] or start_stop[0] == coding_parts[exon_number-1]['stop'] + 1:
                coding_parts[exon_number-1]['stop'] = start_stop[1]
            #Check if spectre exon (mostly stop codon) is completely in the previous exon
            elif coding_parts[exon_number-1]['start'] <= start_stop[0] + 1 <= coding_parts[exon_number-1]['stop'] and \
                    coding_parts[exon_number - 1]['start'] <= start_stop[1] <= coding_parts[exon_number - 1]['stop']:
                continue
            #Else: just add another exon
            else:
                coding_parts[exon_number]['start'] = start_stop[0] + 1 #0based to 1based coordinates
                coding_parts[exon_number]['stop'] = start_stop[1]
                exon_number += 1

    #extra sorting step as SPECtre does not always put the exons in the right order in its output
    #Sort keys based on start of the coding parts
    sorted_keys = sorted(coding_parts.keys(), key=lambda x: coding_parts[x]['start'])
    #Place coding parts in a sorted nested defaultdict
    sorted_coding_parts = defaultdict(lambda: defaultdict())
    for i in range(0, len(coding_parts)):
        sorted_coding_parts[i] = coding_parts[sorted_keys[i]]

    return sorted_coding_parts

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

## Load SPECtre results file
def load_spectre_results(spectre_tmp, merged_results_file_with_header):

    #Init
    results = defaultdict(lambda: defaultdict())

    merged_results_file_without_header = spectre_tmp + "/merged_results_without_header.txt"
    #Remove header lines with SPECtre run info
    os.system("grep -v '^#' "+merged_results_file_with_header+" > "+merged_results_file_without_header)

    with open(merged_results_file_without_header, 'r') as FR:
        #Read line with column headers
        first_line = FR.readline()
        first_line.rstrip("\n")
        col_names = first_line.split("\t")
        #Read the other lines
        for line in FR:
            line.rstrip("\n")
            line_data = line.split("\t")
            for i in range(1, len(col_names)):
                results[line_data[0]][col_names[i]] = line_data[i]

    return results

## Run SPECtre
def run_spectre(installation_dir, bam, spectre_tmp, isoform_file, gtf, offsets, chrs, cores, threads_per_chromosome, starttime):

    if not os.path.isfile(bam+".bai"):
        #Index bam file
        command_index = "samtools index "+bam
        print "\tIndex untreated BAM file"
        print "\t"+command_index
        sys.stdout.flush()
        os.system(command_index)

    #Make results dir
    results_dir = spectre_tmp+"/spectre_results"
    if not os.path.isdir(results_dir):
        os.system("mkdir "+results_dir)

    #Calculate how many chrs can run in parallel
    nr_of_parallel_chrs = int(math.floor(cores/threads_per_chromosome))

    #Start multiprocessing over chromosomes
    print "\tRun SPECtre on "+str(cores)+" cores. "+str(nr_of_parallel_chrs)+" chrs parallel with "+str(threads_per_chromosome)+\
        " threads per chromosome"
    sys.stdout.flush()

    pool = Pool(processes=nr_of_parallel_chrs)
    [pool.apply_async(run_spectre_chr, args=(chrm, installation_dir, bam, spectre_tmp, isoform_file, gtf, offsets, threads_per_chromosome, starttime)) for chrm in chrs.keys()]
    pool.close()
    pool.join()

    #Merge chromosomal results into one big file. SPECtre score is in the 15th column if no FLOSS, full or ORF score is calculated
    merged_results_file = spectre_tmp + "/merged_results.txt"
    if not os.path.isfile(merged_results_file):
        command_merge = "python "+installation_dir+"/scripts/merge_split_spectre_files.py --dir "+results_dir+"/ --spec " \
                        "21 --log "+spectre_tmp+"/spectre_merge.log --full 25 --floss 29 --orfscore 30 > "+merged_results_file
        print "\t[%s]: Merge chromosomal result files" % (convert_time(time.time() - starttime))
        print "\t"+command_merge
        sys.stdout.flush()
        os.system(command_merge)
        print "\t[%s]: COMPLETE: Merge" % (convert_time(time.time() - starttime))
        sys.stdout.flush()
    else:
        print "\tMerged results file already exists"
        sys.stdout.flush()

    return merged_results_file

## Run SPECtre (chromosomal)
def run_spectre_chr(chrm, installation_dir, bam, spectre_tmp, isoform_file, gtf, offsets, threads_per_chromosome, starttime):
    try:
        if not os.path.isfile(spectre_tmp+"/spectre_results/spectre_results_"+chrm+".txt"):
            command_spectre = "python " + installation_dir + "/SPECtre.py --input " + bam + " --output " + spectre_tmp + \
                              "/spectre_results/spectre_results_" + chrm + ".txt --fpkm " + isoform_file + " --gtf " + gtf + " --log " + \
                              spectre_tmp + "/spectre_" + chrm + ".log --target " + chrm + " --offsets " + offsets + " --nt " + \
                              str(threads_per_chromosome) + " --floss --orfscore --full"
            print "\t\t[%s]: Run SPECtre for chromosome %s" % ((convert_time(time.time() - starttime)), chrm)
            print "\t\t"+command_spectre
            sys.stdout.flush()
            os.system(command_spectre)
            print "\t\t\t ** [%s]: SPECtre for chromosome %s finished **" % ((convert_time(time.time() - starttime)), chrm)
            sys.stdout.flush()
        else:
            print "\t\tSPECtre results file for chromosome "+chrm+" already exists"
            sys.stdout.flush()
    except:
        traceback.print_exc()

    return

## Download and install SPECtre
def download_spectre(spectre_tmp, link):

    installation_dir = spectre_tmp + "/spectre-master"
    if not os.path.isfile(spectre_tmp+"/spectre-master/SPECtre.py"):
        #Download from github
        os.system("wget -q --no-check-certificate "+link)
        #unzip
        os.system("unzip -q master.zip")
        os.system("rm -rf master.zip")
        #Move to tmp folder
        os.system("mv spectre-master "+spectre_tmp)
    else:
        print "SPECtre installation already present"

    return installation_dir

## Generate isoform file from translated transcripts table
def generate_iso_file(db, spectre_tmp, tr_calling):

    #Init
    transcript_table=""

    #Delete existing file
    iso_file = spectre_tmp+"/isoforms.fpkm_tracking"
    if os.path.isfile(iso_file):
        os.system("rm -rf "+iso_file)

    #Translated transcript table nomenclature
    if tr_calling=="rule_based":
        transcript_table="tr_translation"
    elif tr_calling=="ribozinb":
        transcript_table="tr_translation_ribozinb"
    else:
        print "Transcript calling is neither 'rule_based' nor 'ribozinb', but '"+tr_calling+"'"
        sys.exit()

    #Get translated transcript info
    query = "SELECT stable_id, gene_stable_id, chr, seq_region_start, seq_region_end, FPKM FROM "+transcript_table+";"
    tr_transcripts = fetchDict(db, query)

    #Open isoform file
    with open(iso_file, 'w') as FW:
        #Write header
        header = "tracking_id\tclass_code\tnearest_ref_id\tgene_id\tgene_short_name\ttss_id\tlocus\tFPKM\n"
        FW.write(header)
        #Write transcript lines
        for transcript_id in tr_transcripts:
            locus = tr_transcripts[transcript_id]['chr']+":"+str(int(round(tr_transcripts[transcript_id]['seq_region_start'], 0)))+"-"+str(int(round(tr_transcripts[transcript_id]['seq_region_end'],0)))
            line = transcript_id+"\t-\t-\t"+tr_transcripts[transcript_id]['gene_stable_id']+"\t-\t-\t"+locus+"\t"+str(tr_transcripts[transcript_id]['FPKM'])+"\n"
            FW.write(line)

    return iso_file

## Fetch hash from sqlite query ##
# The first column should contain the ID
def fetchDict(dbpath, query):
    #Init
    output = defaultdict(lambda: defaultdict())

    #Make connection to DB
    try:
        con = sqlite.connect(dbpath)
        con.text_factory = str
    except:
        print "Could not connect to "+dbpath
        sys.exit()

    with con:
        cur = con.cursor()
        if cur.execute(query):
            #Fetch all data from query
            query_output = cur.fetchall()
            colnames = list(map(lambda x: x[0], cur.description))
        #Restructure data into a nested dictionary
        for i in range(0, len(query_output)):
            id = query_output[i][0]
            for j in range(1,len(colnames)):
                col = colnames[j]
                value = query_output[i][j]
                output[id][col]=value
    return output

def get_offsets_as_list(db):

    offset_str=""

    #Make result db connection
    try:
        con = sqlite.connect(db)
    except:
        print "ERROR: could not connect to "+db
        sys.exit()

    with con:
        cur = con.cursor()

        query = "SELECT * FROM p_offsets_untreated ORDER BY length;"
        if cur.execute(query):
            output = cur.fetchall()
            for i in range(0, len(output)):
                length = output[i][0]
                offset = output[i][1]
                if offset_str == "":
                    offset_str = str(length)+":"+str(offset)
                else:
                    offset_str = offset_str+","+str(length)+":"+str(offset)
        else:
            print "ERROR: could not fetch the untreated offsets from "+db
            print query
            sys.exit()

    return offset_str

## Get number of cores from arguments table
def get_cores(db):

    #Make result db connection
    try:
        con = sqlite.connect(db)
    except:
        print "ERROR: could not connect to "+db
        sys.exit()

    with con:
        cur = con.cursor()

        query = "SELECT value FROM arguments WHERE variable='nr_of_cores';"
        if cur.execute(query):
            cores = int(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the number of cores argument from "+db
            print query
            sys.exit()

    return cores

## Get untreated BAM file from arguments table
def get_bam_from_arguments(db):

    #Make result db connection
    try:
        con = sqlite.connect(db)
    except:
        print "ERROR: could not connect to "+db
        sys.exit()

    with con:
        cur = con.cursor()

        query = "SELECT value FROM arguments WHERE variable='out_bam_untr';"
        if cur.execute(query):
            bam = str(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the output untreated BAM argument from "+db
            print query
            sys.exit()

    return bam

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

        query = "SELECT value FROM arguments WHERE variable='tr_calling';"
        if cur.execute(query):
            tr_calling_method = str(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the tr_calling argument in "+db
            print query
            sys.exit()

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

    return tr_calling_method, igenomes_root, species, version, ens_db

def convert_time(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)

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
    pool = Pool(processes=cores)
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
                if re.search("^>", line):
                    continue
                else:
                    line.rstrip(chars='\n')
                    FW.write(line)

    return

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