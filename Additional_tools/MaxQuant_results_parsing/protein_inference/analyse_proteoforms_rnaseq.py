import traceback
import sys
import os
import time
import argparse
import re
from collections import defaultdict
import sqlite3 as sqlite
#from dict_functions import dict_funcs
import pandas as pd
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import ConnectionPatch, ArrowStyle

__author__ = 'Steven Verbruggen'
#Execute 'python analyse_proteoforms.py -h' for more information

def main():

    starttime = time.time()

    print
    print "#######################"
    print "# Analyze proteoforms #"
    print "#######################"
    print
    print "This program is part of the PROTEOFORMER pipeline"
    print

    #Parse arguments from command line
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS,
                                     description="This tool is part of the PROTEOFORMER pipeline. It analyses "
                                                 "the different proteoforms found after MS analysis of a combined database "
                                                 "of PROTEOFORMER and UniProt. It searches for the classification of all "
                                                 "proteoforms found by PROTEOFORMER, not yet in UniProt and with confiramtion"
                                                 " of MS with MaxQuant.",
                                     add_help=False)

    man_args = parser.add_argument_group("Mandatory parameters")
    man_args.add_argument("--max_protein_db", "-m", action="store", required=True, nargs="?", metavar="DB",
                          type=str, help="Max protein DB (mandatory)")
    man_args.add_argument("--ens_db", "-e", action="store", required=True, nargs="?", metavar="DB",
                          type=str, help="Ensembl database (mandatory)")
    man_args.add_argument("--mapping_file", "-M", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="Uniprot mapping file. You can get that from the Uniprot website, downloads, "
                                         "ID mapping (mandatory)")
    man_args.add_argument("--fasta_file", "-f", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="Combined fasta file of proteoformer and uniprot (mandatory)")
    man_args.add_argument("--protein_groups", "-g", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="MaxQuant output file of protein groups (mandatory)")
    man_args.add_argument("--peptides", "-p", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="MaxQuant output file of peptides(mandatory)")

    opt_args = parser.add_argument_group("Optional parameters")
    opt_args.add_argument("--help", "-h", action="help", help="Show help message and exit")
    opt_args.add_argument("--workdir", "-w", action="store", required=False, nargs="?", metavar="FOLDER", default="",
                        type=str, help="Working directory (default: CWD)")
    opt_args.add_argument("--csv_output", "-x", action="store", required=False, nargs="?", metavar="PATH",
                          type=str, default="proteoform_classifications.csv", help="CSV output file of classifications (default: proteoform_classifications.csv)")

    args = parser.parse_args()

    #default workdir is CWD
    if args.workdir == "":
        args.workdir = os.getcwd()

    #List parameters
    print "Parameters:"
    for arg in vars(args):
        print '    %-15s\t%s' % (arg, getattr(args, arg))
    print
    sys.stdout.flush()

    ########
    # MAIN #
    ########

    #Load proteoformer-only identifications from max protein MS results. Herein, possible new proteoforms could be included.
    identifications = load_proteoformer_identifications(args.max_protein_db)

    #Make a Ensembl transcript id - Uniprot id mapping dict
    mapping_dict = construct_mapping_dict(args.mapping_file)

    #Get all sequences from the fasta file
    fasta_sequences = load_fasta(args.fasta_file)

    #Read in protein groups
    protein_groups = read_proteingroups(args.protein_groups)

    #Read in peptides
    peptides = read_peptides(args.peptides)

    ##############
    # Test cases #
    ##############
    #N-terminal extension
    #test_id = "ENST00000373382_10_69124218_5UTR_110db1;ENST00000395098_10_69124218_5UTR_100db1;ENST00000263559_10_69124230_5UTR_100db1;O75436;S4R3Q6;ENST00000489794_10_69132940_5UTR_100db1;ENST00000395098_10_69124278_aTIS_101db1;ENST00000395098_10_69124230_5UTR_100db1;ENST00000490696_10_69124218_ntr_100db1;ENST00000263559_10_69157162_CDS_100db1"
    #test_id = "ENST00000409021_4_6112850_aTIS_100db1;A0A1B0GUE0;F2Z2K5;Q96N16"
    #test_id="ENST00000446526_17_76353703_aTIS_101db1;ENST00000446526_17_76353814_5UTR_100db1;ENST00000446526_17_76353856_5UTR_100db1;Q14558;ENST00000446526_17_76332389_CDS_100db1;B4DP31;ENST00000324684_17_76344678_5UTR_110db1;C9J168;ENST00000446526_17_76328792_CDS_100db1;ENST00000472686_17_76328792_ntr_100db1;C9JNQ3;C9JUN4;C9JKT9"
    #N-terminal truncation
    #test_id = "ENST00000338193_12_56751067_CDS_100db1;ENST00000338193_12_56743008_CDS_100db1"
    #test_id = "ENST00000281456_4_185143475_CDS_010db2"
    #N-terminal truncation without confirming peptide
    #test_id = "ENST00000593845_19_15397845_ntr_100db1;ENST00000397410_19_15397854_CDS_100db1;ENST00000397410_19_15401290_CDS_010db2;ENST00000595465_19_15418923_aTIS_101db1;ENST00000397410_19_15418923_aTIS_101db1;Q9ULX6"
    #NTR
    #test_id = "ENST00000506634_4_170605301_ntr_100db1"
    #test_id = "ENST00000428881_1_171683767_ntr_100db1"
    #3UTR
    #test_id = "ENST00000425660_7_5528626_3UTR_100db1"
    #Exon exclusion
    #test_id = "ENST00000233596_19_1491270_aTIS_111db1;ENST00000233596_19_1491222_5UTR_100db1;ENST00000233596_19_1491213_5UTR_100db1;ENST00000233596_19_1496294_CDS_100db1;ENST00000233596_19_1495575_CDS_100db1;Q96HR9;ENST00000395479_19_1496294_CDS_100db1;ENST00000395479_19_1495575_CDS_100db1;A8MXN1;A8MWX0;ENST00000591735_19_1491222_ntr_100db1;ENST00000591735_19_1491213_ntr_100db1"
    #uORF
    #test_id = "ENST00000419065_6_3258968_5UTR_100db1;ENST00000438998_6_3258968_5UTR_100db1;ENST00000473000_6_3258968_5UTR_100db1"
    #out of frame ORF
    #test_id = "ENST00000244020_20_43460226_CDS_100db1"
    #test_id = "ENST00000412249_3_51674477_CDS_100db1;ENST00000614067_3_51674477_CDS_100db1"
    #C terminal truncation
    #test_id="ENST00000398970_17_1456117_aTIS_101db1"
    #C terminal extension
    #test_id="ENST00000338435_2_190881085_aTIS_101db1;O94925;H7BZD1;B8ZZC5;ENST00000469774_2_190900619_ntr_100db1;B8ZZA8;ENST00000409626_2_190927321_5UTR_100db1"
    #test_id="ENST00000372798_1_40059347_aTIS_101db1;ENST00000340450_1_40040797_5UTR_110db1;ENST00000372798_1_40040797_5UTR_100db1;ENST00000340450_1_40060135_CDS_100db1;ENST00000479759_1_40067567_ntr_100db1;ENST00000479759_1_40069737_ntr_100db1;Q5T0R9;ENST00000479759_1_40069851_ntr_100db1;ENST00000476047_3_76434856_ntr_100db1;ENST00000417455_10_43605301_ntr_100db1;ENST00000476047_3_76434832_ntr_100db1;ENST00000479759_1_40070450_ntr_100db1;Q5T0S3;Q5T0R8;ENST00000493172_6_17421556_aTIS_101db1;B7Z385;A0A087X0J3;A0A087WZ15;ENST00000465994_6_17421556_aTIS_101db1;E9PDI2;P40123"
    #N terminal splice variant
    #test_id="ENST00000319410_16_75644428_aTIS_101db1;ENST00000302445_16_75641719_CDS_100db1;ENST00000302445_16_75647609_CDS_010db2;Q15046;ENST00000566560_16_75640313_ntr_100db1;H3BVA8;ENST00000566560_16_75641719_ntr_100db1;ENST00000569298_16_75628619_ntr_100db1;J3KRL2;H3BPV7;ENST00000566772_16_75636073_aTIS_101db1;H3BMR9;ENST00000566772_16_75635974_CDS_100db1;ENST00000566772_16_75635989_CDS_100db1;ENST00000566560_16_75635974_ntr_100db1;H3BQK5;ENST00000566560_16_75635989_ntr_100db1"
    #N-terminal truncation with SAV -> multiple variations
    #test_id = "ENST00000593845_19_15397845_ntr_100db1;ENST00000397410_19_15397854_CDS_100db1;ENST00000397410_19_15401290_CDS_010db2;ENST00000595465_19_15418923_aTIS_101db1;ENST00000397410_19_15418923_aTIS_101db1;Q9ULX6-2;Q9ULX6"

    #More complex variation
    #test_id = "ENST00000345136_8_143939461_aTIS_101db1;ENST00000398774_8_143944663_aTIS_101db1;ENST00000354958_8_143953771_aTIS_101db1;ENST00000356346_8_143973472_aTIS_101db1;ENST00000354589_8_143943890_aTIS_101db1;ENST00000357649_8_143942515_aTIS_101db1;ENST00000436759_8_143975369_aTIS_101db1;ENST00000527096_8_143975369_aTIS_101db1;E9PMV1;H0YDN1;E9PKG0;E9PIA2;E9PQ28;REV__A0A0U1RR03"
    #test_id="ENST00000301329_17_782255_aTIS_101db1;F6TLX2;Q9HC38;I3L3Q4;ENST00000571073_17_776957_ntr_110db1;ENST00000571073_17_771387_ntr_100db1;ENST00000571073_17_771393_ntr_100db1;ENST00000301329_17_771435_CDS_100db1;ENST00000576239_17_776957_ntr_100db1;I3L1F4;ENST00000573137_17_771393_CDS_100db1;I3L1I0;I3NI27;I3NI24;I3L2C2;I3L277"


    #identifications_copy = defaultdict(lambda: defaultdict())
    #identifications_copy[test_id] = identifications[test_id]
    #identifications = identifications_copy

    #dict_funcs.print_dict(identifications[test_id])

    #Analyze identifications and parse them in the right categories
    identifications = analyze_identifications(args.workdir, identifications, mapping_dict, fasta_sequences,protein_groups, peptides, args.ens_db)

    #Output identifications and counts as csv
    output_csv(identifications, args.csv_output)


    # End of program message
    print
    print "-----------------------"
    print "[%s]: PROGRAM COMPLETE" % (convert_time(time.time() - starttime))
    print "-----------------------"
    print
    sys.stdout.flush()

    return

##########
#  SUBS  #
##########

#Output identifications
def output_csv(identifications, csv_file):

    id=0
    with open(csv_file, 'w') as FW:
        keys = ["base_proteoform","proving_peptides","peptide_PEP","max_proteins","protein_group"]
        header_string= "id"+","+(','.join(keys))+"\n"
        FW.write(header_string)
        for protein_group in identifications:
            value_string = str(identifications[protein_group][keys[0]])
            for key in keys[1:]:
                if type(identifications[protein_group][key]) is list:
                    value_string = value_string+','+str("|".join(identifications[protein_group][key]))
                else:
                    if re.search(', ', str(identifications[protein_group][key])):
                        value = re.sub(', ', '|', str(identifications[protein_group][key]))
                        value_string = value_string+','+value
                    else:
                        value_string = value_string+','+str(identifications[protein_group][key])
            FW.write(str(id)+','+value_string+"\n")
            id+=1

    return

#Analyze identifications
def analyze_identifications(workdir, identifications, mapping_dict, fasta_sequences, protein_groups, peptides, ens_db):

    #For all protein_groups
    for protein_group in identifications:
        #Get the most important proteoform
        base_proteoform = find_base_proteoform(identifications, protein_group)
        identifications[protein_group]['base_proteoform'] = base_proteoform
        #Get the most important canonical protein
        base_canonical, canonical_found_in_ensembl = find_base_canonical(protein_group, base_proteoform, mapping_dict)
        #Get sequences
        base_proteoform_seq = fasta_sequences[base_proteoform]
        #Get peptide info
        identifications[protein_group]['proving_peptides'], identifications[protein_group]['peptide_PEP'] = \
                peptide_info(identifications, protein_group, protein_groups, peptides, base_proteoform_seq, "")


    return identifications


#Get peptide info for in csv of proving peptides
def peptide_info(identifications, protein_group, protein_groups, peptides, proteoform_seq, base_canonical):

    proving_peptides = []
    peptide_PEP = []
    non_canonical_peptides = []
    non_canonical_PEP = []

    print("Base canonical: "+base_canonical)

    peptide_ids = []
    peptide_ids.extend(protein_groups[protein_group]['Peptide IDs'].split(';'))
    for peptide_id in peptide_ids:
        print("Peptide ID: "+peptide_id)
        if peptide_id in peptides: #reverses and contaminants were not loaded from the files
            peptide_seq = peptides[peptide_id]['Sequence']
            #Check if peptide surely in proteoform sequence
            if peptide_seq in proteoform_seq:
                #Check if peptide also in other canonical sequences
                peptide_proteins = []
                peptide_proteins.extend(peptides[peptide_id]['Proteins'].split(';'))
                print(peptide_proteins)
                uniprot_proteins_present=False
                for peptide_protein in peptide_proteins:
                    m = re.search('^ENST', peptide_protein)
                    if not m:
                        uniprot_proteins_present = True
                if uniprot_proteins_present==False:
                    proving_peptides.append(peptide_seq)
                    peptide_PEP.append(peptides[peptide_id]['PEP'])
                if base_canonical not in peptide_proteins:
                    non_canonical_peptides.append(peptide_seq)
                    non_canonical_PEP.append(peptides[peptide_id]['PEP'])
    if len(proving_peptides)==0: #Sometimes a new proteoform overlaps partly with one canonical and another, then search for the base canonical and compare with this one
        proving_peptides = non_canonical_peptides
        peptide_PEP = non_canonical_PEP

    return proving_peptides, peptide_PEP

#Search longest identical stretch
def search_longest_stretch(proteoform_seq, canonical_seq):

    #Init
    longest_stretch = ""
    cur_stretch = ""


    for pos in range(0, len(proteoform_seq)):
        #Check if identical so that stretch should be expanded
        if proteoform_seq[pos]==canonical_seq[pos]:
            cur_stretch = cur_stretch + proteoform_seq[pos]
        else:
            cur_stretch = ""
        #Check if it is the longest stretch so far
        if len(cur_stretch)>len(longest_stretch):
            longest_stretch = cur_stretch

    #Calculate percentage covered by longest stretch
    percentage = float(len(longest_stretch))/float(len(canonical_seq))

    return longest_stretch, percentage

#Check only SAV
def check_only_sav(aligned_proteoform_seq, aligned_canonical_seq):

    #Init
    only_sav=True
    prev_not_matched = False

    for pos in range(0, len(aligned_proteoform_seq)):
        if prev_not_matched==False:
            if aligned_proteoform_seq[pos]!=aligned_canonical_seq[pos]:
                prev_not_matched=True
        else:
            if aligned_proteoform_seq[pos]==aligned_canonical_seq[pos]:
                prev_not_matched=False
            else:
                only_sav=False
                break

    return only_sav

#Check if proteoform start is in phase
def check_if_in_phase(base_proteoform_start, exon_info, strand):

    #Init
    in_phase=""

    if strand==1:
        for rank in sorted(exon_info.keys()):
            #Search the exon in which the start is in
            if base_proteoform_start>=exon_info[rank]['start'] and base_proteoform_start<=exon_info[rank]['stop']:
                #Check if in same phase
                phase_diff = (base_proteoform_start - exon_info[rank]['start']) % 3
                if phase_diff==0:
                    in_phase='Y'
                elif phase_diff==1 or phase_diff==2:
                    in_phase='N'
    elif strand==-1:
        for rank in sorted(exon_info.keys()):
            if base_proteoform_start<=exon_info[rank]['start'] and base_proteoform_start>=exon_info[rank]['stop']:
                phase_diff = (exon_info[rank]['start'] - base_proteoform_start) % 3
                if phase_diff==0:
                    in_phase='Y'
                elif phase_diff==1 or phase_diff==2:
                    in_phase='N'

    return in_phase

#Get exon info of canonical ORF
def get_canonical_exon_info(transcript_id, canonical_start, canonical_stop, strand, ens_db):

    #Init
    exon_info = defaultdict(lambda: defaultdict())
    transcript_exons = defaultdict(lambda: defaultdict())

    #Connect to ensembl
    try:
        conn = sqlite.connect(ens_db)
    except:
        print "Could not connect to "+ens_db
        sys.exit()

    with conn:
        cur = conn.cursor()

        #Get all transcript exons
        query = "SELECT et.rank, e.seq_region_start, e.seq_region_end, e.phase, e.end_phase FROM exon AS e " \
                "JOIN exon_transcript AS et ON e.exon_id = et.exon_id JOIN transcript AS tr ON "\
                "et.transcript_id = tr.transcript_id WHERE tr.stable_id = '"+transcript_id+"' ORDER BY et.rank;"
        cur.execute(query)
        output = cur.fetchall()

        #Parse into defaultdict
        for i in range(0,len(output)):
            transcript_exons[output[i][0]]['seq_region_start'] = output[i][1]
            transcript_exons[output[i][0]]['seq_region_end'] = output[i][2]
            transcript_exons[output[i][0]]['phase'] = output[i][3]
            transcript_exons[output[i][0]]['end_phase'] = output[i][4]

        #Only keep the translated exon parts. Go from absolute ensembl exon stradness to relative
        if strand==1:
            for rank in sorted(transcript_exons.keys()):
                #Both in this exon: trim both sides
                if canonical_start>=transcript_exons[rank]['seq_region_start'] and canonical_stop<=transcript_exons[rank]['seq_region_end']:
                    exon_info[rank]['start'] = canonical_start
                    exon_info[rank]['phase'] = 0
                    exon_info[rank]['stop'] = canonical_stop
                    exon_info[rank]['stop_phase'] = 2
                #Canonical start in exon: trim
                elif canonical_start>=transcript_exons[rank]['seq_region_start'] and canonical_start<=transcript_exons[rank]['seq_region_end']:
                    exon_info[rank]['start'] = canonical_start
                    exon_info[rank]['phase'] = 0
                    exon_info[rank]['stop'] = transcript_exons[rank]['seq_region_end']
                    exon_info[rank]['stop_phase'] = transcript_exons[rank]['end_phase']
                #Canonical stop in exon: trim
                elif canonical_stop>=transcript_exons[rank]['seq_region_start'] and canonical_stop<=transcript_exons[rank]['seq_region_end']:
                    exon_info[rank]['start'] = transcript_exons[rank]['seq_region_start']
                    exon_info[rank]['phase'] = transcript_exons[rank]['phase']
                    exon_info[rank]['stop'] = canonical_stop
                    exon_info[rank]['stop_phase'] = 2
                #No canonical start or stop in exon, but inside canonical translation zone: just copy
                elif canonical_start<transcript_exons[rank]['seq_region_start'] and canonical_stop>transcript_exons[rank]['seq_region_end']:
                    exon_info[rank]['start'] = transcript_exons[rank]['seq_region_start']
                    exon_info[rank]['phase'] = transcript_exons[rank]['phase']
                    exon_info[rank]['stop'] = transcript_exons[rank]['seq_region_end']
                    exon_info[rank]['stop_phase'] = transcript_exons[rank]['end_phase']
        elif strand==-1:
            for rank in sorted(transcript_exons.keys()):
                #Both in this exon: trim both sides
                if canonical_start<=transcript_exons[rank]['seq_region_end'] and canonical_stop>=transcript_exons[rank]['seq_region_start']:
                    exon_info[rank]['start'] = canonical_start
                    exon_info[rank]['phase'] = 0
                    exon_info[rank]['stop'] = canonical_stop
                    exon_info[rank]['stop_phase'] = 2
                #Canonical start in exon: trim
                elif canonical_start<=transcript_exons[rank]['seq_region_end'] and canonical_start>=transcript_exons[rank]['seq_region_start']:
                    exon_info[rank]['start'] = canonical_start
                    exon_info[rank]['phase'] = 0
                    exon_info[rank]['stop'] = transcript_exons[rank]['seq_region_start']
                    exon_info[rank]['stop_phase'] = transcript_exons[rank]['phase']
                #Canonical stop in exon: trim
                elif canonical_stop<=transcript_exons[rank]['seq_region_end'] and canonical_stop>=transcript_exons[rank]['seq_region_start']:
                    exon_info[rank]['start'] = transcript_exons[rank]['seq_region_end']
                    exon_info[rank]['phase'] = transcript_exons[rank]['end_phase']
                    exon_info[rank]['stop'] = canonical_stop
                    exon_info[rank]['stop_phase'] = 2
                #No canonical start or stop in exon but in canonical translation zone: just copy
                elif canonical_start>transcript_exons[rank]['seq_region_end'] and canonical_stop<transcript_exons[rank]['seq_region_start']:
                    exon_info[rank]['start'] = transcript_exons[rank]['seq_region_end']
                    exon_info[rank]['phase'] = transcript_exons[rank]['end_phase']
                    exon_info[rank]['stop'] = transcript_exons[rank]['seq_region_start']
                    exon_info[rank]['stop_phase'] = transcript_exons[rank]['phase']

    return exon_info

#Get start of proteoform
def get_start(proteoform):

    #Init
    start = ""

    m = re.search('^ENST\d+?\_.+?\_(\d+?)\_', proteoform)
    if m:
        start = int(m.group(1))

    return start

#Get canonical start and stop of a transcript
def get_canonical_start_stop(transcript_id, ens_db):

    #Init
    can_start = 0
    can_stop  = 0

    #Connect to ensembl
    conn = sqlite.connect(ens_db)

    with conn:
        cur = conn.cursor()

        #Get canonical translation id
        query_transcript = "SELECT canonical_translation_id, seq_region_strand FROM transcript WHERE stable_id='"+transcript_id+"';"
        cur.execute(query_transcript)
        (canonical_translation_id, strand) = cur.fetchone()

        #Get start and stop info of canonical translation
        query_translation = "SELECT seq_start, seq_end, start_exon_id, end_exon_id FROM translation WHERE translation_id='"+str(canonical_translation_id)+"';"
        cur.execute(query_translation)
        (seq_start, seq_end, start_exon_id, end_exon_id) = cur.fetchone()

        #Get starts and ends of start and end exons
        query_start_exon = "SELECT seq_region_start, seq_region_end FROM exon where exon_id='"+str(start_exon_id)+"';"
        cur.execute(query_start_exon)
        (start_exon_start, start_exon_end) = cur.fetchone()
        query_end_exon = "SELECT seq_region_start, seq_region_end FROM exon where exon_id='"+str(end_exon_id)+"';"
        cur.execute(query_end_exon)
        (end_exon_start, end_exon_end) = cur.fetchone()

        #Calculate canonical start and end
        if strand==1:
            can_start = start_exon_start + seq_start - 1
            can_stop  = end_exon_start + seq_end - 1
        elif strand==-1:
            can_start = start_exon_end - seq_start + 1
            can_stop = end_exon_end - seq_end + 1

    return (can_start, can_stop, strand)

#Find annotation of accession
def get_annotation(proteoform):

    #Init
    annotation = ""

    m = re.search('^ENST\d+?\_.+?\_\d+?\_(.+?)\_.+', proteoform)
    if m:
        annotation = m.group(1)

    return annotation

#Find biotype of proteoform
def find_biotype(proteoform, ens_db):

    #Get transcript_id
    transcript_id = get_transcript_id(proteoform)

    #Connect to ensembl
    conn = sqlite.connect(ens_db)
    with conn:
        cur = conn.cursor()
        query = "select biotype from transcript where stable_id='"+transcript_id+"';"
        cur.execute(query)
        biotype = cur.fetchone()[0]

    return biotype

#Get transcript id out of proteoform accession
def get_transcript_id(proteoform):

    #Init
    transcript_id=""

    m = re.search('^(ENST\d+?)\_', proteoform)
    if m:
        transcript_id = m.group(1)

    return transcript_id

#Correct position based on alignment
def correct_pos(aligned_seq, org_pos):

    #subtract each time the amount of position which are not in an insert until i=0
    i = int(org_pos)
    rest_seq = aligned_seq
    tmp_seq = ""
    corrected = 0
    while i!=0:
        tmp_seq = rest_seq[:i]
        rest_seq = rest_seq[i:]
        i = tmp_seq.count('-')
        corrected = corrected + len(tmp_seq)

    return corrected

#Get peptides and their positions on the fulle protein sequence
def get_peptides(protein_groups, peptides, protein_group):

    #Init
    proteoform_peptides = defaultdict(lambda: defaultdict())

    #Get all peptides of protein group
    peptide_list = protein_groups[protein_group]['Peptide IDs'].split(';')

    #For all peptides
    for peptide in peptide_list:
        #Get start, stop position and sequence
        if peptide in peptides.keys():#Otherwise it is a contaminant or reverse peptide
            proteoform_peptides[peptide]['sequence'] = peptides[peptide]['Sequence']
            proteoform_peptides[peptide]['start'] = peptides[peptide]['Start position']
            proteoform_peptides[peptide]['end'] = peptides[peptide]['End position']
            proteoform_peptides[peptide]['proteins'] = peptides[peptide]['Proteins']

    return proteoform_peptides

#Calculate mapped percentage
def calc_perc_mapped(aligned_base_proteoform, aligned_base_canonical):

    #Init
    perc_mapped = 0
    identical_positions = 0
    length_without_indels = 0

    #Go over alignment per base
    for base_pos in range(0, len(aligned_base_proteoform)):
        #Check for indels
        if aligned_base_proteoform[base_pos]!='-' and aligned_base_canonical[base_pos]!='-':
            length_without_indels+=1
            #Check if bases are identical
            if aligned_base_proteoform[base_pos]==aligned_base_canonical[base_pos]:
                #If identical, +1 identical positions
                identical_positions+=1

    #Calculate percentage identical positions
    perc_mapped = float(identical_positions)/float(length_without_indels)

    return perc_mapped

#Read clustal output
def read_clustal_output(fasta, base_proteoform, base_canonical):

    #Init
    clustal_seqs = defaultdict()

    #Read output
    with open(fasta, 'r') as FR:
        lines = FR.readlines()
        lines = map(lambda x: x.rstrip("\n"), lines)
        lines = map(lambda x: x.rstrip("\r"), lines)
        # Init
        accession = ""
        sequence = ""
        for i in range(0, len(lines)):
            # Check if accession -> new entry
            if ((re.search('^>', lines[i]))):
                if (i != 0):
                    # Save previous entry
                    clustal_seqs[accession[1:]] = sequence
                # Get new entry accession and empty sequence
                accession = lines[i]
                sequence = ""
            else:
                # Concat next line of sequence
                sequence += lines[i]
        # Save the last entry
        if (sequence != ""):
            clustal_seqs[accession[1:]] = sequence
            accession = ""
            sequence = ""

    aligned_base_proteoform = clustal_seqs[base_proteoform]
    aligned_base_canonical = clustal_seqs[base_canonical]

    return (aligned_base_proteoform, aligned_base_canonical)


#Align with ClustalO
def clustal_align(workdir, base_proteoform, input_fa, clustalo_path):

    #Init output
    output_fa = workdir+"/"+base_proteoform+"aligned_output.fa"

    #Execute
    command = clustalo_path+" -i "+input_fa+" -o "+output_fa+" --outfmt fa --force --wrap=1000000000" #Force to overwrite, output in fasta format
    os.system(command)

    return output_fa

#Prepare tmp fasta
def prepare_fasta(workdir, base_proteoform, base_proteoform_seq, base_canonical, base_canonical_seq):

    #Define path of tmp fasta file
    tmp_fasta = workdir+"/tmp_input.fa"

    with open(tmp_fasta, 'w') as FW:
        FW.write(">"+base_proteoform+"\n")
        FW.write(base_proteoform_seq+"\n")
        FW.write(">"+base_canonical+"\n")
        FW.write(base_canonical_seq+"\n")

    return tmp_fasta

#Find the base canonical protein
def find_base_canonical(protein_group, base_proteoform, mapping_dict):

    #Init
    base_canonical = ""

    #try to find it in the protein group
    #Split by ;
    proteins = map(lambda x: str(x), protein_group.split(';'))
    for protein in proteins:
        if not re.search('^ENST', protein):
            return protein, False
    #Check if base canonical already found in protein group
    if base_canonical=="":
        #If not found there, search the base canonical protein in ensembl db with the transcript id of the base proteoform
        base_canonical = search_base_canonical_with_mapping_info(base_proteoform, mapping_dict)

    return base_canonical, True

#Try to find the base canonical with ensembl info of the base proteoform
def search_base_canonical_with_mapping_info(base_proteoform, mapping_dict):

    #Init
    base_canonical=""

    #Get the transcript id of the base proteoform
    transcript_id = ""
    m = re.search('^(ENST.+?)\_', base_proteoform)
    if m:
        transcript_id=m.group(1)

    if transcript_id in mapping_dict.keys():
        base_canonical = mapping_dict[transcript_id]

    return base_canonical

#Find the base proteoform
def find_base_proteoform(identifications, protein_group):

    #If multiple max proteins
    if isinstance(identifications[protein_group]['max_proteins'], list):
        base_proteoform = identifications[protein_group]['max_proteins'][0] #Take for now the first proteoform as base proteoform
    elif isinstance(identifications[protein_group]['max_proteins'], str) or isinstance(identifications[protein_group]['max_proteins'], unicode):
        base_proteoform = identifications[protein_group]['max_proteins']

    #print "Base proteoform: "+base_proteoform

    return base_proteoform

#Read peptides input info from maxquant
def read_peptides(in_file):

    #Init
    peptides = defaultdict(lambda: defaultdict())

    #Open file
    with open(in_file, 'r') as FR:

        #Read lines
        lines = FR.readlines()
        lines = map(lambda x: x.rstrip("\n"), lines)
        lines = map(lambda x: x.rstrip("\r"), lines)

        #Parse column headers
        col_headers = re.split('\t', lines[0])

        #Search index of contaminants and reverses
        idx_reverses = col_headers.index("Reverse")
        idx_cont = col_headers.index("Potential contaminant")

        #Parse data
        for line_nr in range(1,len(lines)):
            data = re.split('\t', lines[line_nr])
            #Remove contaminants and reverse identifications
            if(data[idx_reverses]!='+' and data[idx_cont]!='+'):
                for val_nr in range(0, len(data)):
                    peptides[data[48]][col_headers[val_nr]] = data[val_nr] #data[48] is peptide ID

    return peptides

#Read protein group input info from maxquant
def read_proteingroups(in_file):

    #Init
    protein_groups = defaultdict(lambda: defaultdict())

    #Open file
    with open(in_file, 'r') as FR:

        #Read lines
        lines = FR.readlines()
        lines = map(lambda x: x.rstrip("\n"), lines)
        lines = map(lambda x: x.rstrip("\r"), lines)

        #Parse column headers
        col_headers = re.split('\t', lines[0])

        #Search index of contaminants and reverses
        idx_reverses = col_headers.index("Reverse")
        idx_cont = col_headers.index("Potential contaminant")

        #Parse data
        for line_nr in range(1,len(lines)):
            data = re.split('\t', lines[line_nr])
            #Remove contaminants and reverse identifications
            if(data[idx_reverses]!='+' and data[idx_cont]!='+'):
                for val_nr in range(0, len(data)):
                    protein_groups[data[0]][col_headers[val_nr]] = data[val_nr]

    return protein_groups

#Load fasta file
def load_fasta(fasta):

    #Init
    fasta_sequences = defaultdict()

    with open(fasta, 'r') as FR:
        lines = FR.readlines()
        lines = map(lambda x: x.rstrip("\n"), lines)
        lines = map(lambda x: x.rstrip("\r"), lines)
        # Init
        accession = ""
        sequence = ""
        for i in range(0, len(lines)):
            # Check if accession -> new entry
            if ((re.search('^>', lines[i]))):
                if (i != 0):
                    # Save previous entry
                    id = get_id_from_accession(accession)
                    fasta_sequences[id] = sequence
                # Get new entry accession and empty sequence
                accession = lines[i]
                sequence = ""
            else:
                # Concat next line of sequence
                sequence += lines[i]
        # Save the last entry
        if (sequence != ""):
            id = get_id_from_accession(accession)
            fasta_sequences[id] = sequence
            accession = ""
            sequence = ""

    return fasta_sequences

#Get proteoform id out of accession
def get_id_from_accession(accession):

    #Init
    id=""

    m = re.search('^.+?\|(.+?)\|.+', accession)
    if m:
        id = m.group(1)

    return id

#Construct ID mapping dict
def construct_mapping_dict(mapping_file):

    #Init
    #mapping_dict[ensembl transcript id] = uniprot accession
    mapping_dict = defaultdict()

    with open(mapping_file, 'r') as FR:
        for line in FR:
            values = line.split('\t')
            transcript_ids = values[19].split('; ')
            for transcript_id in transcript_ids:
                mapping_dict[transcript_id] = values[0]


    return mapping_dict

#Load proteoformer-only identifications
def load_proteoformer_identifications(max_protein_db):

    #Init
    #identifications->{protein_group}->{column} = value
    identifications = defaultdict(lambda: defaultdict())

    #DB connection
    try:
        conn = sqlite.connect(max_protein_db)
    except:
        print "Could not connect to "+max_protein_db
        sys.exit()

    with conn:
        cur = conn.cursor()

        #Get all protein groups only identified with a max protein by Proteoformer
        query = "SELECT protein_group, max_proteins FROM max_proteins WHERE sources='Proteoformer';"
        cur.execute(query)

        colnames = list(map(lambda x: x[0], cur.description))

        data = cur.fetchall()
        for i in range(0, len(data)):
            id = data[i][0]
            for j in range(0, len(colnames)):
                col = colnames[j]
                value = data[i][j]
                if re.search('\|', value):
                    identifications[id][col] = map(lambda x: str(x), value.split('|'))
                else:
                    identifications[id][col] = value

    return identifications

def convert_time(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)

###### DIRECT TO MAIN ############
if __name__ == "__main__":
    try:
        main()
    except Exception, e:
        traceback.print_exc()
##################################