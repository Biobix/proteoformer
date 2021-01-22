#!/usr/bin/env python
import traceback
import sys
import os
import time
import argparse
from collections import defaultdict
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import sqlite3 as sqlite

#from dict_functions import dict_funcs


__author__ = 'Steven Verbruggen'
#Execute 'python parse_maxquant_uniprot.py -h' for more information

def main():

    starttime = time.time()

    print
    print "###################################"
    print "# Parse MaxQuant output (UniProt) #"
    print "###################################"
    print
    print "This program is part of the PROTEOFORMER pipeline"
    print

    #Parse arguments from command line
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS,
                                     description="This tool is part of the PROTEOFORMER pipeline. It parses the"
                                                 "output of a MaxQuant MS validation (uniprot).",
                                     add_help=False)

    man_args = parser.add_argument_group("Mandatory parameters")
    man_args.add_argument("--protein_groups", "-g", action="store", required=True, nargs="?", metavar="PATH",
                          default="", type=str, help="Input file with all MaxQuant identified protein groups (mandatory)")
    man_args.add_argument("--peptides", "-p", action="store", required=True, nargs="?", metavar="PATH",
                           default="", type=str, help="Input file with all MaxQuant identified peptides (mandatory)")
    man_args.add_argument("--psms", "-x", action="store", required=True, nargs="?", metavar="PATH",
                          default="", type=str, help="Input file with all MaxQuant MS/MS spectra (mandatory)")
    man_args.add_argument("--combined_fasta", "-f", action="store", required=True, nargs="?", metavar="PATH",
                          default="", type=str, help="Combined fasta file on which MaxQuant was run. Underlying "
                                                     "fasta files should be generated with verbose output (mandatory)")

    opt_args = parser.add_argument_group("Optional parameters")
    opt_args.add_argument("--help", "-h", action="help", help="Show help message and exit")
    opt_args.add_argument("--workdir", "-w", action="store", required=False, nargs="?", metavar="FOLDER", default="",
                        type=str, help="Working directory (default: CWD)")
    opt_args.add_argument("--score_file", "-s", action="store", required=False, nargs="?", metavar="PATH",
                          default="scores_uniprot.txt", type=str, help="Output tab separated score file"
                                                                              " (default: level_scores_uniprot.txt)")
    opt_args.add_argument("--venn_file", "-v", action="store", required=False, nargs="?", metavar="PATH",
                          default="venn_diagram_uniprot.png", type=str, help="Path to a png file where the Venn diagram will "
                                                                     "be stored (default: level_venn_diagram_uniprot.png)")
    opt_args.add_argument("--max_protein_db", "-d", action="store", required=False, nargs="?", metavar="PATH",
                          default="max_proteins_uniprot.db", type=str, help="SQLite db path for storing max proteins from "
                                                                    "MaxQuant protein groups with their resp. source "
                                                                    "(default: max_proteins_uniprot.db)")
    opt_args.add_argument("--level", "-l", action="store", required=False, nargs="?", metavar="String",
                          default="protein", type=str, help="Level on which counts need to be determined: protein, "
                                                            "peptide or PSM (default: protein).")
    opt_args.add_argument("--ens_db", "-e", action="store", required=False, nargs="?", metavar="PATH",
                          default="ENS_hsa_86.db", type=str, help="Ensembl database for transcript annotation "
                                                                  "(defaut: ENS_hsa_86.db)")

    args = parser.parse_args()

    #default workdir is CWD
    if args.workdir == "":
        args.workdir = os.getcwd()
    if args.score_file == "scores_uniprot.txt":
        args.score_file = args.workdir+"/"+args.level+"_"+args.score_file
    if args.venn_file == "venn_diagram_uniprot.png":
        args.venn_file = args.workdir+"/"+args.level+"_"+args.venn_file
    if args.max_protein_db == "max_proteins_uniprot.db":
        args.max_protein_db = args.workdir+"/"+args.max_protein_db
    if (args.level!="protein" and args.level!="peptide" and args.level!="PSM"):
        print "Level argument should be 'protein', 'peptide' or 'PSM'!"
        sys.exit()

    #List parameters
    print "Parameters:"
    for arg in vars(args):
        print '    %-15s\t%s' % (arg, getattr(args, arg))
    print
    sys.stdout.flush()

    ########
    # MAIN #
    ########

    # Parse accessions of combined fasta file
    accessions = parse_accessions(args.combined_fasta)

    #Init
    counts_per_source = defaultdict()
    conv_counts_per_source = defaultdict()

    #Protein group level
    if(args.level == "protein"):

        #Read in protein groups (except rev's and contaminants)
        protein_groups = read_proteingroups(args.protein_groups)

        #Get transcript info
        transcript_info = get_transcript_info(args.ens_db)

        #Count identified protein groups by max proteins
        (counts_per_source, max_protein_dict) = count_protein_groups(protein_groups, accessions, transcript_info)
        conv_counts_per_source = conv_counts(counts_per_source)

        #Create sqlite
        store_in_db(args.max_protein_db, max_protein_dict)

    #Peptide level
    elif(args.level == "peptide"):

        #read in peptides
        peptides = read_peptides(args.peptides)

        #Count identified peptides over sources
        counts_per_source = count_peptides(peptides, accessions)
        conv_counts_per_source = conv_counts(counts_per_source)

    #PSM level
    elif(args.level == "PSM"):

        #Read in PSMs
        psms = read_PSMs(args.psms)

        #Count identified peptides over sources
        counts_per_source = count_psms(psms, accessions)
        conv_counts_per_source = conv_counts(counts_per_source)

    # Create tab separated score file
    create_score_file(counts_per_source, conv_counts_per_source, args.score_file, args.level)

    # Create venn
    construct_venn(conv_counts_per_source, args.venn_file, args.level)



    #End of program message
    print
    print "-----------------------"
    print "[%s]: PROGRAM COMPLETE" % (convert_time(time.time()-starttime))
    print "-----------------------"
    print
    sys.stdout.flush()

    return

##########
#  SUBS  #
##########

def count_psms(psms, accessions):

    #Init
    counts_per_source = defaultdict()

    for psm in psms:
        #Filter contaminants
        m_con = re.search('CON\_\_', psms[psm]['Proteins'])
        if m_con:
            continue

        #Get all proteins in which the PSM is found
        found_in_proteins = re.split(';', psms[psm]['Proteins'])

        #Sources of proteins in which PSM is found
        found_sources_list = []
        for protein in found_in_proteins:
            sources_prot = re.split('\+', accessions[protein]['source'])
            found_sources_list.extend(sources_prot)

        #Get unique sources
        found_sources_list = list(set(found_sources_list))
        found_sources_list = sorted(found_sources_list, cmp=cmp_groups)
        found_sources = '+'.join(found_sources_list)

        #Add to counts
        if found_sources not in counts_per_source.keys():
            counts_per_source[found_sources] = 1
        else:
            counts_per_source[found_sources] += 1

    return counts_per_source

#Count source score per peptide
def count_peptides(peptides, accessions):

    #Init
    counts_per_source = defaultdict()

    for peptide in peptides:
        #Get all proteins in which the peptide is found
        found_in_proteins = re.split(';', peptides[peptide]['Proteins'])

        #Skip peptides which also map to a decoy or a contaminant
        m1 = re.search('CON\_\_', peptides[peptide]['Proteins'])
        m2 = re.search('REV\_\_', peptides[peptide]['Proteins'])
        if (m1 or m2):
            continue

        found_sources_list = []
        for protein in found_in_proteins:
            sources_prot = re.split('\+', accessions[protein]['source'])
            found_sources_list.extend(sources_prot)

        #Get unique sources
        found_sources_list = list(set(found_sources_list))
        found_sources_list = sorted(found_sources_list, cmp=cmp_groups)
        found_sources = '+'.join(found_sources_list)

        #Add +1 to the combination of sources found
        if found_sources not in counts_per_source.keys():
            counts_per_source[found_sources] = 1
        else:
            counts_per_source[found_sources]+=1

    return counts_per_source

#Construct sqlite db
def store_in_db(db, max_protein_dict):

    #Init
    table_name = "max_proteins"

    try:
        con = sqlite.connect(db)
    except:
        print "Could not connect to "+db
        sys.exit()

    with con:
        cur = con.cursor()

        # Delete existing table
        drop_query = "DROP TABLE IF EXISTS '" + table_name + "';"
        cur.execute(drop_query)

        # Create new table
        create_query = "CREATE TABLE IF NOT EXISTS '" + table_name + "' (" \
                                                                     "'protein_group' varchar(512) NOT NULL default '',"\
                                                                     "'sources' varchar(512) NOT NULL default '',"\
                                                                     "'max_proteins' varchar(512) NOT NULL default '');"
        cur.execute(create_query)

        #Insert into table
        for protein_group in max_protein_dict.keys():
            insert_query = "INSERT INTO '"+table_name+"' VALUES ('"+protein_group+"', '"+max_protein_dict[protein_group]['sources']+\
                           "', '"+max_protein_dict[protein_group]['max_proteins']+"');"
            #print insert_query
            cur.execute(insert_query)

    return

#Construct Venn diagram of counts per bincode
def construct_venn(counts, venn_file, level):

    #2group venn
    if(len(counts.keys())==3):

        subsets = (counts['Proteoformer'], counts['Uniprot'], counts['Uniprot+Proteoformer'])
        set_labels = ('1', '2')
        textstr = "1: proteoformer\n2: uniprot"

        fig, ax = plt.subplots(nrows=1, ncols=1)
        v = venn2(subsets=subsets, set_labels=set_labels, ax=ax)

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(-0.03, 0.07, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
        plt.tight_layout()

        fig.savefig(venn_file)
        plt.close(fig)

    return

#Create score file
def create_score_file(counts_per_source, conv_counts_per_source, score_file, level):

    with open(score_file, 'w') as FW:
        FW.write("Counts over different sources on "+level+" level after MaxQuant analysis:\n")
        line = '{:30} {:20}'.format("Files", "Counts")
        FW.write(line+"\n")
        for sources in counts_per_source.keys():
            line='{:30} {:20}'.format(sources, str(counts_per_source[sources]))
            FW.write(line+"\n")

        FW.write("\n")
        line = '{:30} {:20}'.format("Files", "Counts")
        FW.write(line + "\n")
        for sources in conv_counts_per_source.keys():
            line = '{:30} {:20}'.format(sources, str(conv_counts_per_source[sources]))
            FW.write(line + "\n")


    return

#Convert Uniprot count labels
def conv_counts(counts_per_source):

    #Init
    conv_counts_per_source = defaultdict()

    #Convert source labels (TrEMBL and SwissProt should be Uniprot)
    for sources in counts_per_source.keys():
        m1 = re.search('TrEMBL\+SwissProt(.*)', sources)
        m2 = re.search('SwissProt(.*)', sources)
        m3 = re.search('TrEMBL(.*)', sources)
        if m1:
            if 'Uniprot'+m1.group(1) in conv_counts_per_source.keys():
                conv_counts_per_source['Uniprot'+m1.group(1)] += counts_per_source[sources]
            else:
                conv_counts_per_source['Uniprot'+m1.group(1)] = counts_per_source[sources]
        elif m2:
            if 'Uniprot'+m2.group(1) in conv_counts_per_source.keys():
                conv_counts_per_source['Uniprot'+m2.group(1)] += counts_per_source[sources]
            else:
                conv_counts_per_source['Uniprot'+m2.group(1)] = counts_per_source[sources]
        elif m3:
            if 'Uniprot'+m3.group(1) in conv_counts_per_source.keys():
                conv_counts_per_source['Uniprot'+m3.group(1)] += counts_per_source[sources]
            else:
                conv_counts_per_source['Uniprot'+m3.group(1)] = counts_per_source[sources]
        else:
            if sources in conv_counts_per_source.keys():
                conv_counts_per_source[sources] += counts_per_source[sources]
            else:
                conv_counts_per_source[sources] = counts_per_source[sources]

    return conv_counts_per_source

#Construct counts per source and save the maximum scoring proteins per protein group
def count_protein_groups(protein_groups, accessions, transcript_info):

    #Init
    counts_per_source = defaultdict()
    max_protein_dict = defaultdict(lambda: defaultdict())

    #Go over all identified protein groups
    for protein_group in protein_groups:
        #Get all proteins and their (all) counts
        all_proteins = re.split(';', protein_group)
        counts = re.split(';', protein_groups[protein_group]['Peptide counts (all)'])
        protein_counts = defaultdict()
        for i in range(0, len(counts)):
            protein_counts[all_proteins[i]] = int(counts[i])

        #Select proteins with max counts
        max_count = max(protein_counts.values())
        max_proteins = []
        for protein in all_proteins:
            if (protein_counts[protein]==max_count):
                max_proteins.append(protein)

        #Get sources of all max proteins and make a combined source label
        sources_list = []
        for protein in max_proteins:
            m1 = re.search('CON\_\_', protein)
            m2 = re.search('REV\_\_', protein)
            if (m1 or m2):
                continue
            to_add_sources = re.split('\+', accessions[protein]['source'])
            sources_list.extend(to_add_sources)
            transcript_id = accessions[protein]['tr_stable_id']

        #Get unique sources
        sources_list = list(set(sources_list))
        sources_list = sorted(sources_list, cmp=cmp_groups)
        sources = '+'.join(sources_list)

        #Add +1 to the combination of sources found
        if sources not in counts_per_source.keys():
            counts_per_source[sources] = 1
        else:
            counts_per_source[sources]+=1

        #Add max protein info to dict
        max_protein_dict[protein_group]['max_proteins'] = '|'.join(max_proteins)
        max_protein_dict[protein_group]['sources'] = sources

    return (counts_per_source, max_protein_dict)

#Cmp for sorting data groups
def cmp_groups(x,y):
    if x=='TrEMBL' or x=="SwissProt":
        return -1
    elif y=='TrEMBL' or y=='SwissProt':
        return 1
    else:
        if x<y:
            return -1
        else:
            return 1

#Parse accessions of combined fasta file
def parse_accessions(combined_fasta):

    #Init
    #accessions->{main_acc}->{feature} = value
    accessions = defaultdict(lambda: defaultdict())

    with open(combined_fasta, 'r') as FR:
        for line in FR:
            #Get only the header lines, not the sequences
            m1 = re.search('^>(.+?)\|(.+?)\|(.+)$', line)
            if m1:
                kind = m1.group(1)
                main_acc = m1.group(2)
                descr = m1.group(3)
                part = main_acc+"|"+descr

                if kind=="sp":
                    accessions[main_acc]["source"] = "SwissProt"
                elif kind=="tr":
                    accessions[main_acc]["source"] = "TrEMBL"
                elif kind=="generic":
                    accessions[main_acc]["source"] = "Proteoformer"

                #For Proteoformer generated main accessions
                if accessions[main_acc]["source"] == "Proteoformer":
                    m2= re.search('^(ENST\d+?)\_(\S+?)\|(ENST\d+?)',part)
                    if m2:
                        accessions[main_acc]['tr_stable_id'] = m2.group(1)
                        accessions[main_acc]['protein_name'] = ""
                        accessions[main_acc]['description'] = ""
                #For Uniprot main accessions
                else:
                    m3 = re.search('^.+?\|(.+)OS=', part)
                    if m3:
                        accessions[main_acc]['protein_name'] = main_acc
                        accessions[main_acc]['description'] = m3.group(1).replace(',', '')
                        accessions[main_acc]['tr_stable_id'] = ''
                        #Parse side accessions
                        m4 = re.search('OS=.+?\[(.+?)\]', descr)
                        if m4:
                            accessions[main_acc]['side_accessions'] = m4.group(1)
                            #Check if there is PROTEOFORMER in the side accessions
                            m5 = re.search('ENST', accessions[main_acc]['side_accessions'])
                            if m5:
                                accessions[main_acc]["source"] = accessions[main_acc]["source"] + "+Proteoformer"
                            #Search for bincode of the first accessions
                            m6 = re.search('^(ENST\d+?)\_(\S+?)\|(ENST\d+?)', accessions[main_acc]['side_accessions'])
                            if m6:
                                accessions[main_acc]['tr_stable_id'] = m6.group(1)
                        else:
                            accessions[main_acc]['side_accessions'] = ""

    return accessions

#Read MS MS spectra from maxquant
def read_PSMs(in_file):

    #Init
    psms = defaultdict(lambda: defaultdict())

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

        #Parse data
        for line_nr in range(1, len(lines)):
            data = re.split('\t', lines[line_nr])
            #Remove reverse identifications
            if(data[idx_reverses]!='+'):
                for val_nr in range(0, len(data)):
                    psms[data[52]][col_headers[val_nr]] = data[val_nr]

    return psms

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

# Get transcript-gene info
def get_transcript_info(ens_db):

    #Init
    transcript_info = defaultdict(lambda: defaultdict())

    #Connect to db
    try:
        con = sqlite.connect(ens_db)
    except:
        print "Could not connect to "+ens_db
        sys.exit()

    with con:
        cur = con.cursor()

        query = "SELECT tr.stable_id, g.stable_id, g.description FROM transcript AS tr JOIN gene AS g ON tr.gene_id=g.gene_id;"
        cur.execute(query)

        output = cur.fetchall()

        for i in range(0, len(output)):
            transcript_info[output[i][0]]['gene_id'] = output[i][1]
            transcript_info[output[i][0]]['description'] = output[i][2]

    return transcript_info


#Time converter
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
