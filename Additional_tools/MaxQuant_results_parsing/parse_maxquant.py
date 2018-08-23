#!/usr/bin/env python
import traceback
import sys
import os
import time
import argparse
from collections import defaultdict
import re
#from dict_functions import dict_funcs
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import sqlite3 as sqlite

__author__ = 'Steven Verbruggen'
#Execute 'python parse_maxquant.py -h' for more information

def main():

    starttime = time.time()

    print
    print "#########################"
    print "# Parse MaxQuant output #"
    print "#########################"
    print
    print "This program is part of the PROTEOFORMER pipeline"
    print

    #Parse arguments from command line
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS,
                                     description="This tool is part of the PROTEOFORMER pipeline. It parses the"
                                                 "output of a MaxQuant MS validation.",
                                     add_help=False)

    man_args = parser.add_argument_group("Mandatory parameters")
    man_args.add_argument("--protein_groups", "-g", action="store", required=True, nargs="?", metavar="PATH",
                          default="", type=str, help="Input file with all MaxQuant identified protein groups (mandatory)")
    man_args.add_argument("--peptides", "-p", action="store", required=True, nargs="?", metavar="PATH",
                           default="", type=str, help="Input file with all MaxQuant identified peptides (mandatory)")
    man_args.add_argument("--psms", "-x", action="store", required=True, nargs="?", metavar="PATH",
                          default="", type=str, help="Input file with all MaxQuant MS/MS spectra (mandatory)")
    man_args.add_argument("--overview_file", "-t", action="store", required=True, nargs="?", metavar="PATH",
                          default="", type=str, help="Overview file coming out of db combination tool (mandatory)")
    man_args.add_argument("--combined_fasta", "-f", action="store", required=True, nargs="?", metavar="PATH",
                          default="", type=str, help="Combined fasta file on which MaxQuant was run. Underlying "
                                                     "fasta files should be generated with verbose output (mandatory)")
    man_args.add_argument("--ens_db", "-e", action="store", required=True, nargs="?", metavar="PATH",
                          default="", type=str, help="Ensembl databse for transcript and gene info (mandatory)")


    opt_args = parser.add_argument_group("Optional parameters")
    opt_args.add_argument("--help", "-h", action="help", help="Show help message and exit")
    opt_args.add_argument("--workdir", "-w", action="store", required=False, nargs="?", metavar="FOLDER", default="",
                        type=str, help="Working directory (default: CWD)")
    opt_args.add_argument("--score_file", "-s", action="store", required=False, nargs="?", metavar="PATH",
                          default="scores_per_bincode.txt", type=str, help="Output tab separated score file"
                                                                              " (default: level_scores_per_bincode.txt)")
    opt_args.add_argument("--venn_file", "-v", action="store", required=False, nargs="?", metavar="PATH",
                          default="venn_diagram.png", type=str, help="Path to a png file where the Venn diagram will "
                                                                     "be stored (default: level_venn_diagram.png)")
    opt_args.add_argument("--removed_redundancy", "-r", action="store", required=False, nargs="?", metavar="Y/N",
                          default="Y", type=str, help="Whether the redundancy was removed in uncombined fasta files. "
                                                      "Should be Y or N (default: Y).")
    opt_args.add_argument("--max_protein_db", "-d", action="store", required=False, nargs="?", metavar="PATH",
                          default="max_proteins.db", type=str, help="SQLite db path for storing max proteins from "
                                                                    "MaxQuant protein groups with their resp. bin codes "
                                                                    "(default: max_proteins.db)")
    opt_args.add_argument("--level", "-l", action="store", required=False, nargs="?", metavar="String",
                          default="protein", type=str, help="Level on which counts need to be determined: protein, "
                                                            "peptide or PSM (default: protein).")

    args = parser.parse_args()

    #default workdir is CWD
    if args.workdir == "":
        args.workdir = os.getcwd()
    if args.score_file == "ms_scores_per_bincode.txt":
        args.score_file = args.workdir+"/"+args.score_file
    if args.venn_file == "ms_venn_diagram.png":
        args.venn_file = args.workdir+"/"+args.venn_file
    if args.max_protein_db == "max_proteins.db":
        args.max_protein_db = args.workdir+"/"+args.max_protein_db
    if (args.removed_redundancy!="Y" and args.removed_redundancy!="N"):
        print "Removed redundancy argument should be Y or N!"
        sys.exit()
    if (args.level!="protein" and args.level!="peptide" and args.level!="PSM"):
        print "Level argument should be 'protein', 'peptide' or 'PSM'!"
        sys.exit()
    if (args.venn_file=="venn_diagram.png"):
        args.venn_file = args.level+"_"+args.venn_file
    if (args.score_file=="scores_per_bincode.txt"):
        args.score_file = args.level+"_"+args.score_file

    #List parameters
    print "Parameters:"
    for arg in vars(args):
        print '    %-15s\t%s' % (arg, getattr(args, arg))
    print
    sys.stdout.flush()

    ########
    # MAIN #
    ########

    #Protein group level
    if (args.level=="protein"):
        #Read in protein groups (except rev's and contaminants)
        protein_groups = read_proteingroups(args.protein_groups)

        #Read in peptide file
        peptides = read_peptides(args.peptides)

        #Get files and bincodes from db combination overview file
        (files, bincodes) = parse_overview_file(args.overview_file)
        length_bin_code = len(files)

        #Parse accessions in the fasta file
        accessions = parse_accessions(args.combined_fasta)

        #Get transcript-gene info
        transcript_info = get_transcript_info(args.ens_db)

        #Count identified protein groups per bincode
        (counts_per_bincode, max_protein_dict) = count_bincodes_proteins(protein_groups, peptides, accessions, bincodes, length_bin_code, args.removed_redundancy, transcript_info)

        #Create file with score overview
        create_score_file(counts_per_bincode, args.score_file, args.level)

        #Create venn diagram
        construct_venn(counts_per_bincode, args.venn_file, files, args.level)

        #Store max proteins in sqlite db
        store_in_db(args.max_protein_db, max_protein_dict)

    #Peptide level
    elif(args.level=="peptide"):
        #Read peptides
        peptides = read_peptides(args.peptides)

        # Get files and bincodes from db combination overview file
        (files, bincodes) = parse_overview_file(args.overview_file)
        length_bin_code = len(files)

        # Count identified peptides per bincode
        counts_per_bincode = count_bincodes_peptides(peptides, bincodes, length_bin_code)

        #create file with score overview
        create_score_file(counts_per_bincode, args.score_file, args.level)

        #Construct venn
        construct_venn(counts_per_bincode, args.venn_file, files, args.level)

    #MS/MS level
    elif(args.level=="PSM"):
        #Read MS/MS table
        PSMs = read_PSMs(args.psms)

        # Get files and bincodes from db combination overview file
        (files, bincodes) = parse_overview_file(args.overview_file)
        length_bin_code = len(files)

        # Count identified psms per bincode
        counts_per_bincode = count_bincodes_psms(PSMs, bincodes, length_bin_code)

        # create file with score overview
        create_score_file(counts_per_bincode, args.score_file, args.level)

        # Construct venn
        construct_venn(counts_per_bincode, args.venn_file, files, args.level)

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

#Count psms over bincodes
def count_bincodes_psms(psms, bincodes, length_bin_code):

    #Init
    counts_per_bincode = defaultdict()
    for bincode in bincodes.values():
        counts_per_bincode[bincode] = 0

    for psm in psms:
        #Filter empty proteins
        if(psms[psm]['Proteins']!=''):
            #Search for contaminants
            m_con = re.search('CON', psms[psm]['Proteins'])
            if not m_con:
                # Get all proteins in which the psm is found
                found_in_proteins = re.split(';', psms[psm]['Proteins'])

                identified_bin_codes = []
                for protein in found_in_proteins:
                    #Get bin code of all proteins
                    m = re.search('^ENST\d+\_\S+\_\d+\_\S+\_(\d+)db\d+$', protein)
                    if m:
                        identified_bin_codes.append(m.group(1))

                #make a consensus bin code for all proteins in which the psm is found
                score_bin_code = '0' * length_bin_code
                for bincode in identified_bin_codes:
                    new_bin_code = ''
                    for i in range(0, length_bin_code):
                        if(score_bin_code[i]=='1' or bincode[i]=='1'):
                            new_bin_code+='1'
                        else:
                            new_bin_code+='0'
                    score_bin_code = new_bin_code
                # +1 for combined bin code
                counts_per_bincode[score_bin_code] += 1

    return counts_per_bincode

#Count peptides over bincodes
def count_bincodes_peptides(peptides, bincodes, length_bin_code):

    #Init
    counts_per_bincode = defaultdict()
    for bincode in bincodes.values():
        counts_per_bincode[bincode] = 0

    for peptide in peptides:
        #Get all proteins where peptide is found
        found_in_proteins = re.split(';', peptides[peptide]['Proteins'])

        #Get bincode for each protein
        identified_bincodes = []
        for protein in found_in_proteins:
            m = re.search('^ENST\d+\_\S+\_\d+\_\S+\_(\d+)db\d+$', protein)
            if m:
                identified_bincodes.append(m.group(1))

        #Make a consensus bin code for all proteins in which the peptide is found
        score_bin_code = '0' * length_bin_code
        for bincode in identified_bincodes:
            new_bin_code = ''
            for i in range(0, length_bin_code):
                if(score_bin_code[i]=='1' or bincode[i]=='1'):
                    new_bin_code+='1'
                else:
                    new_bin_code+='0'
            score_bin_code = new_bin_code

        #+1 for combined bincode
        counts_per_bincode[score_bin_code] += 1

    return counts_per_bincode

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
                                                                     "'max_proteins' varchar(512) NOT NULL default '',"\
                                                                     "'annotations' varchar(128) NOT NULL default ''," \
                                                                     "'bin_code' varchar(128) NOT NULL default ''," \
                                                                     "'gene_ids' varchar(128) NOT NULL default ''," \
                                                                     "'descriptions' varchar(512) NOT NULL default '');"
        cur.execute(create_query)

        #Insert into table
        for protein_group in max_protein_dict.keys():
            insert_query = "INSERT INTO '"+table_name+"' VALUES ('"+protein_group+"', '"+max_protein_dict[protein_group]['max_proteins']+\
                           "', '"+max_protein_dict[protein_group]['annotations']+"', '"+max_protein_dict[protein_group]['bin_code']+\
                            "', '"+max_protein_dict[protein_group]['genes']+"', '"+max_protein_dict[protein_group]['descriptions']+"');"
            #print insert_query
            cur.execute(insert_query)

    return

#Construct Venn diagram of counts per bincode
def construct_venn(counts, venn_file, files, level):

    #2group venn
    if(len(counts.keys())==3):
        subsets = (counts['10'], counts['01'], counts['11'])
        set_labels = ()
        for i in sorted(files.keys()):
            set_labels+=(files[i])
        textstr = ""
        for i in sorted(files.keys()):
            textstr = textstr+str(i)+": "+files[i]+"\n"

        fig, ax = plt.subplots(figsize=(10,6), nrows=1, ncols=1)
        v = venn2(subsets=subsets, set_labels=set_labels, ax=ax)

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(-0.3, 0.15, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax.set_title(level+" counts over the different bin codes")
        plt.tight_layout()

        fig.savefig(venn_file)
        plt.close(fig)

    #3group venn
    if(len(counts.keys())==7):
        subsets = (counts['100'], counts['010'], counts['110'], counts['001'], counts['101'], counts['011'], counts['111'])
        set_labels = ()
        for i in sorted(files.keys()):
            set_labels += (files[i],)
        textstr = ""
        for i in sorted(files.keys()):
            textstr = textstr+str(i)+": "+files[i]+"\n"

        fig, ax = plt.subplots(figsize=(10,6), nrows=1, ncols=1)
        v = venn3(subsets=subsets, set_labels=set_labels, ax=ax)

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(-0.3, 0.15, textstr, transform=ax.transAxes, fontsize=12, verticalalignment='top', bbox=props)
        ax.set_title(level + " counts over the different bin codes")
        plt.tight_layout()

        fig.savefig(venn_file)
        plt.close(fig)

    return

#Create score file
def create_score_file(counts_per_bincode, score_file, level):

    #Open file
    with open(score_file, 'w') as FW:
        print "Counts of different bin codes on "+level+" level after MaxQuant analysis:"
        FW.write("Counts of different bin codes on "+level+" level after MaxQuant analysis:\n")
        print "\t\tBinCode\t\t"+level+" count"
        FW.write("\t\tBinCode\t\t"+level+" count\n")

        #Sort bincodes and print/write counts
        number=1
        for i in sorted(counts_per_bincode.keys(), cmp=cmp_files):
            print str(number)+":\t\t"+str(i)+"\t\t"+str(counts_per_bincode[i])
            FW.write(str(number) + ":\t\t" + str(i) + "\t\t" + str(counts_per_bincode[i])+"\n")
            number+=1

    return

# Cmp for sorting files
def cmp_files(x, y):
    if x.count('1') > y.count('1'):
        return 1
    elif x.count('1') < y.count('1'):
        return -1
    else:
        if x < y:
            return 1
        else:
            return -1

#Count bincodes over different identified protein groups
def count_bincodes_proteins(protein_groups, peptides, accessions, bincodes, length_bin_code, removed_redundancy, transcript_info):

    #Init
    counts_per_bincode=defaultdict()
    for bincode in bincodes.values():
        counts_per_bincode[bincode] = 0
    max_protein_dict = defaultdict(lambda:defaultdict())


    for protein_group in protein_groups:
    #for protein_group in ['ENST00000217426_20_34303270_aTIS_001db3;ENST00000217426_20_34303285_5UTR_100db1;ENST00000217426_20_34303288_5UTR_010db2;ENST00000430630_9_120720673_ntr_100db1']:

        #Get all proteins of the group and their counts
        all_proteins = re.split(';', protein_group)
        counts = re.split(';', protein_groups[protein_group]['Peptide counts (all)'])
        protein_counts = defaultdict()
        for i in range(0,len(all_proteins)):
            protein_counts[all_proteins[i]] = int(counts[i])

        ## Correct for methionine-clipped peptides which are not accounted to proteins with extensions as their main accession but with the shorter form underneath
        #Only necessary when the redundancy is removed in the uncombined fasta files
        if(removed_redundancy=="Y"):

            #Get protein with highest count of this protein group
            m_all_in = re.search('^(\S+?);', protein_group)
            if m_all_in:
                all_in_protein = m_all_in.group(1)
            else:
                #Only one protein in the protein group
                all_in_protein = protein_group

            #Get all peptides
            all_peptide_ids = re.split(';',protein_groups[protein_group]['Peptide IDs'])

            #Go over all peptides
            for peptide_id in all_peptide_ids:

                #Only for peptides which start on position 2
                if(int(peptides[peptide_id]['Start position'])==2):

                    #Get the proteins in which the peptide is included
                    proteins_with_the_peptide = re.split(';', peptides[peptide_id]['Proteins'])

                    #For proteins of the protein group which not include the peptide
                    not_incl_proteins = list(set(all_proteins) - set(proteins_with_the_peptide))
                    for not_incl_protein in not_incl_proteins:

                        #Check if the start of the all-in protein (=genomic coordinate of the start position) is in a side accession of that protein
                        for side_acc in accessions[not_incl_protein]['side_accessions']:

                            #Has to be in the same db
                            if(accessions[not_incl_protein]['main_db']==accessions[not_incl_protein]['side_accessions'][side_acc]['side_db']):

                                #Check for presence of the shorter form
                                if (accessions[not_incl_protein]['side_accessions'][side_acc]['start']==accessions[all_in_protein]['start']):
                                    #Then do the correction for that specific peptide
                                    protein_counts[not_incl_protein] += 1

        #Select all proteins with max count
        max_count = max(protein_counts.values())
        max_proteins = []
        annotations = []
        transcripts = []
        genes=[]
        descriptions=[]
        for protein in protein_counts.keys():
            if(protein_counts[protein]==max_count):
                max_proteins.append(protein)
                #ENST00000357849_8_67061900_CDS_010db2
                m = re.search('^(ENST\d+?)\_\S+?\_\d+?\_(\S+?)\_\d+?db\d$', protein)
                if m:
                    if m.group(2) not in annotations:
                        annotations.append(m.group(2))
                    if m.group(1) not in transcripts:
                        transcripts.append(m.group(1))
        for transcript in transcripts:
            if transcript_info[transcript]['gene_id'] not in genes:
                genes.append(transcript_info[transcript]['gene_id'])
                descriptions.append(re.sub('\'','',transcript_info[transcript]['description']))
        max_protein_dict[protein_group]['max_proteins'] = '|'.join(max_proteins)
        max_protein_dict[protein_group]['annotations'] = '|'.join(sorted(annotations))
        max_protein_dict[protein_group]['genes'] = '|'.join(genes)
        max_protein_dict[protein_group]['descriptions'] = '|'.join(descriptions)

        #Make a consensus bin code of all max proteins
        score_bin_code = '0' * length_bin_code
        for protein in max_proteins:
            new_bin_code = ''
            for i in range(0, length_bin_code):
                if(score_bin_code[i]=='1' or accessions[protein]['bin_code'][i]=='1'):
                    new_bin_code+='1'
                else:
                    new_bin_code+='0'
            score_bin_code = new_bin_code

        #Give that bincode a +1 in score
        counts_per_bincode[score_bin_code]+=1
        max_protein_dict[protein_group]['bin_code'] = score_bin_code

    return (counts_per_bincode, max_protein_dict)

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

#Parse accessions of the combined fasta file
def parse_accessions(combined_fasta):

    #Init
    accessions = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict())))
    count=0

    #Read fasta file
    with open(combined_fasta, 'r') as FR:
        for line in FR:
            count+=1
            #only accessions, not the sequences
            m_acc = re.search('^>generic\|(ENST\S+?)\|ENSG\d+', line)
            if m_acc:
                main_acc = m_acc.group(1)
                #Get tr_stable_id and start
                m_start = re.search('^(ENST\d+)_\S+?_(\d+?)_\S+?_(\d+?)db(\d+)$', main_acc)
                if m_start:
                    accessions[main_acc]['tr_stable_id'] = m_start.group(1)
                    accessions[main_acc]['start'] = int(m_start.group(2))
                    accessions[main_acc]['bin_code'] = m_start.group(3)
                    accessions[main_acc]['main_db'] = int(m_start.group(4))
                #Get all info of side accessions
                m_side = re.search('\[(.+)\]$', line)
                if m_side:
                    side_accs = re.split('#', m_side.group(1))
                    for side_acc in side_accs:
                        m_side_acc = re.search('^(ENST\d+?)_\S+?_(\d+?)_\S+?_db(\d+)$', side_acc)
                        if m_side_acc:
                            accessions[main_acc]['side_accessions'][side_acc]['tr_stable_id'] = m_side_acc.group(1)
                            accessions[main_acc]['side_accessions'][side_acc]['start'] = int(m_side_acc.group(2))
                            accessions[main_acc]['side_accessions'][side_acc]['side_db'] = int(m_side_acc.group(3))
                else:
                    accessions[main_acc]['side_accessions'] = []

    return accessions

#Get files and bincodes from overview file
def parse_overview_file(file):

    #Init
    files = defaultdict()
    bincodes = defaultdict()

    #Read file
    with open(file, 'r') as FR:
        for line in FR:
            #Strip trailing characters
            line = line.rstrip("\n")
            line = line.rstrip("\r")

            #Parts of file where bincodes are listed
            m_bincode = re.search('^(\d+):\t\t([0,1]+)\t\t', line)
            if m_bincode:
                bincodes[m_bincode.group(1)] = m_bincode.group(2)
            else:
                m_file = re.search('^(\d+):\t\t(\S+)$', line)
                if m_file:
                    files[m_file.group(1)] = m_file.group(2)

    return (files, bincodes)

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

        #Parse data
        for line_nr in range(1, len(lines)):
            data = re.split('\t', lines[line_nr])
            #Remove reverse identifications
            if(data[46]!='+'):
                for val_nr in range(0, len(data)):
                    psms[data[50]][col_headers[val_nr]] = data[val_nr]

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

        #Parse data
        for line_nr in range(1,len(lines)):
            data = re.split('\t', lines[line_nr])
            #Remove contaminants and reverse identifications
            if(data[44]!='+' and data[45]!='+'):
                for val_nr in range(0, len(data)):
                    peptides[data[46]][col_headers[val_nr]] = data[val_nr] #data[46] is peptide ID

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

        #Parse data
        for line_nr in range(1,len(lines)):
            data = re.split('\t', lines[line_nr])
            #Remove contaminants and reverse identifications
            if(data[24]!='+' and data[25]!='+'):
                for val_nr in range(0, len(data)):
                    protein_groups[data[0]][col_headers[val_nr]] = data[val_nr]

    return protein_groups


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