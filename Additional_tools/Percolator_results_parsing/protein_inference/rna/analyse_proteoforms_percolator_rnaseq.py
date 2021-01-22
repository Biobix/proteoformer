import traceback
import time
import argparse
import os
import sys
#import numpy as np
#import matplotlib
#import matplotlib.pyplot as plt
#from matplotlib_venn import venn2, venn3
#import seaborn as sns
import pandas as pd
from collections import defaultdict
import re
#from statistics import mean
#sns.set_style('whitegrid')

__author__ = 'Steven Verbruggen'
#Execute 'python analyse_proteoforms_percolator.py -h' for more

def main():

    starttime = time.time()

    print
    print "####################################"
    print "# Analyze proteoforms (Percolator) #"
    print "####################################"
    print
    print "This program is part of the PROTEOFORMER pipeline"
    print


    #Parse arguments from command line
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS,
                                     description="This tool is part of the PROTEOFORMER pipeline. It analyses "
                                                 "the different proteoforms found after MS analysis of a combined database "
                                                 "of PROTEOFORMER and UniProt. It searches for the classification of all "
                                                 "proteoforms found by PROTEOFORMER, not yet in UniProt and with confiramtion"
                                                 " of MS with Percolator after MaxQuant.",
                                     add_help=False)

    man_args = parser.add_argument_group("Mandatory parameters")
    man_args.add_argument("--protein_groups", "-g", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="MaxQuant output file of protein groups (mandatory)")
    man_args.add_argument("--peptides", "-p", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="MaxQuant output file of peptides(mandatory)")
    man_args.add_argument("--fasta", "-f", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="Path to the combined fasta file with decoys and contaminants included (mandatory)")

    opt_args = parser.add_argument_group("Optional parameters")
    opt_args.add_argument("--help", "-h", action="help", help="Show help message and exit")
    opt_args.add_argument("--workdir", "-w", action="store", required=False, nargs="?", metavar="FOLDER", default="",
                        type=str, help="Working directory (default: CWD)")
    opt_args.add_argument("--qval_threshold_pept", "-q", action="store", required=False, nargs="?", metavar="FLOAT",
                          default=0.01, type=float, help="Q-value threshold on which the peptides will be filtered (default: 0.01).")
    opt_args.add_argument("--PEP_threshold_pept", "-P", action="store", required=False, nargs="?", metavar="FLOAT",
                          default=0.05, type=float, help="Peptide PEP threshold on which the peptides will be filtered (default: 0.05)")
    opt_args.add_argument("--csv_output", "-x", action="store", required=False, nargs="?", metavar="PATH",
                          type=str, default="novel_proteoforms_percolator.csv",
                          help="CSV output file of classifications (default: novel_proteoforms_percolator.csv)")

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

    #Load peptides file
    peptides = pd.read_csv(args.peptides, sep="\t", header=0)
    peptides.loc[:, 'stripped_seqs'] = [elem.strip('_.') for elem in list(peptides.loc[:, 'peptide'])]

    #Load protein groups file
    protein_groups = pd.read_csv(args.protein_groups, sep="\t", header=0)

    ##Load fasta file
    #fasta = load_fasta(args.fasta)

    #Select significant peptides
    sign_peptides = select_significant_peptides(peptides, args.qval_threshold_pept, args.PEP_threshold_pept)

    #Select CustomDB-unique peptides
    customdb_peptides = select_customdb_peptides(sign_peptides)
    print("Amount of custom db unique significant peptides: "+str(len(customdb_peptides)))

    #Analyze identifications
    identifications = determine_identifications(customdb_peptides, protein_groups)

    #Output identifications
    output_csv(identifications, args.csv_output)


    #Analyze the type of identifications




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
        keys = ["base_proteoform","proving_peptides","peptide_PEP","protein_group"]
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

def determine_identifications(customdb_peptides, protein_groups):

    #init
    identifications = defaultdict(lambda: defaultdict())

    #Go over all custom DB unique peptide events
    for idx, row in customdb_peptides.iterrows():

        #Peptide sequence
        pepseq = row['stripped_seqs']

        #Matching protein group
        proteinGroupId = ""
        possibleProteinGroups = protein_groups[[pepseq in elem for elem in protein_groups.loc[:,'peptideIds']]]
        #However, protein groups from which the pep can be a subsequence are also present. Therefore, filter further on full sequence match
        for idx_pg,row_pg in possibleProteinGroups.iterrows():
            pg_peptides = row_pg['peptideIds'].split(' ')
            if pepseq in pg_peptides:
                proteinGroupId = str(row_pg['ProteinId'])

        #Already existing protein groups are thus taken together in following lines
        if proteinGroupId not in identifications.keys():
            identifications[proteinGroupId]['protein_group'] = proteinGroupId
            identifications[proteinGroupId]['proving_peptides'] = []
            identifications[proteinGroupId]['peptide_PEP'] = []
        identifications[proteinGroupId]['proving_peptides'].append(pepseq)
        identifications[proteinGroupId]['peptide_PEP'].append(str(row['posterior_error_prob']))
        identifications[proteinGroupId]['base_proteoform'] = proteinGroupId.split(';')[0]

    return identifications

def load_fasta(in_fasta):

    #Init
    input_seqs = defaultdict(lambda: defaultdict())

    with open(in_fasta, 'r') as FR:
        lines = FR.readlines()
        accession = ''
        protein_name = ''

        for i in range(0, len(lines)):
            lines[i] = lines[i].rstrip('\n')
            m = re.search('^>.+$', lines[i])
            if m:
                accession = lines[i]
                m2 = re.search('^>.+\|(.+)\|.+$', accession)
                if m2:
                    protein_name = m2.group(1)
            else:
                seq = lines[i]
                input_seqs[protein_name]['accession']=accession
                input_seqs[protein_name]['seq'] = seq
                accession = ''
                protein_name = ''

    return input_seqs

def select_customdb_peptides(peptides):

    peps_proteinids = list(peptides['proteinIds'])
    mask = []
    for proteinid in peps_proteinids:
        pep_proteinids = proteinid.split(';')
        m = re.compile('^ENST')
        filtered_pep_proteinids = list(filter(m.match, pep_proteinids))
        if len(pep_proteinids) == len(filtered_pep_proteinids):
            mask.append(True)
        else:
            mask.append(False)
    customdb_peptides = peptides[mask]

    return customdb_peptides

def select_significant_peptides(peptides, qval_th, PEP_th):
    qval_filt_peptides = peptides[peptides['q-value'] < qval_th]
    qval_PEP_filt_peptides = qval_filt_peptides[qval_filt_peptides['posterior_error_prob']<PEP_th]

    return qval_PEP_filt_peptides

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