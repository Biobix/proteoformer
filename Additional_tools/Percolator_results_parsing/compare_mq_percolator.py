import traceback
import time
import sys
import argparse
import numpy as np
import pandas as pd
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3

__author__ = 'Steven Verbruggen'
#Execute 'python compare_mq_percolator.py -h' for more information

def main():

    starttime = time.time()

    print
    print "#########################################"
    print "# Compare MQ and Percolator proteoforms #"
    print "#########################################"
    print
    print "This program is part of the PROTEOFORMER pipeline"
    print

    #Parse arguments from command line
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS,
                                     description="This tool is part of the PROTEOFORMER pipeline. It compares "
                                                 "the different proteoforms found after MS analysis of a combined database "
                                                 "of PROTEOFORMER and UniProt. It compares the results between a MQ and "
                                                 "Percolator search.",
                                     add_help=False)


    man_args = parser.add_argument_group("Mandatory parameters")
    man_args.add_argument("--mq_csv", "-m", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="MaxQuant output file of proteoform classifications (mandatory)")
    man_args.add_argument("--percolator_csv", "-p", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="Percolator output file of proteoform classifications (mandatory)")

    opt_args = parser.add_argument_group("Optional parameters")
    opt_args.add_argument("--help", "-h", action="help", help="Show help message and exit")
    opt_args.add_argument("--venn_file", "-v", action="store", required=False, nargs="?", metavar="PATH",
                          type=str, default="overlap_mq_percolator.png", help="Venn diagram output file of overlap (default: overlap_mq_percolator.png)")

    args = parser.parse_args()

    #List parameters
    print "Parameters:"
    for arg in vars(args):
        print '    %-15s\t%s' % (arg, getattr(args, arg))
    print
    sys.stdout.flush()

    ########
    # MAIN #
    ########

    #Load maxquant proteoform classifications results
    mq_proteoforms = load_csv_file(args.mq_csv)

    #Load Percolator proteoforms
    percolator_proteoforms = load_csv_file(args.percolator_csv)

    #Check the overlap between MQ and Percolator
    proteoform_overlap, mq_pg_ids, percolator_pg_ids = check_overlap(mq_proteoforms, percolator_proteoforms)
    #print_dict(proteoform_overlap)

    print("\n\nOverlap: " + str(len(proteoform_overlap.keys())))
    print("Maxquant protein group IDs in overlap: "+str(len(set(mq_pg_ids))))
    print("Percolator protein group IDs in overlap: "+str(len(set(percolator_pg_ids))))

    #Construct counts
    counts = defaultdict()
    counts['maxquant'] = len(mq_proteoforms) - len(set(mq_pg_ids))
    counts['percolator'] = len(percolator_proteoforms) - len(set(percolator_pg_ids))
    counts['maxquant+percolator'] = len(set(percolator_pg_ids))

    #Construct venn figure
    construct_venn(counts, args.venn_file)

    #Check PEP improvement
    




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

#Construct Venn diagram of counts per datagroup
def construct_venn(counts, venn_file):

    #2group venn
    if(len(counts.keys())==3):
        subsets = (counts['maxquant'], counts['percolator'], counts['maxquant+percolator'])
        set_labels = ('1', '2')
        textstr = "1: MaxQuant\n2: Percolator"

        fig, ax = plt.subplots(nrows=1, ncols=1)
        v = venn2(subsets=subsets, set_labels=set_labels, ax=ax)

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(-0.03, 0.07, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
        plt.tight_layout()

        fig.savefig(venn_file)
        plt.close(fig)

#Check overlap
def check_overlap(mq_proteoforms, percolator_proteoforms):

    #Init
    overlap = defaultdict(lambda: defaultdict())
    percolator_pg_ids = []
    maxquant_pg_ids = []

    overlap_found = False
    #Go over MQ protein groups
    for mq_idx, mq_row in mq_proteoforms.iterrows():
        mq_peptides = mq_row['proving_peptides'].split('|')
        #Per protein group, go over peptides
        for mq_peptide in mq_peptides:
            #Go over Percolator protein groups
            for percolator_idx, percolator_row in percolator_proteoforms.iterrows():
                #Check if the MQ peptide can be found in the peptides of the Percolator protein group
                if mq_peptide in percolator_row['proving_peptides'].split('|'):
                    #Check if one of the proteins from percolator is in the MQ protein group (link on protein level present)
                    for percolator_protein in percolator_row['protein_group'].split(';'):
                        if percolator_protein in mq_row['protein_group'].split(';'):
                            protein_group_id = percolator_row['protein_group']
                            overlap[protein_group_id]['percolator_protein_group'] = percolator_row['protein_group']
                            percolator_pg_ids.append(percolator_row['protein_group'])
                            overlap[protein_group_id]['mq_protein_group'] = mq_row['protein_group']
                            maxquant_pg_ids.append(mq_row['protein_group'])
                            overlap[protein_group_id]['percolator_peptides'] = percolator_row['proving_peptides']
                            overlap[protein_group_id]['mq_proving_peptides'] = mq_row['proving_peptides']
                            overlap[protein_group_id]['percolator_PEP'] = percolator_row['peptide_PEP']
                            overlap[protein_group_id]['mq_PEP'] = mq_row['peptide_PEP']
                            overlap[protein_group_id]['mq_max_proteins'] = mq_row['max_proteins']
                            overlap[protein_group_id]['mq_classification'] = mq_row['classification']
                            overlap[protein_group_id]['mq_base_proteoform'] = mq_row['base_proteoform']
                            overlap[protein_group_id]['mq_gene_ids'] = mq_row['gene_ids']
                            overlap_found = True

                            break #Quit searching proteins if overlapping protein is found
                    if overlap_found==True:
                        break #Quit searching Percolator protein groups if overlap found
            #if overlap_found==True:
                #break #Stop going over MQ peptides if overlap already found  -> Let this step be in comments: as Percolator is more peptide based, a protein group of maxquant can be divided over multiple protein groups (because of differing peptides) in Percolator. Therefore, a maxquant protein group should be able to link to multiple percolator protein groups

    return overlap, maxquant_pg_ids, percolator_pg_ids

#Load csv proteoform results
def load_csv_file(input_file):

    proteoform_results = pd.read_csv(input_file, sep=",", header=0)

    return proteoform_results



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