import traceback
import time
import sys
import argparse
import re
import numpy as np
import pandas as pd
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
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
    man_args.add_argument("--mq_peptides", "-M", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="MaxQuant peptides output file (mandatory)")
    man_args.add_argument("--mq_mod_peptides", "-x", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="MaxQuant modification specific peptides output file (mandatory)")
    man_args.add_argument("--percolator_peptides", "-P", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="Percolator peptides output file (mandatory)")

    opt_args = parser.add_argument_group("Optional parameters")
    opt_args.add_argument("--help", "-h", action="help", help="Show help message and exit")
    opt_args.add_argument("--venn_file", "-v", action="store", required=False, nargs="?", metavar="PATH",
                          type=str, default="overlap_mq_percolator.png", help="Venn diagram output file of overlap (default: overlap_mq_percolator.png)")
    opt_args.add_argument("--percolator_uniq_out_csv", "-u", action="store", required=False, nargs="?", metavar="PATH",
                          type=str, default="percolator_uniq_proteoforms.csv", help="Percolator unique proteoforms output csv file (default: percolator_uniq_proteoforms.csv)")

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
    proteoform_overlap, mq_overlap_pg_ids, percolator_overlap_pg_ids = check_overlap(mq_proteoforms, percolator_proteoforms, args.mq_peptides, args.percolator_peptides, args.mq_mod_peptides)
    #print_dict(proteoform_overlap)

    print("\n\nOverlap: " + str(len(proteoform_overlap.keys())))
    print("Maxquant protein group IDs in overlap: "+str(len(set(mq_overlap_pg_ids))))
    print("Percolator protein group IDs in overlap: "+str(len(set(percolator_overlap_pg_ids)))+"\n")

    #Get MQ-unique protein groups
    mq_uniq_pgs = []
    for mq_idx, mq_row in mq_proteoforms.iterrows():
        if mq_row['protein_group'] not in mq_overlap_pg_ids:
            mq_uniq_pgs.append(mq_row['protein_group'])
    #Get Percolator-uniqe protein groups
    percolator_uniq_pgs = []
    for percolator_idx, percolator_row in percolator_proteoforms.iterrows():
        if percolator_row['protein_group'] not in percolator_overlap_pg_ids:
            percolator_uniq_pgs.append(percolator_row['protein_group'])
    #Construct counts
    counts = defaultdict()
    counts['maxquant'] = len(mq_uniq_pgs)
    counts['percolator'] = len(percolator_uniq_pgs)
    counts['maxquant+percolator'] = len(set(percolator_overlap_pg_ids))

    #Construct venn figure
    construct_venn(counts, args.venn_file)

    #Check PEP improvement
    PEP_improvement(proteoform_overlap, mq_overlap_pg_ids, percolator_overlap_pg_ids, mq_proteoforms,
                    percolator_proteoforms)

    #Export Percolator unique proteoforms
    export_percolator_uniq_proteoforms(percolator_uniq_pgs, percolator_proteoforms, args.percolator_uniq_out_csv)

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

#Export Percolator unique proteoforms
def export_percolator_uniq_proteoforms(uniq_pgs, proteoforms, out_file):

    #Select Percolator unique proteoforms
    mask = []
    for idx, row in proteoforms.iterrows():
        if row['protein_group'] in uniq_pgs:
            mask.append(True)
        else:
            mask.append(False)
    uniq_proteoforms = proteoforms[mask]

    #Export to csv
    uniq_proteoforms.to_csv(out_file, sep=",", header=True, index=False)

    return

#Check PEP improvement
def PEP_improvement(proteoform_overlap, mq_overlap_proteingroups, percolator_overlap_proteingroups, mq_proteoforms,
                    percolator_proteoforms):

    #Init
    mq_overlap_peps = []
    percolator_overlap_peps = []
    mq_non_percolator_overlap_peps = []

    #Check PEP improvement in the overlap
    #For each protein group, check for common peptides
    for pg in proteoform_overlap.keys():
        for mq_peptide in proteoform_overlap[pg]['mq_proving_peptides'].split('|'):
            if mq_peptide in proteoform_overlap[pg]['percolator_peptides'].split('|'):
                #Get index of peptides for both MQ and Percolator
                mq_index = proteoform_overlap[pg]['mq_proving_peptides'].split('|').index(mq_peptide)
                percolator_index = proteoform_overlap[pg]['percolator_peptides'].split('|').index(mq_peptide)
                #Get the peptide PEPs and store in list
                mq_overlap_peps.append(float(proteoform_overlap[pg]['mq_PEP'].split('|')[mq_index]))
                percolator_overlap_peps.append(float(proteoform_overlap[pg]['percolator_PEP'].split('|')[percolator_index]))
            else:
                #Also keep the MQ peptides that are not present in Percolator for the overlapping peptides
                mq_index = proteoform_overlap[pg]['mq_proving_peptides'].split('|').index(mq_peptide)
                mq_non_percolator_overlap_peps.append(float(proteoform_overlap[pg]['mq_PEP'].split('|')[mq_index]))

    print("Median of the MaxQuant peptide PEPs of the overlapping proteoforms: "+str(np.median(mq_overlap_peps)))
    print("Median of the Percolator peptide PEPs of the overlapping proteoforms: "+str(np.median(percolator_overlap_peps))+"\n")
    print("Median of the MaxQuant but non-Percolator peptide PEPs of the overlapping proteoforms: " + str(np.median(mq_non_percolator_overlap_peps)) + "\n")

    #Take log10 of list
    mq_overlap_peps_log10 = np.log10(mq_overlap_peps)
    percolator_overlap_peps_log10 = np.log10(percolator_overlap_peps)
    mq_non_percolator_overlap_peps_log10 = np.log10(mq_non_percolator_overlap_peps)

    #Overlap Plots
    sns.set_style(style="darkgrid")
    binwidth=2
    bins = np.arange(-16-binwidth,0+binwidth,binwidth)
    #bins = np.arange(round(np.min(mq_overlap_peps_log10),-1)-binwidth,round(np.max(mq_overlap_peps_log10),-1)+binwidth,binwidth)
    fig, (ax1, ax3, ax5, ax6) = plt.subplots(nrows=4, ncols=1, figsize=(14,26))
    sns.distplot(mq_overlap_peps_log10, hist=True, bins=bins, kde=True, ax=ax1, kde_kws={'color':'#f58d42', 'lw':2}, hist_kws={'color':'#E0BC99'})
    ax1.set_xlim([-16, 5])
    ax1.set_title('MaxQuant peptide PEPs of the overlapping proteoforms', fontsize=28, fontweight='bold')
    ax1.set_xlabel('Peptide Posterior error probability (log10)', fontsize=24)
    ax1.set_ylabel('Abundance', fontsize=24)
    ax1.tick_params(axis='x', which='major', labelsize=20)
    ax1.tick_params(axis='y', which='major', labelsize=20)
    #ax1.text(-0.09, 1.02, 'A', transform=ax1.transAxes,size=32, weight='bold')

    '''
    bins = np.arange(round(np.min(mq_overlap_peps_log10),-1)-binwidth,round(np.max(mq_overlap_peps_log10),-1)+binwidth,binwidth)
    sns.distplot(mq_overlap_peps_log10, hist=True, bins=bins, kde=True, ax=ax2, kde_kws={'color':'#f58d42', 'lw':2}, hist_kws={'color':'#E0BC99'})
    ax2.set_xlim([-20, 5])
    ax2.set_title('MaxQuant peptide PEPs of the overlapping proteoforms (zoom)', fontsize=28)
    ax2.set_ylabel('Abundance', fontsize=24)
    ax2.tick_params(axis='x', which='major', labelsize=20)
    ax2.tick_params(axis='y', which='major', labelsize=20)
    ax2.text(-0.09, 1.02, 'B', transform=ax2.transAxes, size=32, weight='bold')
    '''

    #bins = np.arange(round(np.min(percolator_overlap_peps_log10),-1)-binwidth,round(np.max(percolator_overlap_peps_log10),-1)+binwidth,binwidth)
    sns.distplot(percolator_overlap_peps_log10, hist=True, bins=bins, kde=True, ax=ax3, kde_kws={'color':'#f58d42', 'lw':2}, hist_kws={'color':'#E0BC99'})
    ax3.set_xlim([-16, 5])
    ax3.set_title('Percolator peptide PEPs of the overlapping proteoforms', fontsize=28, fontweight='bold')
    ax3.set_xlabel('Peptide Posterior error probability (log10)', fontsize=24)
    ax3.set_ylabel('Abundance', fontsize=24)
    ax3.tick_params(axis='x', which='major', labelsize=20)
    ax3.tick_params(axis='y', which='major', labelsize=20)
    #ax3.text(-0.09, 1.02, 'B', transform=ax3.transAxes, size=32, weight='bold')

    '''
    bins = np.arange(round(np.min(mq_non_percolator_overlap_peps_log10),-1)-binwidth,round(np.max(mq_non_percolator_overlap_peps_log10),-1)+binwidth,binwidth)
    sns.distplot(mq_non_percolator_overlap_peps_log10, hist=True, bins=bins, kde=True, ax=ax4, kde_kws={'color':'#2032f7', 'lw':2}, hist_kws={'color':'#4278cf'})
    ax4.set_xlim([-16, 5])
    ax4.set_title('MaxQuant but non-Percolator peptide PEPs of the overlapping proteoforms', fontsize=24, fontweight='bold')
    ax4.set_xlabel('Peptide Posterior error probability (log10)', fontsize=24)
    ax4.set_ylabel('Abundance', fontsize=24)
    ax4.tick_params(axis='x', which='major', labelsize=20)
    ax4.tick_params(axis='y', which='major', labelsize=20)
    ax4.text(-0.09, 1.02, 'D', transform=ax4.transAxes, size=32, weight='bold')
    '''

    ## Non-overlap PEP distributions ##

    #Get MQ-unique peptide PEPs
    mq_uniq_peps=[]
    for mq_idx, mq_row in mq_proteoforms.iterrows():
        if mq_row['protein_group'] not in mq_overlap_proteingroups:
            for PEP_mq in mq_row['peptide_PEP'].split('|'):
                mq_uniq_peps.append(float(PEP_mq))

    #Get Percolator-unique PEPs
    percolator_uniq_peps=[]
    for percolator_idx, percolator_row in percolator_proteoforms.iterrows():
        if percolator_row['protein_group'] not in percolator_overlap_proteingroups:
            for PEP_percolator in percolator_row['peptide_PEP'].split('|'):
                percolator_uniq_peps.append(float(PEP_percolator))

    print("Median of the peptide PEPs of the MaxQuant unique proteoforms: "+str(np.median(mq_uniq_peps)))
    print("Median of the peptide PEPs of the Percolator unique proteoforms: "+str(np.median(percolator_uniq_peps)))

    #Take log10 of lists
    mq_uniq_peps_log10 = np.log10(mq_uniq_peps)
    percolator_uniq_peps_log10 = np.log10(percolator_uniq_peps)

    #bins = np.arange(round(np.min(mq_uniq_peps_log10),-1)-binwidth,round(np.max(mq_uniq_peps_log10),-1)+binwidth,binwidth)
    sns.distplot(mq_uniq_peps_log10, hist=True, bins=bins, kde=True, ax=ax5, kde_kws={'color':'#fa3434', 'lw':2}, hist_kws={'color':'#FF9999'})
    ax5.set_xlim([-16, 5])
    ax5.set_title('Peptide PEPs of MaxQuant unique proteoforms', fontsize=24, fontweight='bold')
    ax5.set_xlabel('Peptide Posterior error probability (log10)', fontsize=24)
    ax5.set_ylabel('Abundance', fontsize=24)
    ax5.tick_params(axis='x', which='major', labelsize=20)
    ax5.tick_params(axis='y', which='major', labelsize=20)
    #ax5.text(-0.09, 1.02, 'E', transform=ax5.transAxes, size=32, weight='bold')

    #bins = np.arange(round(np.min(percolator_uniq_peps_log10),-1)-binwidth,round(np.max(percolator_uniq_peps_log10),-1)+binwidth,binwidth)
    sns.distplot(percolator_uniq_peps_log10, hist=True, bins=bins, kde=True, ax=ax6, kde_kws={'color':'#36cf36', 'lw':2}, hist_kws={'color':'#99CC99'})
    ax6.set_xlim([-16, 5])
    ax6.set_title('Peptide PEPs of Percolator unique proteoforms', fontsize=24, fontweight='bold')
    ax6.set_xlabel('Peptide Posterior error probability (log10)', fontsize=24)
    ax6.set_ylabel('Abundance', fontsize=24)
    ax6.tick_params(axis='x', which='major', labelsize=20)
    ax6.tick_params(axis='y', which='major', labelsize=20)
    #ax6.text(-0.09, 1.02, 'F', transform=ax6.transAxes, size=32, weight='bold')

    plt.tight_layout()
    fig.savefig('overlap_PEPs.png')
    #fig.savefig('overlap_PEPs.svg', format='svg')

    return


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
        #fig.savefig(venn_file, format='svg')
        plt.close(fig)

#Check overlap
def check_overlap(mq_proteoforms, percolator_proteoforms, mq_peptides_file, percolator_peptides_file, mq_mod_peptides_file):

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
                            overlap[protein_group_id]['mq_base_proteoform'] = mq_row['base_proteoform']
                            overlap_found = True

                            break #Quit searching proteins if overlapping protein is found
                    if overlap_found==True:
                        break #Quit searching Percolator protein groups if overlap found
            #if overlap_found==True:
                #break #Stop going over MQ peptides if overlap already found  -> Let this step be in comments: as
                        # Percolator is more peptide based, a protein group of maxquant can be divided over multiple
                        # protein groups (because of differing peptides) in Percolator. Therefore, a maxquant protein
                        # group should be able to link to multiple percolator protein groups

    #Check extra for overlapping peptides
    #Due to the difference between the PI strategies of MQ (max scoring protein system, differing peptide needed at
    # both sides to have separate PGs) and Percolator (differing peptide needed at one side for different PG),
    #some PGs do not occur in the overlap
    #We should check what the peptide PEPs are at the opposite side and if one of the prooving peptides is identified
    # at the other side (PEP low enough), bring them to the overlap portion

    #Get MQ-unique protein groups after first overlap run
    mq_uniq_pgs = []
    for mq_idx, mq_row in mq_proteoforms.iterrows():
        if mq_row['protein_group'] not in maxquant_pg_ids:
            mq_uniq_pgs.append(mq_row['protein_group'])
    #Get Percolator-uniqe protein groups after first overlap run
    percolator_uniq_pgs = []
    for percolator_idx, percolator_row in percolator_proteoforms.iterrows():
        if percolator_row['protein_group'] not in percolator_pg_ids:
            percolator_uniq_pgs.append(percolator_row['protein_group'])

    #Read peptides of both runs
    peptides_mq = read_peptides_maxquant(mq_peptides_file)
    mod_peptides_mq = read_mod_peptides_maxquant(mq_mod_peptides_file)
    peptides_percolator = read_peptides_percolator(percolator_peptides_file)

    #Then for each of these protein groups, check the peptide PEPs in the other dataset and if they are all <0.01, put them in the overlap
    #print("Total amount of initial MQ unique protein groups: "+str(len(mq_uniq_pgs))+"\n")
    for mq_pg in mq_uniq_pgs:
        #print(mq_pg)
        #Check if one of the prooving MQ peptides was identified in Percolator and if so if it is in Percolator under the PEP threshold
        a_pep_under_th = False
        all_percolator_protein_groups = []
        all_percolator_peptides=[]
        all_percolator_PEPs=[]
        for mq_peptide in mq_proteoforms[mq_proteoforms['protein_group']==mq_pg]['proving_peptides'].values[0].split('|'):
            #print(mq_peptide)
            if len(mq_peptide)>30: #Too long peptides could not be identified by Prosit. Put them in the overlap as we cannot know for sure if Percolator+Prosit would miss them or not
                all_percolator_protein_groups.append(mq_pg)
                all_percolator_peptides.append("TooLongPeptide")
                all_percolator_PEPs.append("TooLongPeptide")
                a_pep_under_th = True
            else:
                if mq_peptide in list(peptides_percolator['stripped_seqs']):
                    #print(peptides_percolator[peptides_percolator['stripped_seqs']==mq_peptide]['posterior_error_prob'].values[0])
                    #print(peptides_percolator[peptides_percolator['stripped_seqs']==mq_peptide]['proteinIds'].values[0])
                    if peptides_percolator[peptides_percolator['stripped_seqs']==mq_peptide]['posterior_error_prob'].values[0]<0.01:
                        all_percolator_protein_groups.append(peptides_percolator[peptides_percolator['stripped_seqs'] == mq_peptide]['proteinIds'].values[0])
                        all_percolator_peptides.append(mq_peptide)
                        all_percolator_PEPs.append(str(peptides_percolator[peptides_percolator['stripped_seqs'] == mq_peptide]['posterior_error_prob'].values[0]))
                        a_pep_under_th=True
                else:
                    if mq_peptide in list(mod_peptides_mq['Sequence']):
                        if mod_peptides_mq[mod_peptides_mq['Sequence'] == mq_peptide]['Modifications'].values[0] == "Acetyl (Protein N-term)":
                            #Acetylated peptides could also not be identified by Prosit.
                            all_percolator_protein_groups.append(mq_pg)
                            all_percolator_peptides.append("AcetylatedPeptide")
                            all_percolator_PEPs.append("AcetylatedPeptide")
                            a_pep_under_th = True
        if a_pep_under_th==True:
            #Then save in overlap
            protein_group_id = all_percolator_protein_groups[0] #Select the first as key
            overlap[protein_group_id]['percolator_protein_group'] = all_percolator_protein_groups[0]
            overlap[protein_group_id]['all_percolator_protein_groups'] = all_percolator_protein_groups
            percolator_pg_ids.append(all_percolator_protein_groups[0])
            overlap[protein_group_id]['percolator_peptides'] = '|'.join(all_percolator_peptides)
            overlap[protein_group_id]['percolator_PEP'] = '|'.join(all_percolator_PEPs)
            overlap[protein_group_id]['mq_protein_group'] = mq_pg
            maxquant_pg_ids.append(mq_pg)
            overlap[protein_group_id]['mq_proving_peptides'] = mq_proteoforms[mq_proteoforms['protein_group']==mq_pg]['proving_peptides'].values[0]
            overlap[protein_group_id]['mq_PEP'] = mq_proteoforms[mq_proteoforms['protein_group']==mq_pg]['peptide_PEP'].values[0]
            overlap[protein_group_id]['mq_max_proteins'] = mq_proteoforms[mq_proteoforms['protein_group'] == mq_pg]['max_proteins'].values[0]
            overlap[protein_group_id]['mq_base_proteoform'] = mq_proteoforms[mq_proteoforms['protein_group'] == mq_pg]['base_proteoform'].values[0]
            #print_dict(overlap[protein_group_id])

    #print("Amount of protein groups to bring to overlap: "+str(count))

    return overlap, maxquant_pg_ids, percolator_pg_ids

def read_mod_peptides_maxquant(in_file):

    mod_peptides = pd.read_csv(in_file, sep="\t", header=0)
    mod_peptides.columns = [strip_non_ascii(x) for x in mod_peptides.columns]

    return mod_peptides

def strip_non_ascii(string):
    #Returns the string without non ASCII characters
    stripped = (c for c in string if 0 < ord(c) < 127)
    return ''.join(stripped)

#Read peptides input info from Percolator
def read_peptides_percolator(in_file):

    #Load peptides file
    peptides = pd.read_csv(in_file, sep="\t", header=0)
    peptides.loc[:, 'stripped_seqs'] = [elem.strip('_.') for elem in list(peptides.loc[:, 'peptide'])]

    return peptides

#Read peptides input info from maxquant
def read_peptides_maxquant(in_file):

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