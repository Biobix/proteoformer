#!/usr/bin/env python
import traceback
import sys
import os
import time
import argparse
import re
from collections import defaultdict
from dict_functions import dict_funcs

__author__ = 'Steven Verbruggen'
#Execute 'python parse_percolator_PI.py -h' for more information

def main():

    starttime = time.time()

    print
    print "######################################"
    print "# Parse percolator protein inference #"
    print "######################################"
    print
    print "This program is part of the PROTEOFORMER pipeline"
    print

    #Parse arguments from command line
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS,
                                     description="This tool parses the Percolator protein inference output file.",
                                     add_help=False)

    man_args = parser.add_argument_group("Mandatory parameters")
    man_args.add_argument("--target_proteins", "-p", action="store", required=True, nargs="?", metavar="tab",
                          default="", type=str, help="Percolator target proteins results file (mandatory) ")
    man_args.add_argument("--fasta", "-f", action="store", required=True, nargs="?", metavar="tab",
                          default="", type=str, help="Fasta file used in Percolator (mandatory)")

    opt_args = parser.add_argument_group("Optional parameters")
    opt_args.add_argument("--help", "-h", action="help", help="Show help message and exit")
    opt_args.add_argument("--workdir", "-w", action="store", required=False, nargs="?", metavar="FOLDER", default="",
                        type=str, help="Working directory (default: CWD)")
    opt_args.add_argument("--q_val_threshold", "-q", action="store", required=False, nargs="?", metavar="FLOAT", default=0.01,
                          type=float, help="Q value filtering threshold for protein groups (default: 0.01)")
    opt_args.add_argument("--filtered_results", "-r", action="store", required=False, nargs="?", metavar="PATH",
                          default="filtered_target_proteins.tsv", type=str, help="TSV file of filtered results (default: filtered_target_proteins.tsv)")

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

    #Read protein results file
    print "Read proteins results file"
    sys.stdout.flush()
    protein_results, headers = read_protein_results(args.target_proteins, args.fasta)
    print str(len(protein_results.keys()))+" protein groups in Percolator output file."
    print

    #Filter protein groups from which all peptides are found in other protein groups
    print "Filter peptide redundancy and remove contaminant protein groups"
    sys.stdout.flush()
    protein_results = filter_redundancy(protein_results)
    print str(len(protein_results.keys()))+" protein groups found after peptide redundancy removal."
    print

    #Filter for q value
    print "Filter based on q value"
    sys.stdout.flush()
    protein_results = filter_q_val(protein_results, args.q_val_threshold)
    print str(len(protein_results.keys()))+" protein groups after q value filtering."
    print


    #####FOR DEV purposes###
    #protein_results = load_filtered_proteins(args.filtered_results)
    ####

    print "Print filtered results"
    print
    sys.stdout.flush()
    print_results(protein_results, args.filtered_results, headers)

    #Count per source
    print "Count protein groups per source of proteins"
    sys.stdout.flush()
    counts_per_source = count_protein_groups(protein_results)
    dict_funcs.print_dict(counts_per_source)
    print

    # End of program message
    print
    print "-----------------------"
    print "[%s]: PROGRAM COMPLETE" % (convert_time(time.time() - starttime))
    print "-----------------------"
    print
    sys.stdout.flush()

##########
#  SUBS  #
##########

#FOR DEV PURPOSES
def load_filtered_proteins(filtered_results):

    #Init
    protein_results = defaultdict(lambda: defaultdict(lambda: defaultdict()))
    headers = []

    with open(filtered_results, 'r') as FR:
        lines = FR.readlines()
        lines = map(lambda x: x.rstrip('\n'), lines)

        # Parse column headers
        headers = re.split('\t', lines[0])

        for line_nr in range(1, len(lines)):
            #Parse data per line
            elements = re.split('\t', lines[line_nr])
            groupId = int(elements[1])
            proteinName = elements[0]
            m_con = re.search('^CON', proteinName)
            if m_con:
                protein_results[groupId] = "Contaminants"
            if protein_results[groupId]=="Contaminants":
                continue
            for i in range(0, len(elements)-1):
                protein_results[groupId][proteinName][headers[i]] = elements[i]
            #parse peptide IDs
            protein_results[groupId][proteinName]['peptideIds'] = re.split('\s', elements[len(elements)-1].rstrip())

    return protein_results

def print_results(protein_results, filtered_results, headers):

    with open(filtered_results, 'w') as FW:
        FW.write(("\t".join(headers))+"\n")
        for protein_group in sorted(protein_results.keys()):
            for protein in protein_results[protein_group]:
                line = ""
                for header in headers[:-1]:
                    if line=="":
                        line = protein_results[protein_group][protein][header]
                    else:
                        if 'source' not in protein_results[protein_group][protein].keys():
                            print "PRINT ERROR: source not found: "+str(protein_group)+", "+str(protein)
                        line = line + "\t" + protein_results[protein_group][protein][header]
                line = line + "\t" + (" ".join(protein_results[protein_group][protein]['peptideIds']))
                line += "\n"
                FW.write(line)

    return

def count_protein_groups(protein_results):

    #Init
    counts = defaultdict()

    for protein_group in protein_results:
        if protein_results[protein_group]=="Contaminants":
            continue
        sources_list = []
        for protein in protein_results[protein_group]:
            sources_prot = re.split('\+', protein_results[protein_group][protein]['source'])
            sources_list.extend(sources_prot)

        #Make sources list unique
        sources_list = list(set(sources_list))
        sources_list = sorted(sources_list, cmp=cmp_groups)
        sources = '+'.join(sources_list)

        #Add +1 to the combination of sources found
        if sources not in counts.keys():
            counts[sources] = 1
        else:
            counts[sources]+=1

    return counts

def cmp_groups(x,y):
    if x < y:
        return -1
    else:
        return 1

def filter_q_val(protein_results, threshold):

    #Init
    filtered_protein_results = defaultdict(lambda: defaultdict(lambda: defaultdict()))

    for protein_group in protein_results.keys():
        if protein_results[protein_group] == "Contaminants":
            continue
        for protein in protein_results[protein_group].keys():
            if float(protein_results[protein_group][protein]['q-value'])<=threshold:
                filtered_protein_results[protein_group][protein] = protein_results[protein_group][protein]

    return filtered_protein_results

def filter_redundancy(red_protein_results):

    #Init
    non_red_protein_results = defaultdict(lambda: defaultdict(lambda: defaultdict()))

    #Make a copy
    red_protein_results2 = red_protein_results.copy()

    #For all protein groups
    for red_protein_group in red_protein_results.keys():
        if red_protein_results[red_protein_group] == "Contaminants":
            continue #Skip contaminants
        #For each protein of the protein group
        for red_protein in red_protein_results[red_protein_group].keys():
            redundant_protein_found = "F"
            save_under_group_number = -1
            #Compare against all protein groups
            for red_protein_group2 in red_protein_results2.keys():
                #Skip contaminants
                if red_protein_results2[red_protein_group2]=="Contaminants":
                    continue
                #Do not compare identical protein group IDs
                if red_protein_group==red_protein_group2:
                    continue
                #And all proteins
                for red_protein2 in red_protein_results2[red_protein_group2].keys():
                    #Check if all peptides (make unique list first) are present in the peptide IDs of protein 2
                    all_peptides_found = "T"
                    for peptide in list(set(red_protein_results[red_protein_group][red_protein]['peptideIds'])):
                        #Check if peptide can be found in peptide of protein 2
                        if peptide not in list(set(red_protein_results2[red_protein_group2][red_protein2]['peptideIds'])):
                            all_peptides_found = "F"
                            break #Stop search against that protein 2 if a peptide is not found
                    if all_peptides_found=="T":
                        if len(list(set(red_protein_results[red_protein_group][red_protein]['peptideIds'])))<\
                                len(list(set(red_protein_results2[red_protein_group2][red_protein2]['peptideIds']))):
                            #All peptides were found and both protein has lower amount of peptides -> redundant protein group
                            redundant_protein_found = "T"
                            break #Stop search for redundant protein groups as a redundant group was already found
                        elif len(list(set(red_protein_results[red_protein_group][red_protein]['peptideIds'])))==\
                                len(list(set(red_protein_results2[red_protein_group2][red_protein2]['peptideIds']))):
                            if red_protein_group2<red_protein_group:
                                #Save together under the lowest protein group number
                                if save_under_group_number!=-1 and red_protein_group2<save_under_group_number:
                                    save_under_group_number = red_protein_group2
                                elif save_under_group_number==-1:
                                    save_under_group_number = red_protein_group2
                        else:
                            print "Protein group "+str(red_protein_group)+" was found overlapping with protein group "+str(red_protein_group2)
                            sys.stdout.flush()
                if redundant_protein_found=="T":
                    break #Stop search (also over protein groups) as a redundant group was already found
            if redundant_protein_found=="F":
                if save_under_group_number==-1:
                    #If no redundant protein was found -> add to non-redundant proteins
                    non_red_protein_results[red_protein_group][red_protein] = red_protein_results[red_protein_group][red_protein]
                else:
                    #There is a lower protein group with exactly the same peptides. Transfer to that protein group
                    non_red_protein_results[save_under_group_number][red_protein] = red_protein_results[red_protein_group][red_protein]

    return non_red_protein_results

def read_protein_results(file, fasta):

    #Init
    protein_results = defaultdict(lambda: defaultdict(lambda: defaultdict()))
    headers = []

    #Read fasta file
    accessions = parse_accessions(fasta)

    with open(file, 'r') as FR:
        lines = FR.readlines()
        lines = map(lambda x: x.rstrip('\n'), lines)

        # Parse column headers
        headers = re.split('\t', lines[0])

        for line_nr in range(1, len(lines)):
            #Parse data per line
            elements = re.split('\t', lines[line_nr])
            groupId = int(elements[1])
            proteinName = elements[0]
            m_con = re.search('^CON', proteinName)
            m_rev = re.search('^REV', proteinName)
            if m_con:
                protein_results[groupId] = "Contaminants"
            if protein_results[groupId]=="Contaminants":
                continue
            if m_rev:
                continue
            for i in range(0, len(elements)-1):
                protein_results[groupId][proteinName][headers[i]] = elements[i]
            #parse peptide IDs
            protein_results[groupId][proteinName]['peptideIds'] = re.split('\s', elements[len(elements)-1].rstrip())
            protein_results[groupId][proteinName]['source'] = accessions[proteinName]['source']

    #Adapt headers
    last = headers[-1]
    headers[-1] = 'source'
    headers.append(last)

    return protein_results, headers

#Parse accessions of combined fasta file
def parse_accessions(combined_fasta):

    #Init
    #accessions->{main_acc}->{feature} = value
    accessions = defaultdict(lambda: defaultdict())

    with open(combined_fasta, 'r') as FR:
        for line in FR:

            #Filter out decoys and contaminants
            m_con = re.search('^>CON', line)
            m_rev = re.search('^>REV', line)
            if m_con or m_rev:
                continue

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
                    m2_snp = re.search('^(ENST\d+?)\_(\S+?)\_(\d+?)\_(\S+?)\_(\d+?)\_(\d+?)db(\d+?)\|(ENSG\d+?) (\S+?) (\S+?) ', part)
                    m2 = re.search('^(ENST\d+?)\_(\S+?)\_(\d+?)\_(\S+?)\_(\d+?)db(\d+?)\|(ENSG\d+?) (\S+?) (\S+?) ', part)
                    if m2_snp:
                        accessions[main_acc]['tr_stable_id'] = m2_snp.group(1)
                        accessions[main_acc]['chr'] = m2_snp.group(2)
                        accessions[main_acc]['start'] = m2_snp.group(3)
                        accessions[main_acc]['annotation'] = m2_snp.group(4)
                        accessions[main_acc]['snp_id'] = m2_snp.group(5)
                        accessions[main_acc]['bincode'] = m2_snp.group(6)
                        accessions[main_acc]['main_db'] = m2_snp.group(7)
                        accessions[main_acc]['gene_stable_id'] = m2_snp.group(8)
                        accessions[main_acc]['start_codon'] = m2_snp.group(9)
                        accessions[main_acc]['aTIS_call'] = m2_snp.group(10)
                        accessions[main_acc]['protein_name'] = ""
                        accessions[main_acc]['description'] = ""
                        #Search possible side accessions
                        m2_snp = re.search('\[(.+?)\]$', descr)
                        if m2_snp:
                            accessions[main_acc]['side_accessions'] = m2_snp.group(1)
                        else:
                            accessions[main_acc]['side_accessions'] = ""
                    elif m2:
                        accessions[main_acc]['tr_stable_id'] = m2.group(1)
                        accessions[main_acc]['chr'] = m2.group(2)
                        accessions[main_acc]['start'] = m2.group(3)
                        accessions[main_acc]['annotation'] = m2.group(4)
                        accessions[main_acc]['snp_id'] = ''
                        accessions[main_acc]['bincode'] = m2.group(5)
                        accessions[main_acc]['main_db'] = m2.group(6)
                        accessions[main_acc]['gene_stable_id'] = m2.group(7)
                        accessions[main_acc]['start_codon'] = m2.group(8)
                        accessions[main_acc]['aTIS_call'] = m2.group(9)
                        accessions[main_acc]['protein_name'] = ""
                        accessions[main_acc]['description'] = ""
                        #Search possible side accessions
                        m2 = re.search('\[(.+?)\]$', descr)
                        if m2:
                            accessions[main_acc]['side_accessions'] = m2.group(1)
                        else:
                            accessions[main_acc]['side_accessions'] = ""

                #For Uniprot main accessions
                else:
                    m3 = re.search('^.+?\|(.+)OS=', part)
                    if m3:
                        accessions[main_acc]['protein_name'] = main_acc
                        accessions[main_acc]['description'] = m3.group(1).replace(',', '')
                        accessions[main_acc]['tr_stable_id'] = ''
                        accessions[main_acc]['chr'] = ''
                        accessions[main_acc]['start'] = ''
                        accessions[main_acc]['annotation'] = ''
                        accessions[main_acc]['snp_id']=''
                        accessions[main_acc]['bincode'] = ''
                        accessions[main_acc]['main_db'] = ''
                        accessions[main_acc]['gene_stable_id'] = ''
                        accessions[main_acc]['start_codon'] = ''
                        accessions[main_acc]['aTIS_call'] = ''
                        #Parse side accessions
                        m4 = re.search('OS=.+?\[(.+?)\]', descr)
                        if m4:
                            accessions[main_acc]['side_accessions'] = m4.group(1)
                            #Check if there is PROTEOFORMER in the side accessions
                            m5 = re.search('ENST', accessions[main_acc]['side_accessions'])
                            if m5:
                                accessions[main_acc]["source"] = accessions[main_acc]["source"] + "+Proteoformer"
                            #Search for bincode of the first accession
                            m6_snp = re.search('^(ENST\d+?)\_(\S+?)\_(\d+?)\_(\S+?)\_(\d+?)\_(\d+?)db(\d+?)', accessions[main_acc]['side_accessions'])
                            m6 = re.search('^(ENST\d+?)\_(\S+?)\_(\d+?)\_(\S+?)\_(\d+?)db(\d+?)', accessions[main_acc]['side_accessions'])
                            if m6_snp:
                                accessions[main_acc]['tr_stable_id'] = m6_snp.group(1)
                                accessions[main_acc]['chr'] = m6_snp.group(2)
                                accessions[main_acc]['start'] = m6_snp.group(3)
                                accessions[main_acc]['annotation'] = m6_snp.group(4)
                                accessions[main_acc]['snp_id'] = m6_snp.group(5)
                                accessions[main_acc]['bincode'] = m6_snp.group(6)
                                accessions[main_acc]['main_db'] = m6_snp.group(7)
                            elif m6:
                                accessions[main_acc]['tr_stable_id'] = m6.group(1)
                                accessions[main_acc]['chr'] = m6.group(2)
                                accessions[main_acc]['start'] = m6.group(3)
                                accessions[main_acc]['annotation'] = m6.group(4)
                                accessions[main_acc]['snp_id'] = ''
                                accessions[main_acc]['bincode'] = m6.group(5)
                                accessions[main_acc]['main_db'] = m6.group(6)
                        else:
                            accessions[main_acc]['side_accessions'] = ""

    return accessions

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
