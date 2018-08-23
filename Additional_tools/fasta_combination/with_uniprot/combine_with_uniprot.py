import traceback
import sys
import os
import time
import argparse
import re
from collections import defaultdict
import itertools
#from dict_functions import dict_funcs
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3

__author__ = 'Steven Verbruggen'
#Execute 'python combine_with_uniprot.py -h' for more information

def main():

    starttime = time.time()

    print
    print "#######################"
    print "# Combine wit UniProt #"
    print "#######################"
    print
    print "This program is part of the PROTEOFORMER pipeline"
    print

    #Parse arguments from command line
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS,
                                     description="This tool is part of the PROTEOFORMER pipeline. It combines "
                                                 "a (combined) fasta file of PROTEOFORMER with UniProt, ready for MS inspection.",
                                     add_help=False)

    man_args = parser.add_argument_group("Mandatory parameters")
    man_args.add_argument("--fasta", "-f", action="store", required=True, nargs="?", metavar="PATH",
                          default="", type=str, help="Input (combined) fasta file (mandatory) ")
    man_args.add_argument("--uniprot", "-u", action="store", required=True, nargs="?", metavar="PATH",
                          default="", type=str, help="Input UniProt fasta file (mandatory)")

    opt_args = parser.add_argument_group("Optional parameters")
    opt_args.add_argument("--help", "-h", action="help", help="Show help message and exit")
    opt_args.add_argument("--workdir", "-w", action="store", required=False, nargs="?", metavar="FOLDER", default="",
                        type=str, help="Working directory (default: CWD)")
    opt_args.add_argument("--output_fasta", "-o", action="store", required=False, nargs="?", metavar="PATH",
                          default="comb_fasta_uniprot_proteoformer.fasta", type=str, help="Output fasta file with combination of "
                                                                                            "data of Uniprot and PROTEOFORMER "
                                                                                            "(default: comb_fasta_uniprot_proteoformer.fasta)")
    opt_args.add_argument("--overview_file", "-t", action="store", required=False, nargs="?", metavar="PATH",
                            default="uniprot_proteoformer_overview.txt", type=str, help="Path to a file where information about "
                                                                            "the input files and the database combinations are stored "
                                                                            "(default: uniprot_proteoformer_overview.txt)")
    opt_args.add_argument("--venn_diagram", "-d", action="store", required=False, nargs="?", metavar="PATH",
                          default="venn_diagram.png", type=str, help="Path to a png file where the Venn diagram will "
                                                                     "be stored (default: venn_diagram.png)")

    args = parser.parse_args()

    #default workdir is CWD
    if args.workdir == "":
        args.workdir = os.getcwd()
    if args.output_fasta == "comb_fasta_uniprot_proteoformer.fasta":
        args.output_fasta = args.workdir + "/" + args.output_fasta
    if args.overview_file == "uniprot_proteoformer_overview.txt":
        args.overview_file = args.workdir + "/" + args.overview_file
    if args.venn_diagram == "venn_diagram.png":
        args.venn_diagram = args.workdir + "/" + args.venn_diagram

    #List parameters
    print "Parameters:"
    for arg in vars(args):
        print '    %-15s\t%s' % (arg, getattr(args, arg))
    print
    sys.stdout.flush()

    ########
    # MAIN #
    ########

    #Init
    input_data=defaultdict(lambda: defaultdict())

    #Make combination list of data sources
    sources = ["proteoformer", "uniprot"]
    inputs=[args.fasta, args.uniprot]

    #Construct combined data groups
    iter_sources = []
    for i in range(0, len(sources)):
        iter_sources.append(sources[i])
    iter_sources.append("")
    data_groups = list(itertools.combinations(iter_sources, 2))
    data_groups_conv = []
    for comb in data_groups:
        comb = sorted(comb, cmp=cmp_groups)
        data_groups_conv.append('+'.join(filter(None, comb)))
    data_groups = data_groups_conv

    #Read in fasta's
    for source_nr in range(0, len(sources)):
        input_data = read_fasta(inputs[source_nr], input_data, sources[source_nr])

    #Combine the different sources of data
    (comb_data, counts) = combine_data(input_data, data_groups)

    #Write output file
    write_output(comb_data, args.output_fasta)

    #Construct tab separated table of results
    write_overview(counts, args.overview_file)

    #Construct venn
    construct_venn(counts,args.venn_diagram)

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
        subsets = (counts['proteoformer'], counts['uniprot'], counts['uniprot+proteoformer'])
        set_labels = ('1', '2')
        textstr = "1: proteoformer\n2: uniprot"

        fig, ax = plt.subplots(nrows=1, ncols=1)
        v = venn2(subsets=subsets, set_labels=set_labels, ax=ax)

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(-0.03, 0.07, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
        plt.tight_layout()

        fig.savefig(venn_file)
        plt.close(fig)

    #3group venn
    '''
    if(len(counts.keys())==7):
        subsets = tuple(counts.values(0))
        set_labels = tuple(counts.keys())
        textstr = ""
        for i in counts.keys():
            textstr = i+"\n"

        fig, ax = plt.subplots(nrows=1, ncols=1)
        v = venn3(subsets=subsets, set_labels=set_labels, ax=ax)

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(-0.07, 0.15, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
        plt.tight_layout()

        fig.savefig(venn_file)
        plt.close(fig)
    '''

    return

#Write overview file
def write_overview(counts, overview_file):

    with open(overview_file, 'w') as FW:
        line = '{:30} {:20}'.format("Files", "Counts")
        FW.write(line+"\n")
        for data_group in counts.keys():
            line='{:30} {:20}'.format(data_group, str(counts[data_group]))
            FW.write(line+"\n")

    return

#Write output fasta file
def write_output(comb_data, fasta):

    with open(fasta, 'w') as FW:
        for acc in comb_data:
            FW.write(acc+"\n")
            FW.write(comb_data[acc]+"\n")

    return

#Combine data on sequence level
def combine_data(input_data, data_groups):

    #Init
    #comb_data_per_file->{sequence}->{proteoformer/uniprot} = accession
    comb_data_per_file = defaultdict(lambda: defaultdict())
    #comb_data->{general accession} = sequence
    comb_data = defaultdict()
    #Init counts dict
    counts={}
    for data_group in data_groups:
        counts[data_group] = 0

    #Combine on sequence level
    for source in input_data:
        for accession in input_data[source]:
            comb_data_per_file[input_data[source][accession]][source] = accession

    #Get the data group and count
    for sequence in comb_data_per_file:
        classified_in = ""
        found_in = []
        for source in sorted(comb_data_per_file[sequence].keys(), cmp=cmp_groups):
            classified_in = classified_in + source + '+'
            found_in.append(source)
        classified_in = classified_in.strip('+')
        counts[classified_in]+=1
        #Construct general accession
        main_acc=""
        side_acc=""
        i=1
        for found_data_group in sorted(found_in, cmp=cmp_groups):
            #Search main accession in data group 1
            if i==1:
                m_main = re.search('^>(.+?)\|(.+?)\|.+$',comb_data_per_file[sequence][found_data_group])
                if m_main:
                    main_acc = ">"+m_main.group(1)+"|"+m_main.group(2)
                    if m_main.group(1)!='generic':
                        m_desc = re.search('^>.+?\|.+?\|(.+)$', comb_data_per_file[sequence][found_data_group])
                        if m_desc:
                            main_acc = main_acc+"|"+m_desc.group(1)
                    else:
                        m_desc = re.search('^>.+?\|.+?\|([^\[]+)', comb_data_per_file[sequence][found_data_group])
                        if m_desc:
                            main_acc = main_acc + "|" + m_desc.group(1)
                #Save extra accessions with side accessions
                m_side = re.search('\[(.+)\]', comb_data_per_file[sequence][found_data_group])
                if m_side:
                    new_side_accs = m_side.group(1).split('#')
                    for new in new_side_accs:
                        side_acc = side_acc + new + "#"
            else:
                #For other data groups: main accessions becomes part of the side accessions
                m_main = re.search('^>.+?\|(.+?)\|.+$', comb_data_per_file[sequence][found_data_group])
                if m_main:
                    side_acc = side_acc+m_main.group(1)+"#"
                # Save extra accessions with side accessions
                m_side = re.search('\[(.+)\]', comb_data_per_file[sequence][found_data_group])
                if m_side:
                    new_side_accs = m_side.group(1).split('#')
                    for new in new_side_accs:
                        side_acc = side_acc + new + "#"
            i+=1
        #Strip last hashtag
        side_acc = side_acc.strip('#')
        #Construct total accession
        general_accession = main_acc
        if side_acc!="":
            general_accession = general_accession + " [" + side_acc + "]"
        #Save
        comb_data[general_accession] = sequence

    return (comb_data, counts)

#read fasta
def read_fasta(fasta, input_data, source):


    #Open fasta file
    with open(fasta, 'r') as FR:
        lines = FR.readlines()
        lines = map(lambda x: x.rstrip("\n"), lines)
        lines = map(lambda x: x.rstrip("\r"), lines)
        #Init
        accession = ""
        sequence = ""
        #for i in range(0,50):
        for i in range(0, len(lines)):
            # Check if accession -> new entry
            if ((re.search('^>', lines[i]))):
                if(i!=0):
                    # Save previous entry
                    input_data[source][accession] = sequence
                # Get new entry accession and empty sequence
                accession = lines[i]
                sequence = ""
            else:
                # Concat next line of sequence
                sequence += lines[i]
        # Save the last entry
        if (sequence != ""):
            input_data[source][accession] = sequence
            accession = ""
            sequence = ""

    return input_data

#Cmp for sorting data groups
def cmp_groups(x,y):
    if x=='uniprot':
        return -1
    elif y=='uniprot':
        return 1
    else:
        if x<y:
            return -1
        else:
            return 1

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