#!/usr/bin/env python
import traceback
import sys
import os
import time
import argparse
import re
from collections import defaultdict
import itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
#from print_collection import print_dict

__author__ = 'Steven Verbruggen'
#Execute 'python combine_dbs.py -h' for more information

def main():

    starttime = time.time()

    print
    print "######################"
    print "# Combine FASTA DB's #"
    print "######################"
    print
    print "This program is part of the PROTEOFORMER pipeline"
    print

    #Parse arguments from command line
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS,
                                     description="This tool is part of the PROTEOFORMER pipeline. It combines "
                                                 "different FASTA DB's into one big file, ready for MS inspection.",
                                     add_help=False)

    man_args = parser.add_argument_group("Mandatory parameters")
    man_args.add_argument("--in_files", "-i", action="store", required=True, nargs="?", metavar="COMMA-SEPARATED LIST",
                          default="", type=str, help="Comma-separated list of all paths to all input fasta files (mandatory) ")

    opt_args = parser.add_argument_group("Optional parameters")
    opt_args.add_argument("--help", "-h", action="help", help="Show help message and exit")
    opt_args.add_argument("--workdir", "-w", action="store", required=False, nargs="?", metavar="FOLDER", default="",
                        type=str, help="Working directory (default: CWD)")
    opt_args.add_argument("--output_file", "-o", action="store", required=False, nargs="?", metavar="PATH",
                          default="output_combined_dbs.fa", type=str, help="Path of the combined output fasta file"
                                                                           "(default: outpub_combined_dbs.fa)")
    opt_args.add_argument("--overview_file", "-t", action="store", required=False, nargs="?", metavar="PATH",
                            default="fasta_file_overview.txt", type=str, help="Path to a file where information about "
                                                                            "the input files and the database combinations "
                                                                              "with their binary codes will be stored "
                                                                              "(default: fasta_file_overview.txt)")
    opt_args.add_argument("--venn_diagram", "-d", action="store", required=False, nargs="?", metavar="PATH",
                          default="venn_diagram.png", type=str, help="Path to a png file where the Venn diagram will "
                                                                     "be stored (default: venn_diagram.png)")
    opt_args.add_argument("--verbose_output", "-v", action="store", required=False, nargs="?", metavar="Y/N",
                          default="Y", type=str, help="Whether the output fasta should give the full list of accessions per "
                                                      "combination. If not, redundant accessions between input files, will be "
                                                      "omitted (default: Y)")

    args = parser.parse_args()

    #default workdir is CWD
    if args.workdir == "":
        args.workdir = os.getcwd()
    #Check mandatory input
    if args.in_files == "":
        print "Do not forget to specify all input fasta files!"
        sys.exit()
    if(args.verbose_output!='Y' and args.verbose_output!='N'):
        print "Verbose output argument should be Y or N!"
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

    #Check input files
    input_files = args.in_files.split(',')
    for file in input_files:
        m = re.search('\.[fa$|fasta$]', file)
        if not m:
            print file+" is not recognised as a fasta file"
            sys.exit()

    #Make input file dict
    input_file_dict = make_input_file_dict(input_files, args.overview_file)

    #Read input files
    input_data = read_input_files(input_file_dict)

    #Construct presence possibilities matrix
    comb_dict = construct_combinations(input_file_dict, args.overview_file)

    #Put together based on sequence
    (comb_data, counts_per_bincode) = combine_data_on_seq_level(input_data, comb_dict, input_file_dict, args.verbose_output)

    #Print to output fasta file
    print_to_output(comb_data, args.output_file)

    #Add counts per bin code to the overview file
    add_counts_to_overview_file(counts_per_bincode, comb_dict, args.overview_file)

    #Make Venn diagram
    construct_venn(counts_per_bincode, args.venn_diagram, input_file_dict)

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

#Construct Venn diagram of counts per bincode
def construct_venn(counts, venn_file, input_file_dict):

    #2group venn
    if(len(counts.keys())==3):
        subsets = (counts['10'], counts['01'], counts['11'])
        set_labels = tuple(input_file_dict.keys())
        textstr = ""
        for i in input_file_dict.keys():
            textstr = textstr+str(i)+": "+input_file_dict[i]+"\n"

        fig, ax = plt.subplots(nrows=1, ncols=1)
        v = venn2(subsets=subsets, set_labels=set_labels, ax=ax)

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(-0.07, 0.15, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
        plt.tight_layout()

        fig.savefig(venn_file)
        plt.close(fig)

    #3group venn
    if(len(counts.keys())==7):
        subsets = (counts['100'], counts['010'], counts['110'], counts['001'], counts['101'], counts['011'], counts['111'])
        set_labels = tuple(input_file_dict.keys())
        textstr = ""
        for i in input_file_dict.keys():
            textstr = textstr+str(i)+": "+input_file_dict[i]+"\n"

        fig, ax = plt.subplots(nrows=1, ncols=1)
        v = venn3(subsets=subsets, set_labels=set_labels, ax=ax)

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(-0.07, 0.15, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
        plt.tight_layout()

        fig.savefig(venn_file)
        plt.close(fig)

    return

#Add count to overview file
def add_counts_to_overview_file(counts, comb_dict, overview_file):

    print
    print "Overview of ORF counts per bin code:"
    print "\t\tBinCode\t\tORFs\t\tFiles"
    for i in sorted(counts.keys(), cmp=cmp_files):
        str_files = ", ".join(comb_dict[i]['files'])
        print str(comb_dict[i]['sort_number'])+":\t\t"+i+"\t\t"+str(counts[i])+"\t\t"+str_files
    print

    with open(overview_file, 'a') as FW:
        FW.write("\n")
        FW.write("Overview of ORF counts per bin code:\n")
        FW.write("\t\tBinCode\t\tORFs\t\tFiles\n")
        for i in sorted(counts.keys(), cmp=cmp_files):
            str_files = ", ".join(comb_dict[i]['files'])
            FW.write(str(comb_dict[i]['sort_number']) + ":\t\t" + i + "\t\t" + str(counts[i]) + "\t\t" + str_files+"\n")
        FW.write("\n")

    return

#Print to output fasta file
def print_to_output(comb_data, output_file):

    with open(output_file, 'w') as FW:
        for accession in comb_data:
            FW.write(accession+"\n")
            FW.write(comb_data[accession]+"\n")

    return

#Combine the data of all fasta files. We want to combine on sequence level so that this results in unique sequences.
def combine_data_on_seq_level(input_data, comb_dict, input_file_dict, verbose_output):

    #Init
    #Structure comb_data:
    #   comb_data->{general accession} = sequence
    comb_data = defaultdict()
    #Structure comb_data_file:
    #   comb_data_per_file->{sequence}->{file_number in which present} = accession for that file and sequence
    comb_data_per_file = defaultdict(lambda: defaultdict())
    #Keep counts per bincode
    counts_per_bincode = defaultdict()
    for bincode in comb_dict:
        counts_per_bincode[bincode] = 0

    #Group per sequence
    for file_number in sorted(input_data.keys()):
        for new_acc in input_data[file_number]:
            #Accession will be added to the sequence key if it already exists, else a new will be made
            comb_data_per_file[input_data[file_number][new_acc]][file_number] = new_acc

    for sequence in comb_data_per_file:
        #Construct bin code based on sequence presence over the different files
        bincode=""
        for file_number in sorted(input_file_dict.keys()):
            if file_number in comb_data_per_file[sequence].keys():
                bincode+="1"
            else:
                bincode+="0"
        counts_per_bincode[bincode]+=1

        # Init
        main_accession = ""
        description = ""
        main_suffix = ""
        side_accession = []
        side_suffix = []
        first_file=True
        #Parse accessions to obtain main and side accessions
        #Take randomly the first accession as the main accession of the combined entry
        for file_number in sorted(comb_data_per_file[sequence].keys()):
            #For the first file where the sequence was found
            if first_file==True:
                #Keep main accession as main accession
                m_main = re.search('^>generic\|(.+?)\|([^\[]+)', comb_data_per_file[sequence][file_number])
                if m_main:
                    main_accession = m_main.group(1)
                    description = m_main.group(2).strip(" ")
                    main_suffix = bincode+"db"+str(file_number)
                first_file = False #Keep the fact that the file has been seen already
            else:
                #Not the first file: main accession becomes a side accession as well
                m_main = re.search('^>generic\|(.+?)\|', comb_data_per_file[sequence][file_number])
                if m_main:
                    new_side_accession = m_main.group(1)
                    if ((new_side_accession!=main_accession) and not (new_side_accession in side_accession) or verbose_output=='Y'):
                        side_accession.append(new_side_accession)
                        side_suffix.append(file_number)
            # Other accessions as side accessions
            m_side = re.search('\[(.+)\]', comb_data_per_file[sequence][file_number])
            if m_side:
                new_side_accessions = m_side.group(1).split('#')
                for new_side_accession in new_side_accessions:
                    if ((new_side_accession!=main_accession) and not (new_side_accession in side_accession) or verbose_output=='Y'):
                        side_accession.append(new_side_accession)
                        side_suffix.append(file_number)

        #Construct general accession for this sequence
        general_accession = ">generic|"+main_accession+"_"+main_suffix+"|"+description
        to_add=""
        for i in range(0, len(side_accession)):
            to_add += side_accession[i]+"_db"+str(side_suffix[i])+"#"
        to_add = to_add.strip('#')
        if to_add!="":
            general_accession += " ["+to_add+"]"

        #Save to comb_dict
        comb_data[general_accession] = sequence

    return (comb_data, counts_per_bincode)

#Consturct presence possibilities matrix
def construct_combinations(input_file_dict, overview_file):

    #Init
    #Structure comb_dict: binaire code = comma separated list of files in which present
    comb_dict = defaultdict(lambda: defaultdict())

    #Use itertools to make these codes and lists
    bin_codes = list(itertools.product([0,1], repeat=len(input_file_dict.keys())))
    for i in bin_codes:
        if not sum(i)==0:
            fileset = list(itertools.compress(input_file_dict.values(), i)) #link to overview of files
            i = ''.join(map(str,i)) #Restructure binary code into string
            comb_dict[i]['files'] = fileset

    #Sort combinations
    sort_number=1
    print "Overview of file combinations:"
    print "\t\tBinCode\t\tFiles"
    for i in sorted(comb_dict.keys(), cmp=cmp_files):
        comb_dict[i]['sort_number'] = sort_number
        sort_number+=1
        str_files = ", ".join(comb_dict[i]['files'])
        print str(comb_dict[i]['sort_number'])+":\t\t"+i+"\t\t"+str_files
    print

    return comb_dict

#Cmp for sorting files
def cmp_files(x,y):
    if x.count('1')>y.count('1'):
        return 1
    elif x.count('1')<y.count('1'):
        return -1
    else:
        if x<y:
            return 1
        else:
            return -1

#Read input
def read_input_files(input_file_dict):

    #Init
    input_data = defaultdict(lambda: defaultdict())
    accession = ""
    sequence = ""

    #Read per input fasta file
    for file_number in input_file_dict.keys():
        with open(input_file_dict[file_number], 'r') as FR:
            lines = FR.readlines()
            lines = map(lambda x: x.rstrip("\n"), lines)
            for i in range(0, len(lines)):
                #Check if accession -> new entry
                if((re.search('^>', lines[i]))):
                    if(i!=0):
                        #Save previous entry
                        input_data[file_number][accession] = sequence
                    #Get new entry accession and empty sequence
                    accession = lines[i]
                    sequence = ""
                else:
                    #Concat next line of sequence
                    sequence += lines[i]
            #Save the last entry
            if(sequence!=""):
                input_data[file_number][accession] = sequence
                accession=""
                sequence=""

    return input_data

#Input file dict
def make_input_file_dict(input_files, overview_file):

    #Init
    input_file_dict = defaultdict()

    #Construct dict
    i=1
    for file in input_files:
        input_file_dict[i] = file
        i+=1

    #Print overview to output
    print "Overview of files:"
    print "\t\tFile"
    for i in sorted(input_file_dict.keys()):
        print str(i)+":\t\t"+input_file_dict[i]
    print

    #Print overview to file
    with open(overview_file, 'w') as FW:
        FW.write("Overview of files:\n")
        FW.write("\t\tFile\n")
        for i in sorted(input_file_dict.keys()):
            FW.write(str(i) + ":\t\t" + input_file_dict[i]+"\n")
        FW.write("\n")

    return input_file_dict

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