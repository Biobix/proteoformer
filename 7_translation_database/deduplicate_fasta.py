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
    print "########################"
    print "# Deduplicate FASTA DB #"
    print "########################"
    print
    print "This program is part of the PROTEOFORMER pipeline"
    print

    #Parse arguments from command line
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS,
                                     description="This tool is part of the PROTEOFORMER pipeline. It deduplicates "
                                                 "a fasta file, i.e. combining full sequence exact matches as 1 record."
                                                 "This can be important for file combination to work with correct numbers.",
                                     add_help=False)

    man_args = parser.add_argument_group("Mandatory parameters")
    man_args.add_argument("--input_files", "-i", action="store", required=True, nargs="?", metavar="COMMA-SEPARATED LIST",
                          default="", type=str, help="Comma-separated list of all paths to all input fasta files (mandatory) ")

    opt_args = parser.add_argument_group("Optional parameters")
    opt_args.add_argument("--help", "-h", action="help", help="Show help message and exit")
    opt_args.add_argument("--workdir", "-w", action="store", required=False, nargs="?", metavar="FOLDER", default="",
                        type=str, help="Working directory (default: CWD)")
    opt_args.add_argument("--output_files", "-o", action="store", required=False, nargs="?", metavar="PATH",
                          default="dedupl", type=str, help="Comma-separated list of the output fasta files"
                                                                           "(default: *_dedupl.fa)")

    args = parser.parse_args()

    #default workdir is CWD
    if args.workdir == "":
        args.workdir = os.getcwd()
    #Check mandatory input
    if args.input_files == "":
        print "Do not forget to specify all input fasta files!"
        sys.exit()

    input_files = re.split(',', args.input_files)
    output_files = []
    if args.output_files=="dedupl":
        for i in input_files:
            m_input = re.search('^(.+)\.fa', i)
            if m_input:
                output_file = m_input.group(1)+"_dedupl.fa"
                output_files.append(output_file)
    else:
        output_files = re.split(',', args.output_files)
    if len(input_files)!=len(output_files):
        print("Comma-seperated lists of input and output files has not the same number of elements!")
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

    for i in range(0,len(input_files)):
        input_file = input_files[i]
        output_file = output_files[i]
        print ("File: "+input_file)
        sys.stdout.flush()
        input_data = read_input_file(input_file)
        print("Input: "+str(len(input_data))+" records")
        sys.stdout.flush()
        output_data = deduplicate(input_data)
        print("Output: " + str(len(output_data)) + " records")
        sys.stdout.flush()
        print_to_output(output_data, output_file)
        print("Output written to "+output_file)
        print("\n")
        sys.stdout.flush()


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

#Deduplicate fasta data
def deduplicate(input_data):

    #Init
    output_data = defaultdict()
    per_seq_data = defaultdict()

    #Deduplicate on sequence level
    for acc in input_data:
        seq = input_data[acc]
        if seq not in per_seq_data:
            per_seq_data[seq] = acc
        else:
            m_id = re.search('^>(generic|hspv|sp|tr)\|(.+?)\|.+', acc)
            if m_id:
                id = m_id.group(2)
                updated_acc = per_seq_data[seq]+" dupl_"+id
                per_seq_data[seq] = updated_acc

    #Construct output data
    for seq in per_seq_data:
        acc = per_seq_data[seq]
        output_data[acc] = seq

    return output_data

#Print to output fasta file
def print_to_output(output_data, output_file):

    with open(output_file, 'w') as FW:
        for accession in output_data:
            FW.write(accession+"\n")
            FW.write(output_data[accession]+"\n")

    return

#Read input
def read_input_file(input_file):

    #Init
    input_data = defaultdict()
    accession = ""
    sequence = ""

    #Read fasta file
    with open(input_file, 'r') as FR:
        lines = FR.readlines()
        lines = map(lambda x: x.rstrip("\n"), lines)
        for i in range(0, len(lines)):
            #Check if accession -> new entry
            if((re.search('^>', lines[i]))):
                if(i!=0):
                    #Save previous entry
                    input_data[accession] = sequence
                #Get new entry accession and empty sequence
                accession = lines[i]
                sequence = ""
            else:
                #Concat next line of sequence
                sequence += lines[i]
        #Save the last entry
        if(sequence!=""):
            input_data[accession] = sequence
            accession=""
            sequence=""

    return input_data

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