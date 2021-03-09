#!/usr/bin/env python
import traceback
import sys
import os
import time
import argparse
import re
from collections import defaultdict
import itertools
#from print_collection import print_dict

__author__ = 'Steven Verbruggen'
#Execute 'python combine_dbs.py -h' for more information

def main():

    starttime = time.time()

    print
    print "##################"
    print "# Split FASTA DB #"
    print "##################"
    print
    print "This program is part of the PROTEOFORMER pipeline"
    print

    #Parse arguments from command line
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS,
                                     description="This tool is part of the PROTEOFORMER pipeline. It splits "
                                                 "a fasta file in equal parts.",
                                     add_help=False)

    man_args = parser.add_argument_group("Mandatory parameters")
    man_args.add_argument("--input_file", "-i", action="store", required=True, nargs="?", metavar="PATH",
                          default="", type=str, help="Input fasta file (mandatory) ")
    man_args.add_argument("--parts", "-p", action="store", required=True, nargs="?", metavar="INTEGER",
                          default="", type=int, help="Number of parts to divide the input in (mandatory)")

    opt_args = parser.add_argument_group("Optional parameters")
    opt_args.add_argument("--help", "-h", action="help", help="Show help message and exit")
    opt_args.add_argument("--workdir", "-w", action="store", required=False, nargs="?", metavar="FOLDER", default="",
                        type=str, help="Working directory (default: CWD)")
    opt_args.add_argument("--output_file_prefix", "-o", action="store", required=False, nargs="?", metavar="PATH",
                          default="output_", type=str, help="Prefix for the output files (default: output_)")

    args = parser.parse_args()

    #default workdir is CWD
    if args.workdir == "":
        args.workdir = os.getcwd()
    #Check mandatory input
    if args.input_file == "":
        print "Do not forget to specify the input fasta file!"
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

    (sequences, accessions) = read_input_file(args.input_file)
    total_number_of_sequences = len(accessions)
    print("Input: "+str(total_number_of_sequences)+" records")
    sys.stdout.flush()

    split_in_files(sequences, accessions, args.parts, args.output_file_prefix)


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

#Split in files
def split_in_files(sequences, accessions, parts, output_file_prefix):

    first_accession_index = 0
    for i in range(1, parts):
        last_accession_index = int(len(accessions)*1.0/parts*i)
        select_accessions = accessions[first_accession_index:last_accession_index]
        print_to_output(sequences, select_accessions, output_file_prefix, i)
        print("File "+str(i)+": "+str(len(select_accessions))+" entries written")
        sys.stdout.flush()
        first_accession_index=last_accession_index
    #Also do it for the last part
    last_accession_index=len(accessions)
    select_accessions = accessions[first_accession_index:last_accession_index]
    print_to_output(sequences, select_accessions, output_file_prefix, parts)
    print("File "+str(parts)+": "+str(len(select_accessions))+" entries written")
    sys.stdout.flush()

    return

#Print to output fasta file
def print_to_output(sequences, accessions, prefix, index):

    output_file = prefix+str(index)+".fa"

    with open(output_file, 'w') as FW:
        for accession in accessions:
            FW.write(accession+"\n")
            FW.write(sequences[accession]+"\n")

    return

#Read input
def read_input_file(input_file):

    #Init
    input_data = defaultdict()
    accessions = []
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
                    accessions.append(accession)
                #Get new entry accession and empty sequence
                accession = lines[i]
                sequence = ""
            else:
                #Concat next line of sequence
                sequence += lines[i]
        #Save the last entry
        if(sequence!=""):
            input_data[accession] = sequence
            accessions.append(accession)
            accession=""
            sequence=""

    return(input_data, accessions)

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