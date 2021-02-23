import traceback
import sys
import os
import time
import argparse
import re
from collections import defaultdict
import itertools

__author__ = 'Steven Verbruggen'
#Execute 'python combine_with_uniprot.py -h' for more information

def main():

    starttime = time.time()

    print
    print "##########################"
    print "# Filter uncanonical AAs #"
    print "##########################"
    print
    print "This program is part of the PROTEOFORMER pipeline"
    print

    #Parse arguments from command line
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS,
                                     description="This tool is part of the PROTEOFORMER pipeline. It filters "
                                                 "a (combined) fasta file of PROTEOFORMER for sequences with only the 20 canonical amino acids.",
                                     add_help=False)

    man_args = parser.add_argument_group("Mandatory parameters")
    man_args.add_argument("--fasta", "-f", action="store", required=True, nargs="?", metavar="PATH",
                          default="", type=str, help="Input (combined) fasta file (mandatory) ")

    opt_args = parser.add_argument_group("Optional parameters")
    opt_args.add_argument("--help", "-h", action="help", help="Show help message and exit")
    opt_args.add_argument("--workdir", "-w", action="store", required=False, nargs="?", metavar="FOLDER", default="",
                        type=str, help="Working directory (default: CWD)")
    opt_args.add_argument("--output_fasta", "-o", action="store", required=False, nargs="?", metavar="PATH",
                          default="filtered_db.fasta", type=str, help="Output fasta file "
                                                                                            "(default: filtered_db.fasta)")

    args = parser.parse_args()

    #default workdir is CWD
    if args.workdir == "":
        args.workdir = os.getcwd()
    if args.output_fasta == "filtered_db.fasta":
        args.output_fasta = args.workdir + "/" + args.output_fasta

    #List parameters
    print "Parameters:"
    for arg in vars(args):
        print '    %-15s\t%s' % (arg, getattr(args, arg))
    print
    sys.stdout.flush()

    ########
    # MAIN #
    ########

    #Read in fasta data
    input_data = read_fasta(args.fasta)
    print("Input sequences: "+str(len(input_data)))

    #Filter Fasta
    filtered_data = filter_data(input_data)

    #Write output file
    write_output(filtered_data, args.output_fasta)
    print("Output sequences: "+str(len(filtered_data)))



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

#Filter data
def filter_data(input_data):

    #Init
    output_data = defaultdict()

    filtered_out = 0
    for acc in input_data:
        if not re.search('[^RHKDESTNQCGPAVILMFYW]',input_data[acc]):#If not anything present except for the 20 canonical AAs
            output_data[acc] = input_data[acc]
        else:
            filtered_out+=1

    print("Sequences filtered out: "+str(filtered_out))

    return output_data

#Write output fasta file
def write_output(out_data, fasta):

    with open(fasta, 'w') as FW:
        for acc in out_data:
            FW.write(acc+"\n")
            FW.write(out_data[acc]+"\n")

    return

#read fasta
def read_fasta(fasta):

    input_data=defaultdict()

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
                    input_data[accession] = sequence
                # Get new entry accession and empty sequence
                accession = lines[i]
                sequence = ""
            else:
                # Concat next line of sequence
                sequence += lines[i]
        # Save the last entry
        if (sequence != ""):
            input_data[accession] = sequence
            accession = ""
            sequence = ""

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