#!/usr/bin/env python
import traceback
import sys
import os
import time
import argparse
from collections import defaultdict
import re
from collections import defaultdict


__author__ = 'Steven Verbruggen'

def main():
    
    starttime = time.time()

    print
    print "##########################"
    print "# CNE fasta converter    #"
    print "##########################"
    print

    #Parse arguments from command line
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS,
                                     description="This tool is part of the PROTEOFORMER pipeline. It converts CNE fasta files.",
                                     add_help=False)

    man_args = parser.add_argument_group("Mandatory parameters")
    man_args.add_argument("--input", "-i", action="store", required=True, nargs="?", metavar="fasta", type=str,help="Input fasta")
    man_args.add_argument("--output", "-o", action="store", required=True, nargs="?", metavar="fasta", type=str,help="Output fasta")

    opt_args = parser.add_argument_group("Optional parameters")
    opt_args.add_argument("--help", "-h", action="help", help="Show help message and exit")
    opt_args.add_argument("--workdir", "-w", action="store", required=False, nargs="?", metavar="FOLDER", default="",type=str, help="Working directory (default: CWD)")

    args = parser.parse_args()

    #default workdir is CWD, default tmp, default fdr
    if args.workdir == "":
        args.workdir = os.getcwd()

    ########
    # MAIN #
    ########
    
    input_data = read_fasta(args.input)
    
    output_data = restruc_fasta(input_data)
    
    write_fasta(output_data, args.output)
    
    #End of program message
    print
    print "-----------------------"
    print "[%s]: PROGRAM COMPLETE" % (convert_time(time.time()-starttime))
    print "-----------------------"
    print
    sys.stdout.flush()
    
    
    return

#Write fasta
def write_fasta(output_data, output_path):

    with open(output_path, 'w') as FW:
        for acc in sorted(output_data.keys()):
            line = acc+"\n"+output_data[acc]+"\n"
            FW.write(line)
    
    return

#Restruc fasta
def restruc_fasta(input_data):

    output_data = defaultdict()

    for acc in input_data.keys():
        m=re.search('^>(CNAG.+)\W\|\W(CNAG.+)\W\|\W(.+)$', acc)
        if m:
            tr_id = m.group(1)
            gene_id = m.group(2)
            descr = m.group(3)
            new_acc = ">sp|"+tr_id+"|"+gene_id+" "+descr
            seq = input_data[acc].rstrip('*')
            output_data[new_acc]=seq
    
    return output_data

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

