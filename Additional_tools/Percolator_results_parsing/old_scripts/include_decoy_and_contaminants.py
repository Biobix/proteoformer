#!/usr/bin/env python
import traceback
import sys
import os
import time
import argparse
import re
from collections import defaultdict

__author__ = 'Steven Verbruggen'
#Execute 'python include_decoy_and_contaminants.py -h' for more information

def main():

    starttime = time.time()

    print
    print "#####################################################"
    print "# Include decoys and contaminants for Percolator PI #"
    print "#####################################################"
    print
    print "This program is part of the PROTEOFORMER pipeline"
    print

    #Parse arguments from command line
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS,
                                     description="This tool prepares a FASTA file for Percolator protein inference."
                                                 "It combines the sequence search space with decoys and contaminants.",
                                     add_help=False)

    man_args = parser.add_argument_group("Mandatory parameters")
    man_args.add_argument("--in_fasta", "-i", action="store", required=True, nargs="?", metavar="fasta",
                          default="", type=str, help="Input fasta file (mandatory) ")
    man_args.add_argument("--cont", "-c", action="store", required=True, nargs="?", metavar="fasta",
                          default="", type=str, help="Contaminant fasta file (mandatory)")

    opt_args = parser.add_argument_group("Optional parameters")
    opt_args.add_argument("--help", "-h", action="help", help="Show help message and exit")
    opt_args.add_argument("--workdir", "-w", action="store", required=False, nargs="?", metavar="FOLDER", default="",
                        type=str, help="Working directory (default: CWD)")
    opt_args.add_argument("--out_fasta", "-o", action="store", required=False, nargs="?", metavar="fasta",
                          default="decoy_and_cont_incl.fasta", type=str, help="Output fasta file (default: decoy_and_cont_incl.fasta)")
    opt_args.add_argument("--out_fasta_beforeKR", "-k", action="store", required=False, nargs="?", metavar="fasta",
                          default="decoy_beforeKR_cont_incl.fasta", type=str, help="Output fasta file reversed before K or R (default: decoy_beforeKR_cont_incl.fasta)")

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

    #Load fasta file
    input_seqs = load_fasta(args.in_fasta)

    #Create decoys based on trypsinP
    #TrypsinP splits every K and R after this AA
    decoys = create_decoys(input_seqs)
    decoys_before_KR = create_decoys_beforeKR(input_seqs)

    #Load contaminants
    contaminants = load_contaminants(args.cont)

    #Write total output fasta file
    write_output(args.out_fasta, input_seqs, decoys, contaminants)
    write_output(args.out_fasta_beforeKR, input_seqs, decoys_before_KR, contaminants)

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

def write_output(out, input_seqs, decoys, contaminants):

    with open(out, 'w') as FW:
        for accession in input_seqs:
            FW.write(accession+"\n")
            FW.write(input_seqs[accession]+"\n")
        for accession in decoys:
            FW.write(accession+"\n")
            FW.write(decoys[accession]+"\n")
        for accession in contaminants:
            FW.write(accession+"\n")
            FW.write(contaminants[accession]+"\n")

    return

def load_contaminants(cont_file):

    #Init
    contaminants = defaultdict()

    with open(cont_file, 'r') as FR:
        lines = FR.readlines()

        accession=''
        sequence=''
        for i in range(0, len(lines)):
            lines[i] = lines[i].rstrip('\n')
            lines[i] = lines[i].rstrip('\r')
            m = re.search('^>(.+)$', lines[i])
            if m:
                if accession!='' and sequence!='':
                    contaminants[accession] = sequence
                    accession=''
                    sequence=''
                accession = ">CON__"+m.group(1)
            else:
                sequence += lines[i]
        #Last sequence
        if accession != '' and sequence != '':
            contaminants[accession] = sequence
            accession=''
            sequence=''

    return contaminants

def create_decoys_beforeKR(input_seqs):

    #Init
    decoys = defaultdict()

    for accession in input_seqs.keys():
        #Make accession
        decoy_accession = ''
        m = re.search('^>(.+)$', accession)
        if m:
            decoy_accession = ">REV__"+m.group(1)

        # Split based on K and R
        cur_fragment = ''
        fragments = []
        for AA in input_seqs[accession]:
            if (AA == 'K' or AA == 'R'):
                fragments.append(cur_fragment) #Add to fragments without K or R
                fragments.append(AA) #Add K or R as separate fragment
                cur_fragment = ''
            else:
                cur_fragment += AA
        # add last fragment
        if cur_fragment != '':
            fragments.append(cur_fragment)

        # reverse fragments
        rev_fragments = []
        for frag in fragments:
            rev_fragments.append(frag[::-1])


        #Paste together to decoy sequence
        decoy_seq = ''
        for frag in rev_fragments:
            decoy_seq+=frag

        # Save
        decoys[decoy_accession] = decoy_seq

    return decoys

def create_decoys(input_seqs):

    #Init
    decoys = defaultdict()

    for accession in input_seqs.keys():
        #Make accession
        decoy_accession = ''
        m = re.search('^>(.+)$', accession)
        if m:
            decoy_accession = ">REV__"+m.group(1)

        #Split based on K and R
        cur_fragment = ''
        fragments = []
        for AA in input_seqs[accession]:
            cur_fragment += AA #add to current fragment
            #Split if necessary
            if (AA=='K' or AA=='R'):
                fragments.append(cur_fragment)
                cur_fragment = ''
        #add last fragment
        if cur_fragment != '':
            fragments.append(cur_fragment)

        #reverse fragments
        rev_fragments = []
        for frag in fragments:
            rev_fragments.append(frag[::-1])

        #Paste together to decoy sequence
        decoy_seq = ''
        for frag in rev_fragments:
            decoy_seq+= frag

        #Save
        decoys[decoy_accession]=decoy_seq

    return decoys

def load_fasta(in_fasta):

    #Init
    input_seqs = defaultdict()

    with open(in_fasta, 'r') as FR:
        lines = FR.readlines()
        accession = ''

        for i in range(0, len(lines)):
            lines[i] = lines[i].rstrip('\n')
            m = re.search('^>.+$', lines[i])
            if m:
                accession = lines[i]
            else:
                seq = lines[i]
                input_seqs[accession] = seq
                accession = ''

    return input_seqs

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