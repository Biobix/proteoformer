import traceback
import sys
import os
import time
import argparse
import re
from collections import defaultdict

__author__ = 'Steven Verbruggen'
#Execute 'python convert_hspvdb.py -h' for more information


def main():

    starttime = time.time()

    print
    print "########################"
    print "# Convert HSPVdb fasta #"
    print "########################"
    print
    print "This program is part of the PROTEOFORMER pipeline"
    print

    #Parse arguments from command line
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS,
                                     description="This tool is part of the PROTEOFORMER pipeline. It converts "
                                                 "a fasta file from HSPVdb to be compatible with PROTEOFORMER.",
                                     add_help=False)

    man_args = parser.add_argument_group("Mandatory parameters")
    man_args.add_argument("--input", "-i", action="store", required=True, nargs="?", metavar="PATH",
                          default="", type=str, help="Input HSPVdb fasta file (mandatory) ")

    opt_args = parser.add_argument_group("Optional parameters")
    opt_args.add_argument("--help", "-h", action="help", help="Show help message and exit")
    opt_args.add_argument("--workdir", "-w", action="store", required=False, nargs="?", metavar="FOLDER", default="",
                        type=str, help="Working directory (default: CWD)")
    opt_args.add_argument("--output_fasta", "-o", action="store", required=False, nargs="?", metavar="PATH",
                          default="conv_hspvdb.fasta", type=str, help="Output converted fasta file "
                          "(default: comb_fasta_uniprot_proteoformer.fasta)")

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

    #Read input data
    input_data = read_fasta(args.input)

    #Convert data
    conv_data = convert_data(input_data)

    #Output converted data
    write_fasta(conv_data, args.output_fasta)

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

#Write data to fasta output
def write_fasta(data, fasta):

    with open(fasta, 'w') as FW:
        for acc in data.keys():
            FW.write(acc+"\n")
            FW.write(data[acc]+"\n")
    return

#Convert data to a more proteoformer-like format
def convert_data(input_data):

    #Init
    conv_data = defaultdict()

    #Convert per entry
    for acc in input_data.keys():
        m = re.search('^>(\S+) \|(.+)$', acc)
        if m:
            #Catch accession parts
            hspv_id = 'hspv'+str(m.group(1))
            #print(hspv_id)
            rest_of_acc = str(m.group(2))
            #print(rest_of_acc)

            #Parse and clean up the last part of the accession
            elements = rest_of_acc.split('|')
            if len(elements)==7:
                elements.append('') #Make that all entries have the same amount of elements
            elements = map(str.strip, elements) #Remove leading and trailing whitespaces
            if elements[3]=='..':
                elements[3]=''#Clear peptides without coordinates
            elements = filter(None, elements)
            #print(elements)

            #Construct new accession structure
            acc_info = ' '.join(elements)
            new_acc = ">generic|"+hspv_id+"|"+acc_info

            #Save
            conv_data[new_acc] = input_data[acc]

    return conv_data

#read fasta
def read_fasta(fasta):

    #Init
    input_data = defaultdict()

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


