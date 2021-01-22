import os
import sys
import re
import traceback
import time
import argparse
import numpy as np
import pandas as pd
from collections import defaultdict

__author__ = 'Steven Verbruggen'

def main():

    starttime = time.time()

    print()
    print("################################")
    print("# Refactor tab file for Prosit #")
    print("################################")
    print()

    #Parse arguments from command line
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS,
                                     description="This tool refurbishes a MS2ReScore tab file: M(ox) vocab and protein info.",
                                     add_help=False)

    man_args = parser.add_argument_group("Mandatory parameters")
    man_args.add_argument("--input_pin", "-i", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="The input ms2rescore pin file (mandatory)")

    opt_args = parser.add_argument_group("Optional parameters")
    opt_args.add_argument("--help", "-h", action="help", help="Show help message and exit")
    opt_args.add_argument("--output_pin", "-o", action="store", required=False, nargs="?", metavar="PATH", default="ms2rescore_rearr.pin",
                          type=str, help="The output file: rearranged pin file (default: ms2rescore_rearr.pin)")

    args = parser.parse_args()

    #List parameters
    print("Parameters:")
    for arg in vars(args):
        print('    %-15s\t%s' % (arg, getattr(args, arg)))
    print()
    sys.stdout.flush()
    
    ########
    # MAIN #
    ########
    
    #Parse line per line
    with open(args.output_pin, 'w') as FW:
        with open(args.input_pin, 'r') as FR:
            #Read in input file
            lines = FR.readlines()

            #Parse header
            header_line = lines[0].rstrip()
            FW.write(header_line+"\n")
            header_cols = header_line.split("\t")
            number_of_cols = len(header_cols)

            for line in lines[1:]:
                line = line.rstrip()
                line_vals = line.split("\t")

                #Parse everything except protein names
                normal_vals = line_vals[0:number_of_cols-1]
                normal_vals[-1]=re.sub("M\(ox\)","M", normal_vals[-1])
                normal_vals[-1]=re.sub("^\_","_.", normal_vals[-1])
                normal_vals[-1]=re.sub("\_$","._", normal_vals[-1])
                output_line = "\t".join(normal_vals)
                output_line=output_line+"\t"

                #Parse protein names
                protein_vals = line_vals[number_of_cols-1:]
                i=0
                for protein_val in protein_vals:
                    protein_vals[i]=re.sub("REV_", "REV__", protein_val)
                    i+=1
                protein_field = ";".join(protein_vals)
                output_line=output_line+protein_field+"\n"
                FW.write(output_line)


    
    

    #Refactor proteins with semicolon instead of tab
    #Refactor methionine ox
    
    #Export to new tab file
    #ms2rescore_pin.to_csv(args.output_tab, sep="\t", header=True, index=False)
    
    # End of program message
    print()
    print("-----------------------")
    print("[%s]: PROGRAM COMPLETE" % (convert_time(time.time() - starttime)))
    print("-----------------------")
    print()
    sys.stdout.flush()

    return

##########
#  SUBS  #
##########



## Data Dumper for recursively printing nested dictionaries and defaultDicts ##
def print_dict(dictionary, indent='', braces=0):
    """
    :param dictionary: dict or defaultdict to be printed
    :param indent: indentation
    :param braces:
    :return: void
    """

    for key, value in dictionary.items():
        if isinstance(value, dict):
            print('%s%s%s%s' % (indent, braces * '[', key, braces * ']'))
            print_dict(value, indent + '  ', braces)
        else:
            print(indent + '%s = %s' % (key, value))
    return

def convert_time(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)

###### DIRECT TO MAIN ############
if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        traceback.print_exc()
##################################