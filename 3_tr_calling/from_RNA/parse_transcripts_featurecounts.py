import traceback
import time
import argparse
import os
import sys
import pandas as pd

__author__ = 'Steven Verbruggen'
#Execute 'python analyse_proteoforms_percolator.py -h' for more

def main():

    starttime = time.time()

    print()
    print("#################################################")
    print("# Parse transcripts coming out of featureCounts #")
    print("#################################################")
    print()

    #Parse arguments from command line
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS,
                                     description="A tool to parse the transcripts coming out of featureCounts",
                                     add_help=False)

    man_args = parser.add_argument_group("Mandatory parameters")
    man_args.add_argument("--fc_tab", "-f", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="featureCounts output transcript abundance file")

    opt_args = parser.add_argument_group("Optional parameters")
    opt_args.add_argument("--help", "-h", action="help", help="Show help message and exit")
    opt_args.add_argument("--workdir", "-w", action="store", required=False, nargs="?", metavar="FOLDER", default="",
                        type=str, help="Working directory (default: CWD)")
    opt_args.add_argument("--transcript_output_file", "-x", action="store", required=False, nargs="?", metavar="PATH",
                          type=str, default="called_transcript_fc.csv",
                          help="CSV output file of called transcript (default: called_transcripts_fc.csv)")

    args = parser.parse_args()

    #default workdir is CWD
    if args.workdir == "":
        args.workdir = os.getcwd()

    #List parameters
    print("Parameters:")
    for arg in vars(args):
        print('    %-25s\t%s' % (arg, getattr(args, arg)))
    print()
    sys.stdout.flush()


    ########
    # MAIN #
    ########

    #Read data
    fc_data = pd.read_csv(args.fc_tab, sep="\t", header=0)
    total_transcript_count = fc_data.shape[0]

    #Select transcripts with read count>0
    selected_fc_data = fc_data[fc_data.iloc[:,-1]>0]
    selected_transcript_count = selected_fc_data.shape[0]
    print(str(selected_transcript_count)+" transcripts of the "+str(total_transcript_count)+" total input transcripts have a read count>0 and were selected.\n")

    #Ouput transcript list
    selected_fc_data.to_csv(args.transcript_output_file, index=False, header=True)

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

