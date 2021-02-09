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
    print("##########################################")
    print("# Parse transcripts coming out of Salmon #")
    print("##########################################")
    print()

    #Parse arguments from command line
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS,
                                     description="A tool to parse the transcripts coming out of Salmon",
                                     add_help=False)

    man_args = parser.add_argument_group("Mandatory parameters")
    man_args.add_argument("--salmon_sf", "-s", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="Salmon output transcript abundance file")
    man_args.add_argument("--ambig_info", "-a", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="Ambiguous count info file from Salmon")

    opt_args = parser.add_argument_group("Optional parameters")
    opt_args.add_argument("--help", "-h", action="help", help="Show help message and exit")
    opt_args.add_argument("--workdir", "-w", action="store", required=False, nargs="?", metavar="FOLDER", default="",
                        type=str, help="Working directory (default: CWD)")
    opt_args.add_argument("--transcript_output_file", "-x", action="store", required=False, nargs="?", metavar="PATH",
                          type=str, default="called_transcripts_salmon.csv",
                          help="CSV output file of called transcript (default: called_transcripts_salmon.csv)")
    opt_args.add_argument("--transcript_output_file_merged", "-m", action="store", required=False, nargs="?", metavar="PATH",
                          type=str, default="called_transcripts_mergedCounts.csv",
                          help="CSV output file of called transcript (default: called_transcripts_mergedCounts.csv)")

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
    salmon_data = pd.read_csv(args.salmon_sf, sep="\t", header=0)
    total_transcript_count = salmon_data.shape[0]
    ambig_data = pd.read_csv(args.ambig_info, sep="\t", header=0)
    merged_data = salmon_data.merge(ambig_data, left_index=True, right_index=True)

    #Select transcripts with Salmon read count>0
    selected_salmon_data = salmon_data[salmon_data['NumReads']>0]
    selected_transcript_count = selected_salmon_data.shape[0]
    print("Based on Salmon read counts:")
    print(str(selected_transcript_count)+" transcripts of the "+str(total_transcript_count)+" total input transcripts have a read count>0 and were selected.\n")

    #Select transcripts with ambig read count>0
    merged_data['TotalCount'] = merged_data['UniqueCount'] + merged_data['AmbigCount']
    selected_merged_data = merged_data[merged_data['TotalCount']>0]
    selected_merged_count = selected_merged_data.shape[0]
    print("Based on total raw counts (unique + ambigue counts):")
    print(str(selected_merged_count)+" transcripts of the "+str(total_transcript_count)+" total input transcripts have a read count>0 and were selected.\n")

    #Ouput transcript list
    selected_salmon_data.to_csv(args.transcript_output_file, index=False, header=True)
    selected_merged_data.to_csv(args.transcript_output_file_merged, index=False, header=True)

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

