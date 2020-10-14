import traceback
import time
import argparse
import os
import sys
import re
from multiprocessing import Pool
from functools import partial
import numpy as np
import pandas as pd
import phylopandas as ph

__author__ = 'Steven Verbruggen'
#Execute 'python analyse_proteoforms_percolator.py -h' for more

def main():

    starttime = time.time()

    print()
    print("###################################")
    print("# ORF calling out of RNA-seq data #")
    print("###################################")
    print()

    #Parse arguments from command line
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS,
                                     description="A tool to call ORFs and export a protein FASTA based on a list of called transcripts.",
                                     add_help=False)

    man_args = parser.add_argument_group("Mandatory parameters")
    man_args.add_argument("--called_transcripts", "-t", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="List of called transcripts")
    man_args.add_argument("--transcript_fasta", "-f", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="Path to a FASTA file of all transcripts and their sequences")

    opt_args = parser.add_argument_group("Optional parameters")
    opt_args.add_argument("--help", "-h", action="help", help="Show help message and exit")
    opt_args.add_argument("--workdir", "-w", action="store", required=False, nargs="?", metavar="FOLDER", default="",
                        type=str, help="Working directory (default: CWD)")
    opt_args.add_argument("--output_fasta", "-o", action="store", required=False, nargs="?", metavar="PATH",
                          default="RNA_seq_proteins.fasta", type=str, help="Output FASTA of possible proteins")
    opt_args.add_argument("--cores", "-c", action="store", required=False, nargs="?", metavar="INT", type=int,
                           default=1, help="Amount of cores to use (default: 1)")
    opt_args.add_argument("--min_aa_length", "-m", action="store", required=False, nargs="?", metavar="INTEGER", default=6,
                          type=int, help="Minimum AA length for an ORF to be saved (default: 6)")

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

    #Read in input
    print("Read in the called transcripts")
    print()
    sys.stdout.flush()
    called_transcripts = read_in_called_transcripts_salmon(args.called_transcripts)

    #Search for sequences of the transcripts
    print("Add sequence information to the transcript information")
    print()
    sys.stdout.flush()
    called_transcripts = add_sequences(called_transcripts, args.transcript_fasta)
    #print(called_transcripts)
    #print()

    #Search for all possible ORFs with multiple threads and return the resulting possibly translated proteins
    print("Search for all possible ORFs. Multicore starting.")
    sys.stdout.flush()
    possible_proteins = parallel_search_orfs(called_transcripts, args.cores, args.min_aa_length)
    print("Multicore finished.")
    print()
    sys.stdout.flush()

    #Write output fasta
    print("Write FASTA output.")
    print()
    sys.stdout.flush()
    write_output_fasta(possible_proteins, args.output_fasta)

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

#Write output fasta
def write_output_fasta(possible_proteins, out_fasta):

    #Open file
    with open(out_fasta, 'w') as FW:
        for prot_idx, protein_row in possible_proteins.iterrows():
            header = ">generic|"+protein_row['ORF_ID']+"|"+protein_row['Transcript_ID']+" "+protein_row['StartCodon']
            FW.write(header+"\n")
            FW.write(protein_row['ProteinSeq']+"\n")

    return

#Search for all possible ORFs with multiple threads and return the resulting possibly translated proteins
def parallel_search_orfs(transcripts, cores, min_aa_length):

    #Split the transcripts in N equal chunks based on the available cores
    transcripts_split = np.array_split(transcripts, cores)
    #Start up cores
    pool = Pool(cores)
    #Get the results of the process and concat all dataframes of results per core into one big dataframe
    func = partial(search_orfs, min_aa_length)
    possible_proteins = pd.concat(pool.map(func, transcripts_split))
    #Close multi-processing
    pool.close()
    pool.join()
    #Reset the index of the total results
    possible_proteins = possible_proteins.reset_index(drop=True)

    #print("Full results:")
    #print(possible_proteins)
    #print()

    return possible_proteins

#Search for ORFs (process per core)
def search_orfs(min_aa_length, transcripts):

    #Init
    possible_proteins = []

    #Search for ORFs per sequence and add the results as a list to a nested list of all results for this core
    transcripts.apply(search_orfs_per_seq, axis=1, possible_proteins=possible_proteins, min_aa_length=min_aa_length)
    #Convert the nested list of found ORFs to a dataframe
    possible_proteins = pd.DataFrame(possible_proteins, columns=['ORF_ID', 'Transcript_ID', 'Start', 'End', 'StartCodon', 'SeqLength', 'Sequence',"ProteinSeq"])
    #print("Results of 1 core:")
    #print(possible_proteins)
    #print()

    #Structure output
    #   ORF id (=transcript+start position), transcript id, start position in transcript, end position in transcript, start codon, ORF seq, pep seq
    # Annotation of ORF necessary?
    # Genomic coordinates necessary?


    return possible_proteins

#Search in a transcript sequence for ORFs and add to the list of possible proteins
def search_orfs_per_seq(transcripts_row, possible_proteins, min_aa_length):

    #Rename
    transcript_seq = transcripts_row['Sequence']
    #print(transcript_seq)
    transcript_id = transcripts_row['Name']

    #Init
    deleted_bases = 0

    #Search for possible ORF starts based on NTG start codon
    m = re.search('([ACTG]TG[ACTG]+)',transcript_seq)
    while(m):
        initial_ORF_seq = m.group(1)
        #Search for in-frame stop codon
        final_ORF_seq = ""
        stop_codon_found = False
        #Go over primary ORF sequence per codon
        for i in range(0, len(initial_ORF_seq)-1, 3):
            codon = initial_ORF_seq[i:(i+3)]
            final_ORF_seq+=codon #Add the codon to the tmp ORF sequence
            #Check if codon is a STOP codon
            if (codon=="TAG") or (codon=="TAA") or (codon=="TGA"):
                #Stop the search if a STOP codon is found
                stop_codon_found=True
                break
        #Only add the ORF to the final results if a stop codon was found
        if stop_codon_found==True:
            #Only save the ORF when longer than the min AA length
            if (len(final_ORF_seq)/3)>=min_aa_length:
                start_position = deleted_bases + m.start()+1
                stop_position = start_position + len(final_ORF_seq) - 1
                ORF_id = transcript_id+"_"+str(start_position)
                start_codon = final_ORF_seq[0:3]
                aa_seq = translate_seq(final_ORF_seq, transcript_id)
                found_ORF = [ORF_id,transcript_id, start_position, stop_position, start_codon, len(final_ORF_seq), final_ORF_seq, aa_seq]
                possible_proteins.append(found_ORF)
        #Rest of the sequence to search further in (minus the first base to not find the same ORF twice)
        deleted_bases = deleted_bases + (len(transcript_seq) - len(initial_ORF_seq)) + 1
        transcript_seq = initial_ORF_seq[1:]
        m = re.search('([ACTG]TG[ACTG]+)', transcript_seq) #Search again

    return possible_proteins

## Translate sequence to amino acids
def translate_seq(tr_seq, transcript_id):

    #Init
    aa_seq = ""

    # The codon table represents the corresponding amino acid for every codon.
    codonTable = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
                  'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
                  'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                  'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                  'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                  'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                  'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                  'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
                  'AAN': 'X', 'ATN': 'X', 'ACN': 'T', 'AGN': 'X', 'ANN': 'X', 'ANA': 'X', 'ANG': 'X', 'ANC': 'X',
                  'ANT': 'X', 'TAN': 'X', 'TTN': 'X', 'TCN': 'S', 'TGN': 'X', 'TNN': 'X', 'TNA': 'X', 'TNG': 'X',
                  'TNC': 'X', 'TNT': 'X', 'CAN': 'X', 'CTN': 'L', 'CCN': 'P', 'CGN': 'R', 'CNN': 'X', 'CNA': 'X',
                  'CNG': 'X', 'CNC': 'X', 'CNT': 'X', 'GAN': 'X', 'GTN': 'V', 'GCN': 'A', 'GGN': 'G', 'GNN': 'X',
                  'GNA': 'X', 'GNG': 'X', 'GNC': 'X', 'GNT': 'X', 'NAN': 'X', 'NAA': 'X', 'NAG': 'X', 'NAC': 'X',
                  'NAT': 'X', 'NTN': 'X', 'NTA': 'X', 'NTG': 'X', 'NTC': 'X', 'NTT': 'X', 'NGN': 'X', 'NGA': 'X',
                  'NGG': 'X', 'NGC': 'X', 'NGT': 'X', 'NCN': 'X', 'NCA': 'X', 'NCG': 'X', 'NCC': 'X', 'NCT': 'X',
                  'NNN': 'X', 'NNA': 'X', 'NNG': 'X', 'NNC': 'X', 'NNT': 'X'}

    #Go over sequence per triplet
    codon=""
    for pos in range(0, len(tr_seq)-2, 3):
        codon = tr_seq[pos:pos+3]
        aa_seq+=codonTable[codon]
        #Control if stop codon is only at the end
        if codonTable[codon]=='*':
            if pos!=(len(tr_seq)-3):
                print("Error: stop earlier than expected (transcript ID: "+transcript_id+")")
                aa_seq = "early_stop,"+aa_seq
                print("tr seq: "+tr_seq)
                print("AA seq: "+aa_seq)
                sys.stdout.flush()
                return aa_seq
    #Check if end codon is a STOP codon
    if codon!='TAA' and codon!='TAG' and codon!='TGA':
        print("Error: no stop codon (TAA, TAG or TGA) found at the end (transcript ID: "+transcript_id+")")
        print("tr seq: " + tr_seq)
        print("AA seq: " + aa_seq)
        sys.stdout.flush()
        aa_seq="no_stop"
        return aa_seq
    #Check for near-cognate starts and replace the first amino acid for methionine then
    if re.search('[ACTG]TG|A[ACTG]G|AT[ACTG]', tr_seq[:3]):
        aa_seq = "M" + aa_seq[1:]
    #Clip off the stop codon symbol
    m = re.search('^(\w+)\*$', aa_seq)
    if m:
        aa_seq = m.group(1)
    #print "tr seq: "+tr_seq
    #print "aa seq: "+aa_seq

    return aa_seq

#Add sequences to transcript dataframe
def add_sequences(transcripts, tr_seqs_file):

    #Read in sequence fasta
    tr_seqs = read_in_seqs(tr_seqs_file)
    #Merge called transcripts with transcript sequence data
    transcripts = transcripts.merge(tr_seqs, how="left", left_on="Name", right_on="id")
    #Drop unnecessary columns and rename sequence column header
    transcripts = transcripts.drop(labels=['id','description','label','uid'], axis=1).rename(columns={'sequence': 'Sequence'})
    #Check if Length is as long as length of the sequence: extra control
    entries_with_unsatisfying_lengths = transcripts[transcripts['Length'].astype(int) != transcripts['Sequence'].str.len()].shape[0]
    if(entries_with_unsatisfying_lengths!=0):
        print("WARNING: "+str(entries_with_unsatisfying_lengths)+" transcript entries with unmatching sequence lengths. Further investigation required.\n")

    return transcripts

#Transcript sequences to dataframe
def read_in_seqs(in_file):
    seqs = ph.read_fasta(in_file)
    return seqs

#Called transcripts to dataframe
def read_in_called_transcripts_salmon(in_file):
    called_transcripts = pd.read_csv(in_file, sep=",", header=0)
    return called_transcripts

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
###################################