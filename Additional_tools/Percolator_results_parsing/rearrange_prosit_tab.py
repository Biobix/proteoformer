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
                                     description="This tool refurbishes a Prosit tab file: M(ox) vocab and protein info.",
                                     add_help=False)

    man_args = parser.add_argument_group("Mandatory parameters")
    man_args.add_argument("--input_tab", "-i", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="The input prosit tab file (mandatory)")

    opt_args = parser.add_argument_group("Optional parameters")
    opt_args.add_argument("--help", "-h", action="help", help="Show help message and exit")
    opt_args.add_argument("--msmstxt", "-m", action="store", required=False, nargs="?", metavar="PATH", default="msms.txt",
                         type=str, help="MaxQuant msms.txt file (default:msms.txt)")
    opt_args.add_argument("--proteingroups", "-p", action="store", required=False, nargs="?", metavar="PATH", default="proteinGroups.txt",
                         type=str, help="MaxQuant proteinGroups.txt file (default proteinGroups.txt)")
    opt_args.add_argument("--fasta", '-f', action="store", required=False, nargs="?", metavar="PATH", default='proteoformer.fasta',
                         type=str, help="Fasta file with decoys and contaminants (default: proteoformer.fasta)")
    opt_args.add_argument("--output_tab", "-o", action="store", required=False, nargs="?", metavar="PATH", default="prosit_rearr.tab",
                          type=str, help="The output file: rearranged tab file (default: prosit_rearr.tab)")

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

    #Load in prosit tab file
    print("Reformat M(ox) vocab")
    print()
    sys.stdout.flush()
    prosit_tab = pd.read_csv(args.input_tab, sep="\t")
    
    #Split PSM ID
    prosit_tab = split_prosit_psmid(prosit_tab)
    
    #Convert methionines
    prosit_tab = methionine_convert(prosit_tab)
    
    #Reconstruct the PSM ID
    prosit_tab = reconstruct_psmid(prosit_tab)
    
    #Load maxquant files
    print("Loading MaxQuant files")
    print()
    sys.stdout.flush()
    msmstxt, proteinGroups = load_mq_files(args.msmstxt, args.proteingroups)
    
    #Load fasta file
    print("Loading Fasta file")
    print()
    sys.stdout.flush()
    fasta = load_fasta(args.fasta)
    
    #Get protein info
    print("Get protein info")
    print()
    sys.stdout.flush()
    prosit_tab = get_proteins(prosit_tab, msmstxt, proteinGroups, fasta)

    #Drop added columns
    print("Drop unnecessary columns and export new tab file")
    print()
    sys.stdout.flush()
    prosit_tab = prosit_tab.drop(columns=['Raw file','Scan number','Sequence','Charge','Scan event'], axis=1)
    
    #Export to new tab file
    prosit_tab.to_csv(args.output_tab, sep="\t", header=True, index=False)
    
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

def get_proteins(df, msmstxt, proteinGroups, fasta):
    df['Protein'] = df.apply(row_get_proteins, msmstxt=msmstxt, proteinGroups=proteinGroups, fasta=fasta, axis=1)
    return df

def row_get_proteins(row, msmstxt, proteinGroups, fasta):
    protein_str=''
    proteingroupids = msmstxt[(msmstxt['Raw file']==row['Raw file'])&(msmstxt['Scan number']==int(row['Scan number']))&(msmstxt['Sequence']==row['Sequence'])]['Protein group IDs'].str.split(';').tolist()[0]
    for proteingroupid in proteingroupids:
        proteins = proteinGroups[proteinGroups['id']==int(proteingroupid)]['Protein IDs'].str.split(';').tolist()[0] #This gives the names of the proteins like in the fasta
        for protein in proteins:
            #Check if the PSM sequence is in the full protein sequence (not all members of the protein group contain that PSM sequence)
            try:
                if row['Sequence'] in fasta[protein]:
                    if protein_str=='':
                        protein_str = protein
                    else:
                        protein_str = protein_str+";"+protein
            except:
                print(row['Raw file'])
                print(row['Scan number'])
                print(row['Sequence'])
                print(proteingroupid)
                print(protein)
                print(fasta[protein])
                print()
    #If no protein info is found, return this)
    if protein_str=='':
        protein_str = "NoProteinInfoFound"
    return protein_str

def load_fasta(in_fasta):

    #Init
    input_seqs = defaultdict()

    with open(in_fasta, 'r') as FR:
        lines = FR.readlines()
        accession = ''

        for i in range(0, len(lines)):
            lines[i] = lines[i].rstrip('\n')
            m = re.search('^>.+?\|(.+?)\|.+$', lines[i])
            if m:
                accession = m.group(1)
            else:
                seq = lines[i]
                input_seqs[accession] = seq
                accession = ''

    return input_seqs

def load_mq_files(msmstxt_path, proteingroupstxt_path):
    
    pd.set_option("display.max_colwidth", 10000)
    msmstxt = pd.read_csv(msmstxt_path, sep="\t")
    proteinGroups = pd.read_csv(proteingroupstxt_path, sep="\t")
    
    return(msmstxt, proteinGroups)

def reconstruct_psmid(df):
    df['SpecId'] = df.apply(row_reconstruct_psmid, axis=1)
    return df

def row_reconstruct_psmid(row):
    new_psmid = row['Raw file']+"-"+row['Scan number']+"-"+row["Sequence"]+"-"+row["Charge"]+"-"+row["Scan event"]
    return new_psmid

def methionine_convert(df, id_column="Sequence", id_column_2="Peptide", id_column_3="Protein"):
    df[id_column], df[id_column_2], df[id_column_3] = zip(*df.apply(row_methionine_convert, id_column=id_column, id_column_2=id_column_2, id_column_3=id_column_3, axis=1))
    return df

def row_methionine_convert(row, id_column="Sequence", id_column_2="Peptide", id_column_3="Protein"):
    new_seq = re.sub('m','M',row[id_column])
    new_seq = re.sub('M\(OX\)','M',new_seq)
    new_seq2 = re.sub('m','M',row[id_column_2])
    new_seq2 = re.sub('M\(OX\)','M',new_seq2)
    new_seq3 = re.sub('m','M',row[id_column_3])
    new_seq3 = re.sub('M\(OX\)','M',new_seq3)
    return(new_seq, new_seq2, new_seq3)

def split_prosit_psmid(df, id_column="SpecId"):
    df['Raw file'], df['Scan number'], df['Sequence'], df['Charge'], df['Scan event'] = zip(*df.apply(prosit_psmid_splitter, id_column=id_column, axis=1))
    return df

def prosit_psmid_splitter(row, id_column="SpecId"):    
    els = row[id_column].split("-")
    return(els[0]+"-"+els[1],els[2], els[3], els[4], els[5])  


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