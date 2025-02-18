import traceback
import sys
import os
import time
import argparse
import re
from collections import defaultdict
import itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

def main():
    starttime = time.time()

    print("\n#######################")
    print("# Combine with UniProt #")
    print("#######################\n")
    print("This program is part of the PROTEOFORMER pipeline\n")

    parser = argparse.ArgumentParser(description="This tool is part of the PROTEOFORMER pipeline. It combines "
                                                 "a (combined) fasta file of PROTEOFORMER with UniProt, ready for MS inspection.")

    parser.add_argument("--fasta", "-f", required=True, metavar="PATH", type=str, help="Input (combined) fasta file (mandatory)")
    parser.add_argument("--uniprot", "-u", required=True, metavar="PATH", type=str, help="Input UniProt fasta file (mandatory)")
    parser.add_argument("--workdir", "-w", metavar="FOLDER", default="", type=str, help="Working directory (default: CWD)")
    parser.add_argument("--output_fasta", "-o", metavar="PATH", default="comb_fasta_uniprot_proteoformer.fasta", type=str,
                        help="Output fasta file (default: comb_fasta_uniprot_proteoformer.fasta)")
    parser.add_argument("--overview_file", "-t", metavar="PATH", default="uniprot_proteoformer_overview.txt", type=str,
                        help="Overview file path (default: uniprot_proteoformer_overview.txt)")
    parser.add_argument("--venn_diagram", "-d", metavar="PATH", default="venn_diagram.png", type=str,
                        help="Path to the Venn diagram output file (default: venn_diagram.png)")

    args = parser.parse_args()

    args.workdir = args.workdir or os.getcwd()
    args.output_fasta = os.path.join(args.workdir, args.output_fasta)
    args.overview_file = os.path.join(args.workdir, args.overview_file)
    args.venn_diagram = os.path.join(args.workdir, args.venn_diagram)

    print("Parameters:")
    for arg, value in vars(args).items():
        print(f'    {arg:15}\t{value}')
    print()
    sys.stdout.flush()

    input_data = defaultdict(dict)
    sources = ["custom_pipeline", "uniprot"]
    inputs = [args.fasta, args.uniprot]

    data_groups = ["+".join(sorted(filter(None, comb))) for comb in itertools.combinations(sources + [""], 2)]

    for source, input_file in zip(sources, inputs):
        input_data = read_fasta(input_file, input_data, source)

    comb_data, counts = combine_data(input_data, data_groups)
    write_output(comb_data, args.output_fasta)
    write_overview(counts, args.overview_file)
    construct_venn(counts, args.venn_diagram)
    print(counts)

    print(f"\n-----------------------\n[{convert_time(time.time() - starttime)}]: PROGRAM COMPLETE\n-----------------------\n")
    sys.stdout.flush()


def construct_venn(counts, venn_file):
    if len(counts.keys()) == 3:
        subsets = (
            counts.get('custom_pipeline', 0),
            counts.get('uniprot', 0),
            counts.get('custom_pipeline+uniprot', 0)
        )
        set_labels = ('custom_pipeline', 'uniprot')
        textstr = "1: custom_pipeline\n2: uniprot"

        fig, ax = plt.subplots(nrows=1, ncols=1)
        v = venn2(subsets=subsets, set_labels=set_labels, ax=ax)

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(-0.03, 0.07, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
        plt.tight_layout()

        fig.savefig(venn_file)
        plt.close(fig)
def write_overview(counts, overview_file):
    with open(overview_file, 'w') as fw:
        fw.write(f"{'Files':30} {'Counts':20}\n")
        for group, count in counts.items():
            fw.write(f"{group:30} {count:20}\n")


def write_output(comb_data, fasta):
    with open(fasta, 'w') as fw:
        for acc, seq in comb_data.items():
            fw.write(f"{acc}\n{seq}\n")


def combine_data(input_data, data_groups):
    comb_data_per_file = defaultdict(dict)
    comb_data = {}
    counts = {group: 0 for group in data_groups}

    for source, entries in input_data.items():
        for acc, seq in entries.items():
            comb_data_per_file[seq][source] = acc

    for seq, sources in comb_data_per_file.items():
        classified_in = "+".join(sorted(sources.keys()))
        counts[classified_in] += 1
        main_acc = next(iter(sources.values()))
        side_acc = "#".join(list(sources.values())[1:])
        general_accession = f"{main_acc} [{side_acc}]" if side_acc else main_acc
        comb_data[general_accession] = seq

    return comb_data, counts


def read_fasta(fasta, input_data, source):
    with open(fasta, 'r') as fr:
        accession = ""
        sequence = ""
        for line in fr:
            line = line.strip()
            if line.startswith('>'):
                if accession:
                    input_data[source][accession] = sequence
                accession = line
                sequence = ""
            else:
                sequence += line
        if accession:
            input_data[source][accession] = sequence
    return input_data


def convert_time(seconds):
    h, m, s = map(int, (seconds // 3600, (seconds % 3600) // 60, seconds % 60))
    return f"{h}:{m:02}:{s:02}"


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        traceback.print_exc()
