import traceback
import sys
import os
import time
import argparse
import pandas as pd
import numpy as np
import re
from collections import defaultdict
import sqlite3 as sqlite
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import ConnectionPatch, ArrowStyle

__author__ = 'Steven Verbruggen'
#Execute 'python classify_proteoforms_percolator.py -h' for more information

def main():

    starttime = time.time()

    print
    print "###################################"
    print "# Classify proteoforms Percolator #"
    print "###################################"
    print
    print "This program is part of the PROTEOFORMER pipeline"
    print

    #Parse arguments from command line
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS,
                                     description="This tool is part of the PROTEOFORMER pipeline. It analyses "
                                                 "the different proteoforms found after MS analysis of a combined database "
                                                 "of PROTEOFORMER and UniProt. It searches for the classification of all "
                                                 "proteoforms found by PROTEOFORMER, not yet in UniProt and with confiramtion"
                                                 " of MS with Percolator (and eventual additional features from spectrum "
                                                 "predictors like Prosit and MS2PIP.",
                                     add_help=False)

    man_args = parser.add_argument_group("Mandatory parameters")
    man_args.add_argument("--proteoforms", "-p", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="Percolator proteoform csv file (mandatory)")
    man_args.add_argument("--mapping_file", "-M", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="Uniprot mapping file. You can get that from the Uniprot website, downloads, "
                                         "ID mapping (mandatory)")
    man_args.add_argument("--fasta_file", "-f", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="Combined fasta file of proteoformer and uniprot (mandatory)")
    man_args.add_argument("--ens_db", "-e", action="store", required=True, nargs="?", metavar="DB",
                          type=str, help="Ensembl database (mandatory)")

    opt_args = parser.add_argument_group("Optional parameters")
    opt_args.add_argument("--help", "-h", action="help", help="Show help message and exit")
    opt_args.add_argument("--workdir", "-w", action="store", required=False, nargs="?", metavar="FOLDER", default="",
                        type=str, help="Working directory (default: CWD)")
    opt_args.add_argument("--clustalo", "-c", action="store", required=False, nargs="?", metavar="PATH", default="clustalo",
                          type=str, help="Path to the Clustal Omega executable (default: clustalo)")
    opt_args.add_argument("--csv_output", "-x", action="store", required=False, nargs="?", metavar="PATH",
                          type=str, default="proteoform_classifications_percolator.csv", help="CSV output file of classifications (default: proteoform_classifications_percolator.csv)")
    opt_args.add_argument("--plot_file", "-o", action="store", required=False, nargs="?", metavar="PATH",
                          type=str, default="proteoform_class_percolator_plot.eps", help="Plot output file (default: proteoform_class_percolator_plot.eps)")

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

    #Load Percolator proteoforms
    proteoforms = load_proteoforms(args.proteoforms)

    #Make a Ensembl transcript id - Uniprot id mapping dict
    mapping_dict = construct_mapping_dict(args.mapping_file)

    #Get all sequences from the fasta file
    fasta_sequences = load_fasta(args.fasta_file)

    #Make tmp dir
    tmpdir = args.workdir + "/tmp"
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)

    #Do classification of proteoforms
    proteoforms = classify_proteoforms(proteoforms, mapping_dict, fasta_sequences, args.ens_db, args.workdir, tmpdir, args.clustalo)

    #Output identifications and counts as csv
    output_csv(proteoforms, args.csv_output)

    #Count classifications
    counts = count_classifications(proteoforms)
    print
    print
    print "Counts"
    for classification in sorted(counts.keys()):
        print classification+": "+str(counts[classification])
    print

    #Consturct distribution plot
    construct_plot(counts, args.plot_file)

    ## Protein group ENST00000356142_7_6387234_CDS_100db1;ENST00000356142_7_6374736_aTIS_111db1;ENST00000495499_7_6401880_ntr_100db1;ENST00000348035_7_6400138_CDS_100db1;ENST00000488373_7_6401964_ntr_100db1;ENST00000419413_X_137441779_ntr_100db1;ENST00000419413_X_137441836_ntr_001db3;ENST00000511164_4_46724252_ntr_100db1;ENST00000511164_4_46724351_ntr_110db1;ENST00000511164_4_46724408_ntr_001db3;ENST00000511164_4_46724108_ntr_100db1;ENST00000511164_4_46724156_ntr_100db1;ENST00000511164_4_46724171_ntr_100db1;ENST00000511164_4_46724024_ntr_100db1
    ## Canonical P63000 is probably false positive. Protein P63000 is for some reason not to be found in the 100search space MQ run.

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


# Construct distribution plot
def construct_plot(counts, output_file):
    # Present classification distribution with plot
    # Without confirming peptide: weglaten
    # Multiple variations: maybe print as different lane
    # NTR bundelen
    # Exon inclusion/exclusion bundelen als splicing isoforms
    # Found with ensembl: drop

    # Parse into the counts you want to show
    main_counts = defaultdict()
    main_counts['Splice variants'] = 0
    main_counts['Translation in non-coding region'] = 0
    splice_counts = defaultdict()
    splice_counts['Exon inclusion'] = 0
    splice_counts['Exon exclusion'] = 0
    ntr_counts = defaultdict()
    for classf in counts.keys():
        m_ntr = re.search('^Translation in NTR, (.+)$', classf)
        if re.search('without confirming peptide', classf):
            pass
        elif re.search('splice variant', classf):
            main_counts['Splice variants'] += counts[classf]
            splice_counts[classf] = counts[classf]
        elif re.search('^exon inclusion', classf):
            main_counts['Splice variants'] += counts[classf]
            splice_counts['Exon inclusion'] += counts[classf]
        elif re.search('^exon exclusion', classf):
            main_counts['Splice variants'] += counts[classf]
            splice_counts['Exon exclusion'] += counts[classf]
        elif classf == 'Exon substitution':
            main_counts['Splice variants'] += counts[classf]
            splice_counts[classf] = counts[classf]
        elif classf == 'Found with Ensembl check':
            pass
        elif m_ntr:
            main_counts['Translation in non-coding region'] += counts[classf]
            ntr_counts[m_ntr.group(1)] = counts[classf]
        else:
            main_counts[classf] = counts[classf]

    # Remove empty groups
    if main_counts['Splice variants'] == 0:
        del main_counts['Splice variants']
    if main_counts['Translation in non-coding region'] == 0:
        del main_counts['Translation in non-coding region']
    if splice_counts['Exon inclusion'] == 0:
        del splice_counts['Exon inclusion']
    if splice_counts['Exon exclusion'] == 0:
        del splice_counts['Exon exclusion']

    # Make plot
    fig = plt.figure(figsize=(30, 20))
    # set up subplot grid
    GridSpec(2, 3)

    # Big plot
    # Reorder dataframe
    df = get_ordered_dataframe(main_counts)
    labels_list = df['count'].values.tolist()
    # Plot
    ax1 = plt.subplot2grid((2, 3), (0, 0), colspan=2, rowspan=2)
    # Generate explode tuple
    explode_tuple = tuple()
    for val in df.loc[:, 'classification'].values.tolist():
        if val == 'Translation in non-coding region' or val == 'Splice variants':
            explode_tuple = explode_tuple + (0.2,)
        else:
            explode_tuple = explode_tuple + (0,)

    # Get colors
    colors = gen_color_list(df['classification'].values.tolist())
    print colors
    ax1 = plt.subplot2grid((2, 3), (0, 0), colspan=2, rowspan=2)
    # Generate explode tuple
    explode_tuple = tuple()
    for val in df.loc[:, 'classification'].values.tolist():
        if val == 'Translation in non-coding region' or val == 'Splice variants':
            explode_tuple = explode_tuple + (0.2,)
        else:
            explode_tuple = explode_tuple + (0,)
    patches, texts, autotexts = ax1.pie(df['count'], explode=explode_tuple, autopct='%.1f%%', colors=colors,
                                        textprops={'fontsize': 24}, pctdistance=0.85, labels=labels_list,
                                        labeldistance=1.05)
    # Label size
    for i in range(0, len(texts)):
        texts[i].set_fontsize(18)
    ax1.set_title("Overall proteoform categories", {'fontsize': 36})
    ax1.set_ylabel("")
    box = ax1.get_position()
    ax1.set_position([box.x0 - box.width * 0.1, box.y0,
                      box.width, box.height])
    lgd_labels = df['classification'].values.tolist()
    handles, labels = ax1.get_legend_handles_labels()  # Get legend information of last ax object
    leg = ax1.legend(handles, lgd_labels, fontsize=24, loc='upper center', bbox_to_anchor=(.5, 0), ncol=3)

    # Small plot 1
    df = pd.DataFrame.from_dict(splice_counts, orient="index")
    df.columns = ['count']
    df['classification'] = df.index
    labels_list = df['count'].values.tolist()
    ax2 = plt.subplot2grid((2, 3), (0, 2))
    colors = get_greens(df['classification'].values.tolist())
    df.plot(kind="pie", subplots="True", autopct='%.1f%%', ax=ax2, colors=colors, textprops={'fontsize': 24},
            fontsize=24, pctdistance=0.7, labels=labels_list, labeldistance=1.05)
    ax2.set_ylabel("")
    ax2.set_title("Splice variants", {'fontsize': 36})
    lgd_labels = df['classification'].values.tolist()
    handles, labels = ax2.get_legend_handles_labels()  # Get legend information of last ax object
    leg = ax2.legend(handles, lgd_labels, fontsize=24, loc='upper center', bbox_to_anchor=(.5, 0.075), ncol=2)

    # Small plot 2
    df = pd.DataFrame.from_dict(ntr_counts, orient="index")
    df.columns = ['count']
    df['classification'] = df.index
    labels_list = df['count'].values.tolist()
    ax3 = plt.subplot2grid((2, 3), (1, 2))
    colors = get_reds(df['classification'].values.tolist())
    df.plot(kind="pie", subplots="True", autopct='%.1f%%', ax=ax3, colors=colors, textprops={'fontsize': 24},
            fontsize=24, pctdistance=0.7, labels=labels_list, labeldistance=1.05)
    ax3.set_ylabel("")
    ax3.set_title("Non-coding region proteoforms", {'fontsize': 36})
    lgd_labels = rewrite_labels(df['classification'].values.tolist())
    handles, labels = ax3.get_legend_handles_labels()  # Get legend information of last ax object
    leg = ax3.legend(handles, lgd_labels, fontsize=24, loc='upper center', bbox_to_anchor=(.5, 0.075), ncol=2)

    # Connecting arrows
    (x, y) = texts[0].get_position()
    xy0 = (x + 0.08, y + 0.04)
    xyEnd = (-1.2, 0.5)
    con = ConnectionPatch(xyA=xy0, xyB=xyEnd, coordsA="data", coordsB="data",
                          axesA=ax1, axesB=ax2, arrowstyle=ArrowStyle.CurveFilledB(head_length=2, head_width=1),
                          color='black', linewidth=4)
    ax1.add_artist(con)

    (x, y) = texts[-1].get_position()
    xy0 = (x + 0.06, y)
    xyEnd = (-1.2, 0.7)
    con = ConnectionPatch(xyA=xy0, xyB=xyEnd, coordsA="data", coordsB="data",
                          axesA=ax1, axesB=ax3, arrowstyle=ArrowStyle.CurveFilledB(head_length=2, head_width=1),
                          color='black', linewidth=4)
    ax1.add_artist(con)

    # Output
    fig.suptitle('Classification of PROTEOFORMER-only identified proteoforms', fontsize=54)
    fig.savefig(output_file, dpi=1200)

    return


# Get ordered data frame
def get_ordered_dataframe(main_counts):
    df = pd.DataFrame.from_dict(main_counts, orient="index")
    df.columns = ['count']
    df['classification'] = df.index
    # Replace first and last row for non coding and splice variants
    b, c = df.loc['Translation in non-coding region'].copy(), df.iloc[-1].copy()
    df.loc['Translation in non-coding region'], df.iloc[-1] = c, b
    b, c = df.loc['Splice variants'].copy(), df.iloc[0].copy()
    df.loc['Splice variants'], df.iloc[0] = c, b
    # Sort values, except first and last row
    rest = df.iloc[1:-1].copy()
    rest.sort_values('classification', inplace=True)
    # Remake dataframe with sorted values
    first = pd.DataFrame([[df.iloc[0, 0], df.iloc[0, 1]]], columns=['count', 'classification'])
    last = pd.DataFrame([[df.iloc[-1, 0], df.iloc[-1, 1]]], columns=['count', 'classification'])
    new_df = pd.concat([first, rest, last], ignore_index=True)
    df = new_df
    df.reset_index(drop=True, inplace=True)

    return df


# Generate color list
def gen_color_list(class_list):
    all_labels = ['Splice variants', 'C-terminal extension', 'C-terminal truncation', 'Multiple variations',
                  'N-terminal extension', 'N-terminal truncation', 'Only amino acid substitutions', 'Out of frame ORF',
                  'dORF', 'uORF', 'Translation in non-coding region']

    #cmap = plt.cm.get_cmap('hsv', 10)
    cmap = plt.cm.tab10
    colors = cmap(np.linspace(0., 1., len(all_labels)))

    # Resturcture colors
    b = np.copy(colors[2, :])
    c = np.copy(colors[0, :])
    colors[2, :] = np.copy(c)
    colors[0, :] = np.copy(b)
    b = np.copy(colors[3, :])
    c = np.copy([0.8588, 0.8588, 0.5529, 1])
    colors[3, :] = np.copy(c)
    colors[-1, :] = np.copy(b)

    # Remove unneeded colors
    deleted_labels = 0
    for i in range(0, len(all_labels)):
        if all_labels[i] not in class_list:
            colors = np.delete(colors, i - deleted_labels, 0)
            deleted_labels += 1

    return colors


# Rewrite non coding labels
def rewrite_labels(labels):
    rewritten = list()
    translate = dict({'processed_transcript': 'Processed transcript',
                      'transcribed_processed_pseudogene': 'Transcribed processed pseudogene',
                      'processed_pseudogene': 'Processed pseudogene', 'retained_intron': 'Retained intron'})

    for label in labels:
        rewritten.append(translate[label])

    return rewritten


# Get red colors
def get_reds(class_list):
    all_labels = ['processed_transcript', 'transcribed_processed_pseudogene', 'processed_pseudogene', 'retained_intron']

    cmap = plt.cm.Reds
    colors = cmap(np.linspace(.1, .9, len(all_labels)))

    # Remove unneeded colors
    deleted_labels = 0
    for i in range(0, len(all_labels)):
        if all_labels[i] not in class_list:
            colors = np.delete(colors, i - deleted_labels, 0)
            deleted_labels += 1

    return colors


# Get green colors
def get_greens(class_list):
    all_labels = ['Exon inclusion', 'Exon exclusion', 'Exon substitution', 'C-terminal splice variant',
                  'N-terminal splice variant']

    cmap = plt.cm.Greens
    colors = cmap(np.linspace(.1, .9, len(all_labels)))

    # Remove unneeded colors
    deleted_labels = 0
    for i in range(0, len(all_labels)):
        if all_labels[i] not in class_list:
            colors = np.delete(colors, i - deleted_labels, 0)
            deleted_labels += 1

    return colors

#Count classifications
def count_classifications(proteoforms):

    #Init
    counts = defaultdict()

    for idx, row in proteoforms.iterrows():
        classification = row['classification']
        if classification in counts.keys():
            counts[classification]+=1
        else:
            counts[classification] = 1

    return counts


#Output identifications
def output_csv(proteoforms, csv_file):

    with open(csv_file, 'w') as FW:
        keys = ["id","base_proteoform","classification","proving_peptides","peptide_PEP","annotations","bin_codes","protein_group"]
        header_string= "id"+","+(','.join(keys))+"\n"
        FW.write(header_string)
        for idx, row in proteoforms.iterrows():
            value_string = str(row['id'])
            for key in keys[1:]:
                value_string = value_string+","+row[key]
            FW.write(value_string+"\n")

    return

#Classify proteoforms
def classify_proteoforms(proteoforms, mapping_dict, fasta_sequences, ens_db, workdir, tmpdir, clustalo_path):

    #Init
    proteoforms['classification'] = ""

    #Do the classification, protein group per protein group
    for idx, row in proteoforms.iterrows():
        #Determine base proteoform
        (row['base_proteoform'], row['base_annotation']) = find_base_proteoform(row)
        #Determine base canonical
        row['base_canonical'] = find_base_canonical(row, mapping_dict)
        #Get base proteoform sequence
        row['base_proteoform_seq'] = fasta_sequences[row['base_proteoform']]
        #Init
        row['classification'] = ''

        print row['protein_group']
        print row['proving_peptides']
        #print row['peptide_PEP']
        print row['base_proteoform']
        peptide_in_protein_printer(row['base_proteoform_seq'], row['proving_peptides'].split('|'))
        print row['base_canonical']

        #For non-coding: capture based on beforehand based on annotation
        if row['base_annotation']=="ntr":
            row['classification'] = "Translation in NTR"
            #Get more info about the non-coding transcript's biotype
            biotype = find_biotype(row['base_proteoform'], ens_db)
            row['classification'] = row['classification'] + ", " + biotype
            proteoforms.loc[idx, 'classification'] = row['classification']

        #For dORFs: capture also based on annotation
        elif row['base_annotation']=="3UTR":
            #Do additional check based on coordinates of suspective dORF and canonical ORF
            #Get transcript id of base proteoform
            transcript_id = get_transcript_id(row['base_proteoform'])
            #Get start of base proteoform
            base_proteoform_start = get_start(row['base_proteoform'])
            #Find start and stop of canonical ORF of that transcript
            (canonical_start, canonical_stop, strand) = get_canonical_start_stop(transcript_id, ens_db)
            if strand==1:
                if base_proteoform_start>canonical_stop:
                    row['classification'] = "dORF"
                    proteoforms.loc[idx, 'classification'] = row['classification']
            elif strand==-1:
                if base_proteoform_start<canonical_stop:
                    row['classification'] = "dORF"
                    proteoforms.loc[idx, 'classification'] = row['classification']

        #For extensions and truncations
        else:
            if row['base_canonical'] not in fasta_sequences:
                print "Base canonical not found: "+row['base_canonical']+" for protein group "+row['protein_group']
            row['base_canonical_seq'] = fasta_sequences[row['base_canonical']]
            peptide_in_protein_printer(row['base_canonical_seq'], row['proving_peptides'].split('|'))
            #Align with ClustalO
            #Prepare fasta
            tmp_fa = prepare_fasta(tmpdir, row['base_proteoform'], row['base_proteoform_seq'], row['base_canonical'], row['base_canonical_seq'])
            #Align
            output_fa = clustal_align(workdir, tmpdir, row['base_proteoform'], tmp_fa, clustalo_path)
            # Read aligned sequences
            (aligned_base_proteoform, aligned_base_canonical) = read_clustal_output(output_fa, row['base_proteoform'], row['base_canonical'])
            #Remove fasta
            os.system("rm -rf "+tmp_fa)
            #Calculate mapped percentage
            perc_mapped = calc_perc_mapped(aligned_base_proteoform, aligned_base_canonical)
            print("Perc mapped: "+str(perc_mapped))
            # Get peptides and their position on the full protein sequence
            proteoform_peptides = get_peptides(row['base_proteoform_seq'], row['proving_peptides'].split('|'))
            print_dict(proteoform_peptides)
            print(aligned_base_proteoform)
            print(aligned_base_canonical)
            # Classify proteoform
            row['classification'] = classify_proteoform(proteoform_peptides,
                                                           aligned_base_proteoform,
                                                           aligned_base_canonical, perc_mapped,
                                                           row['base_proteoform'], row['base_canonical'],
                                                           row['base_proteoform_seq'], row['base_canonical_seq'], ens_db)
            proteoforms.loc[idx, 'classification'] = row['classification']

        if row['classification']!='':
            print row['classification']
        else:
            print
            print
            print "-----UNCLASSIFIED!!-------"
            print "Base proteoform: " + row['base_proteoform']
            print "Base canonical: " + row['base_canonical']
            print "Protein group: " + row['protein_group']
            print "--------------------------"
            print
            print
        print

    return proteoforms

#Classify proteoform
def classify_proteoform(proteoform_peptides, aligned_base_proteoform, aligned_base_canonical, perc_mapped, base_proteoform, base_canonical, base_proteoform_seq, base_canonical_seq, ens_db):

    #Init
    classification=""

    #Check if enough identical positions
    if perc_mapped>0.35:
        #Correct start and stop positions based on alignment
        for peptide in proteoform_peptides:
            proteoform_peptides[peptide]['corr_start'] = correct_pos(aligned_base_proteoform, proteoform_peptides[peptide]['start'])
            print("Corrected start: "+str(proteoform_peptides[peptide]['corr_start']))
            proteoform_peptides[peptide]['corr_end'] = correct_pos(aligned_base_proteoform, proteoform_peptides[peptide]['end'])
            print("Corrected end: "+str(proteoform_peptides[peptide]['corr_end']))
            # If N-terminal extension: proteoform starts with AA's and canonical with '-'s. There is also a corr peptide start in this stretch
            m1 = re.search('^[A-Za-z]+', aligned_base_proteoform)
            m2 = re.search('^(\-+)', aligned_base_canonical)
            # If N-terminal truncation: proteoform starts with '-'s and canonical with AA's. Peptide at start or second position of truncation
            m3 = re.search('^(\-+)', aligned_base_proteoform)
            m4 = re.search('^[A-Za-z]+', aligned_base_canonical)
            # If C-terminal extension: proteoform ends with AA's and canonical with '-'s. There is a corr_peptide end in this stretch
            m5 = re.search('[A-Za-z]+$', aligned_base_proteoform)
            m6 = re.search('(\-+)$', aligned_base_canonical)
            # If C-terminal truncation: proteoform ends with '-'s and canonical with AA's. Peptide found at last or second last position of proteoform
            m7 = re.search('(\-+)$', aligned_base_proteoform)
            m8 = re.search('[A-Za-z]+$', aligned_base_canonical)
            # If exon exclusion: gap in proteoform
            m9 = re.search('^([A-Za-z]+)\-+([A-Za-z]+)$', aligned_base_proteoform)
            # multiple
            m9b = re.search('^[A-Za-z]+\-+[A-Za-z]+\-+[A-Za-z]+$', aligned_base_proteoform)
            m9c = re.search('^[A-Za-z]+\-+[A-Za-z]+\-+[A-Za-z]+\-+[A-Za-z]+$', aligned_base_proteoform)
            # If exon inclusion: gap in canonical
            m10 = re.search('^([A-Za-z]+)\-+([A-Za-z]+)$', aligned_base_canonical)
            # multiple
            m10b = re.search('^[A-Za-z]+\-+[A-Za-z]+\-+[A-Za-z]+$', aligned_base_canonical)
            m10c = re.search('^[A-Za-z]+\-+[A-Za-z]+\-+[A-Za-z]+\-+[A-Za-z]+$', aligned_base_canonical)
            # If there are no indels from start to stop in both sequence, then only SAVs
            m11 = re.search('^[A-Za-z]+$', aligned_base_proteoform)
            m12 = re.search('^[A-Za-z]+$', aligned_base_canonical)
            # Sequence absent in Uniprot
            if (aligned_base_proteoform == aligned_base_canonical):
                classification = "Sequence absent in UniProt"

            # N-terminal extension
            if (m1 and m2):
                # Get the rest of the sequences (-initiator methionine) and calculate similarity
                rest_of_proteoform = aligned_base_proteoform[(len(m2.group(1)) + 1):]
                rest_of_canonical = aligned_base_canonical[(len(m2.group(1)) + 1):]
                # Check if the rest of the sequences is identical
                if rest_of_proteoform == rest_of_canonical:
                    peptide_found = False
                    for peptide in proteoform_peptides:
                        if (proteoform_peptides[peptide]['corr_start'] <= len(m2.group(1))):
                            #check if peptide present in novel but not in canonical
                            if (peptide in base_proteoform_seq) and (peptide not in base_canonical_seq):
                                peptide_found = True
                    if peptide_found == True:
                        classification = 'N-terminal extension'
                    else:
                        classification = 'N-terminal extension, without confirming peptide'
                else:
                    if check_only_sav(rest_of_proteoform, rest_of_canonical):
                        classification = "Multiple variations"
                    else:
                        # Search longest identical stretch of characters
                        longest_stretch, perc = search_longest_stretch(rest_of_proteoform, rest_of_canonical)
                        # Check if stretch can be found at the end
                        if rest_of_proteoform.endswith(longest_stretch) and perc > 0.1:
                            classification = "N-terminal splice variant"
                        else:
                            classification = "Multiple variations"
            # C-terminal extension
            elif (m5 and m6):
                # Get the rest of the sequences and calculate similarity
                rest_of_proteoform = aligned_base_proteoform[:(len(aligned_base_proteoform) - len(m6.group(1)))]
                rest_of_canonical = aligned_base_canonical[:(len(aligned_base_canonical) - len(m6.group(1)))]
                # Check if the rest of the sequences is identical
                if rest_of_proteoform == rest_of_canonical:
                    peptide_found = False
                    for peptide in proteoform_peptides:
                        if (proteoform_peptides[peptide]['corr_end'] >= (len(aligned_base_proteoform) - len(m6.group(1)))):
                            #Check if peptide present in novel but not in canonical
                            if (peptide in base_proteoform_seq) and (peptide not in base_canonical_seq):
                                peptide_found=True
                    if peptide_found:
                        classification = "C-terminal extension"
                    else:
                        classification = 'C-terminal extension, without confirming peptide'
                else:
                    if check_only_sav(rest_of_proteoform, rest_of_canonical):
                        classification = "Multiple variations"
                    else:
                        # Search longest identical stretch of characters
                        longest_stretch, perc = search_longest_stretch(rest_of_proteoform, rest_of_canonical)
                        # Check if stretch can be found at the start
                        if rest_of_proteoform.startswith(longest_stretch) and perc > 0.1:
                            classification = "C-terminal splice variant"
                        else:
                            classification = "Multiple variations"
            # N-terminal truncation
            elif (m3 and m4):
                # Get the rest of the sequences (-initiator methionine) and calculate similarity
                rest_of_proteoform = aligned_base_proteoform[(len(m3.group(1)) + 1):]
                rest_of_canonical = aligned_base_canonical[(len(m3.group(1)) + 1):]
                # Check if the rest of the sequences is identical
                if rest_of_proteoform == rest_of_canonical:
                    # Check for peptide
                    peptide_found = False
                    for peptide in proteoform_peptides:
                        if (proteoform_peptides[peptide]['corr_start'] in [len(m3.group(1)) + 1, len(m3.group(1)) + 2]):
                            if (peptide in base_proteoform_seq) and (peptide not in base_canonical_seq):
                                peptide_found=True
                            else:
                                regex = r"^\w*[KR]" + re.escape(peptide) + r"\w*$" #Take in account that there needs to be tryptic cleavage site that can create the peptide
                                m_novel = re.search(regex, base_proteoform_seq, re.IGNORECASE) #Mostly the novel one has a near-cognate methionine that enables the detectability of the new peptide
                                m_can = re.search(regex, base_canonical_seq, re.IGNORECASE) #Mostly the canonical does not have the methionine
                                if m_novel and (not m_can):
                                    peptide_found=True
                                else:
                                    regex = r"^M"+ re.escape(peptide) + r"\w*$"  # Take in account that there needs to be a methionine that can be split off biologically and thus creates the peptide
                                    m_novel = re.search(regex, base_proteoform_seq,re.IGNORECASE)  # Mostly the novel one has a near-cognate methionine that enables the detectability of the new peptide
                                    m_can = re.search(regex, base_canonical_seq,re.IGNORECASE)  # Mostly the canonical does not have the methionine
                                    if m_novel and (not m_can):
                                        peptide_found = True

                    if peptide_found == True:
                        classification = 'N-terminal truncation'
                    else:
                        classification = 'N-terminal truncation, without confirming peptide'
                else:
                    if check_only_sav(rest_of_proteoform, rest_of_canonical):
                        classification = "Multiple variations"
                    else:
                        # Search longest identical stretch of characters
                        longest_stretch, perc = search_longest_stretch(rest_of_proteoform, rest_of_canonical)
                        # Check if stretch can be found at the end
                        if rest_of_proteoform.endswith(longest_stretch) and perc > 0.1:
                            classification = "N-terminal splice variant"
                        else:
                            classification = "Multiple variations"
            # C-terminal truncation
            elif (m7 and m8):
                # Get the rest of the sequences and calculate similarity
                rest_of_proteoform = aligned_base_proteoform[:(len(aligned_base_proteoform) - len(m7.group(1)))]
                rest_of_canonical = aligned_base_canonical[:(len(aligned_base_proteoform) - len(m7.group(1)))]
                # Check if the rest of the sequences is identical
                if rest_of_proteoform == rest_of_canonical:
                    # Check for peptide
                    peptide_found = False
                    for peptide in proteoform_peptides:
                        if (proteoform_peptides[peptide]['corr_end'] in [len(aligned_base_proteoform) - len(m7.group(1)),len(aligned_base_proteoform) - len(m7.group(1)) - 1]):
                            if (peptide in base_proteoform_seq) and (peptide not in base_canonical_seq):
                                peptide_found=True
                    if peptide_found:
                        classification = "C-terminal truncation"
                    else:
                        classification = 'C-terminal truncation, without confirming peptide'
                else:
                    if check_only_sav(rest_of_proteoform, rest_of_canonical):
                        classification = "Multiple variations"
                    else:
                        # Search longest identical stretch of characters
                        longest_stretch, perc = search_longest_stretch(rest_of_proteoform, rest_of_canonical)
                        # Check if stretch can be found at the start
                        if rest_of_proteoform.startswith(longest_stretch) and perc > 0.1:
                            classification = "C-terminal splice variant"
                        else:
                            classification = "Multiple variations"
            # Exon exclusion
            elif (m9):
                regex = r"^" + re.escape(m9.group(1)) + r"[A-Za-z]+?" + re.escape(m9.group(2)) + r"$"
                m_excl = re.search(regex, aligned_base_canonical)
                if m_excl:
                    classification = "exon exclusion"
                else:
                    classification = "exon exclusion, more complex"
            # Exon inclusion
            elif (m10):
                regex = r"^" + re.escape(m10.group(1)) + r"[A-Za-z]+?" + re.escape(m10.group(2)) + r"$"
                m_incl = re.search(regex, aligned_base_proteoform)
                if m_incl:
                    classification = "exon inclusion"
                else:
                    classification = "exon inclusion, more complex"
            elif (m9b):
                classification = "exon exclusion"
            elif (m9c):
                classification = "exon exclusion"
            elif (m10b):
                classification = "exon inclusion"
            elif (m10c):
                classification = "exon inclusion"
            # Only amino acid variations
            elif (m11 and m12):
                if aligned_base_canonical == aligned_base_proteoform:
                    if base_canonical!="non-coding":
                        classification = "Found with Ensembl check"
                else:
                    # Check only sav
                    if check_only_sav(aligned_base_proteoform, aligned_base_canonical):
                        classification = "Only amino acid substitutions"
                    else:
                        classification = "Multiple variations"
    else:
        #No homolog sequence
        # Get transcript id of base proteoform
        transcript_id = get_transcript_id(base_proteoform)
        # Get start of base proteoform
        base_proteoform_start = get_start(base_proteoform)
        # Find start and stop of canonical ORF of that transcript
        (canonical_start, canonical_stop, strand) = get_canonical_start_stop(transcript_id, ens_db)
        #Check for uORF based on start positions
        if strand == 1:
            if base_proteoform_start < canonical_start:
                classification = "uORF"
        elif strand == -1:
            if base_proteoform_start > canonical_start:
                classification = "uORF"

        #Search further: out of frame ORF?
        if classification=='':
            #Get phasing info of the canonical product
            exon_info = get_canonical_exon_info(transcript_id, canonical_start, canonical_stop, strand, ens_db)
            #Check if proteoform start is out of phase with canonical product
            in_phase = check_if_in_phase(base_proteoform_start, exon_info, strand)
            if in_phase=='N':
                classification='Out of frame ORF'

    return classification

#Check if proteoform start is in phase
def check_if_in_phase(base_proteoform_start, exon_info, strand):

    #Init
    in_phase=""

    if strand==1:
        for rank in sorted(exon_info.keys()):
            #Search the exon in which the start is in
            if base_proteoform_start>=exon_info[rank]['start'] and base_proteoform_start<=exon_info[rank]['stop']:
                #Check if in same phase
                phase_diff = (base_proteoform_start - exon_info[rank]['start']) % 3
                if phase_diff==0:
                    in_phase='Y'
                elif phase_diff==1 or phase_diff==2:
                    in_phase='N'
    elif strand==-1:
        for rank in sorted(exon_info.keys()):
            if base_proteoform_start<=exon_info[rank]['start'] and base_proteoform_start>=exon_info[rank]['stop']:
                phase_diff = (exon_info[rank]['start'] - base_proteoform_start) % 3
                if phase_diff==0:
                    in_phase='Y'
                elif phase_diff==1 or phase_diff==2:
                    in_phase='N'

    return in_phase

#Get exon info of canonical ORF
def get_canonical_exon_info(transcript_id, canonical_start, canonical_stop, strand, ens_db):

    #Init
    exon_info = defaultdict(lambda: defaultdict())
    transcript_exons = defaultdict(lambda: defaultdict())

    #Connect to ensembl
    try:
        conn = sqlite.connect(ens_db)
    except:
        print "Could not connect to "+ens_db
        sys.exit()

    with conn:
        cur = conn.cursor()

        #Get all transcript exons
        query = "SELECT et.rank, e.seq_region_start, e.seq_region_end, e.phase, e.end_phase FROM exon AS e " \
                "JOIN exon_transcript AS et ON e.exon_id = et.exon_id JOIN transcript AS tr ON "\
                "et.transcript_id = tr.transcript_id WHERE tr.stable_id = '"+transcript_id+"' ORDER BY et.rank;"
        cur.execute(query)
        output = cur.fetchall()

        #Parse into defaultdict
        for i in range(0,len(output)):
            transcript_exons[output[i][0]]['seq_region_start'] = output[i][1]
            transcript_exons[output[i][0]]['seq_region_end'] = output[i][2]
            transcript_exons[output[i][0]]['phase'] = output[i][3]
            transcript_exons[output[i][0]]['end_phase'] = output[i][4]

        #Only keep the translated exon parts. Go from absolute ensembl exon stradness to relative
        if strand==1:
            for rank in sorted(transcript_exons.keys()):
                #Both in this exon: trim both sides
                if canonical_start>=transcript_exons[rank]['seq_region_start'] and canonical_stop<=transcript_exons[rank]['seq_region_end']:
                    exon_info[rank]['start'] = canonical_start
                    exon_info[rank]['phase'] = 0
                    exon_info[rank]['stop'] = canonical_stop
                    exon_info[rank]['stop_phase'] = 2
                #Canonical start in exon: trim
                elif canonical_start>=transcript_exons[rank]['seq_region_start'] and canonical_start<=transcript_exons[rank]['seq_region_end']:
                    exon_info[rank]['start'] = canonical_start
                    exon_info[rank]['phase'] = 0
                    exon_info[rank]['stop'] = transcript_exons[rank]['seq_region_end']
                    exon_info[rank]['stop_phase'] = transcript_exons[rank]['end_phase']
                #Canonical stop in exon: trim
                elif canonical_stop>=transcript_exons[rank]['seq_region_start'] and canonical_stop<=transcript_exons[rank]['seq_region_end']:
                    exon_info[rank]['start'] = transcript_exons[rank]['seq_region_start']
                    exon_info[rank]['phase'] = transcript_exons[rank]['phase']
                    exon_info[rank]['stop'] = canonical_stop
                    exon_info[rank]['stop_phase'] = 2
                #No canonical start or stop in exon, but inside canonical translation zone: just copy
                elif canonical_start<transcript_exons[rank]['seq_region_start'] and canonical_stop>transcript_exons[rank]['seq_region_end']:
                    exon_info[rank]['start'] = transcript_exons[rank]['seq_region_start']
                    exon_info[rank]['phase'] = transcript_exons[rank]['phase']
                    exon_info[rank]['stop'] = transcript_exons[rank]['seq_region_end']
                    exon_info[rank]['stop_phase'] = transcript_exons[rank]['end_phase']
        elif strand==-1:
            for rank in sorted(transcript_exons.keys()):
                #Both in this exon: trim both sides
                if canonical_start<=transcript_exons[rank]['seq_region_end'] and canonical_stop>=transcript_exons[rank]['seq_region_start']:
                    exon_info[rank]['start'] = canonical_start
                    exon_info[rank]['phase'] = 0
                    exon_info[rank]['stop'] = canonical_stop
                    exon_info[rank]['stop_phase'] = 2
                #Canonical start in exon: trim
                elif canonical_start<=transcript_exons[rank]['seq_region_end'] and canonical_start>=transcript_exons[rank]['seq_region_start']:
                    exon_info[rank]['start'] = canonical_start
                    exon_info[rank]['phase'] = 0
                    exon_info[rank]['stop'] = transcript_exons[rank]['seq_region_start']
                    exon_info[rank]['stop_phase'] = transcript_exons[rank]['phase']
                #Canonical stop in exon: trim
                elif canonical_stop<=transcript_exons[rank]['seq_region_end'] and canonical_stop>=transcript_exons[rank]['seq_region_start']:
                    exon_info[rank]['start'] = transcript_exons[rank]['seq_region_end']
                    exon_info[rank]['phase'] = transcript_exons[rank]['end_phase']
                    exon_info[rank]['stop'] = canonical_stop
                    exon_info[rank]['stop_phase'] = 2
                #No canonical start or stop in exon but in canonical translation zone: just copy
                elif canonical_start>transcript_exons[rank]['seq_region_end'] and canonical_stop<transcript_exons[rank]['seq_region_start']:
                    exon_info[rank]['start'] = transcript_exons[rank]['seq_region_end']
                    exon_info[rank]['phase'] = transcript_exons[rank]['end_phase']
                    exon_info[rank]['stop'] = transcript_exons[rank]['seq_region_start']
                    exon_info[rank]['stop_phase'] = transcript_exons[rank]['phase']

    return exon_info

#Search longest identical stretch
def search_longest_stretch(proteoform_seq, canonical_seq):

    #Init
    longest_stretch = ""
    cur_stretch = ""


    for pos in range(0, len(proteoform_seq)):
        #Check if identical so that stretch should be expanded
        if proteoform_seq[pos]==canonical_seq[pos]:
            cur_stretch = cur_stretch + proteoform_seq[pos]
        else:
            cur_stretch = ""
        #Check if it is the longest stretch so far
        if len(cur_stretch)>len(longest_stretch):
            longest_stretch = cur_stretch

    #Calculate percentage covered by longest stretch
    percentage = float(len(longest_stretch))/float(len(canonical_seq))

    return longest_stretch, percentage

#Check only SAV
def check_only_sav(aligned_proteoform_seq, aligned_canonical_seq):

    #Init
    only_sav=True
    prev_not_matched = False

    for pos in range(0, len(aligned_proteoform_seq)):
        if prev_not_matched==False:
            if aligned_proteoform_seq[pos]!=aligned_canonical_seq[pos]:
                prev_not_matched=True
        else:
            if aligned_proteoform_seq[pos]==aligned_canonical_seq[pos]:
                prev_not_matched=False
            else:
                only_sav=False
                break

    return only_sav

#Correct position based on alignment
def correct_pos(aligned_seq, org_pos):

    #subtract each time the amount of position which are not in an insert until i=0
    i = int(org_pos)
    rest_seq = aligned_seq
    tmp_seq = ""
    corrected = 0
    while i!=0:
        tmp_seq = rest_seq[:i]
        rest_seq = rest_seq[i:]
        i = tmp_seq.count('-')
        corrected = corrected + len(tmp_seq)

    return corrected

#Get peptides and their positions on the full protein sequence
def get_peptides(protein_seq, peptides):

    #Init
    proteoform_peptides = defaultdict(lambda: defaultdict())

    #For all peptides
    for peptide in peptides:
        #Get start, stop position and sequence
        proteoform_peptides[peptide]['sequence'] = peptide
        proteoform_peptides[peptide]['start'] = protein_seq.find(peptide) #0-based coordinates
        proteoform_peptides[peptide]['end'] = proteoform_peptides[peptide]['start'] + len(peptide) - 1 #To correct for end counting

    return proteoform_peptides

#Calculate mapped percentage
def calc_perc_mapped(aligned_base_proteoform, aligned_base_canonical):

    #Init
    perc_mapped = 0
    identical_positions = 0
    length_without_indels = 0

    #Go over alignment per base
    for base_pos in range(0, len(aligned_base_proteoform)):
        #Check for indels
        if aligned_base_proteoform[base_pos]!='-' and aligned_base_canonical[base_pos]!='-':
            length_without_indels+=1
            #Check if bases are identical
            if aligned_base_proteoform[base_pos]==aligned_base_canonical[base_pos]:
                #If identical, +1 identical positions
                identical_positions+=1

    #Calculate percentage identical positions
    perc_mapped = float(identical_positions)/float(length_without_indels)

    return perc_mapped

#Read clustal output
def read_clustal_output(fasta, base_proteoform, base_canonical):

    #Init
    clustal_seqs = defaultdict()

    #Read output
    with open(fasta, 'r') as FR:
        lines = FR.readlines()
        lines = map(lambda x: x.rstrip("\n"), lines)
        lines = map(lambda x: x.rstrip("\r"), lines)
        # Init
        accession = ""
        sequence = ""
        for i in range(0, len(lines)):
            # Check if accession -> new entry
            if ((re.search('^>', lines[i]))):
                if (i != 0):
                    # Save previous entry
                    clustal_seqs[accession[1:]] = sequence
                # Get new entry accession and empty sequence
                accession = lines[i]
                sequence = ""
            else:
                # Concat next line of sequence
                sequence += lines[i]
        # Save the last entry
        if (sequence != ""):
            clustal_seqs[accession[1:]] = sequence
            accession = ""
            sequence = ""

    aligned_base_proteoform = clustal_seqs[base_proteoform]
    aligned_base_canonical = clustal_seqs[base_canonical]

    return (aligned_base_proteoform, aligned_base_canonical)

#Align with ClustalO
def clustal_align(workdir, tmpdir, base_proteoform, input_fa, clustalo_path):

    #Init output
    output_fa = tmpdir+"/"+base_proteoform+"aligned_output.fa"

    #Execute
    if clustalo_path=='clustalo':
        command = clustalo_path+" -i "+input_fa+" -o "+output_fa+" --outfmt fa --force --wrap=1000000000" #Force to overwrite, output in fasta format
        os.system(command)
    else:
        m = re.search('^(.+)clustalo$', clustalo_path)
        if m:
            os.chdir(m.group(1))
            command = "clustalo -i " + input_fa + " -o " + output_fa + " --outfmt fa --force --wrap=1000000000"  # Force to overwrite, output in fasta format
            os.system(command)
            os.chdir(workdir)

    return output_fa

#Prepare tmp fasta
def prepare_fasta(tmpdir, base_proteoform, base_proteoform_seq, base_canonical, base_canonical_seq):

    #Define path of tmp fasta file
    tmp_fasta = tmpdir+"/tmp_"+base_proteoform+".fa"

    with open(tmp_fasta, 'w') as FW:
        FW.write(">"+base_proteoform+"\n")
        FW.write(base_proteoform_seq+"\n")
        FW.write(">"+base_canonical+"\n")
        FW.write(base_canonical_seq+"\n")

    return tmp_fasta

#Get canonical start and stop of a transcript
def get_canonical_start_stop(transcript_id, ens_db):

    #Init
    can_start = 0
    can_stop  = 0

    #Connect to ensembl
    conn = sqlite.connect(ens_db)

    with conn:
        cur = conn.cursor()

        #Get canonical translation id
        query_transcript = "SELECT canonical_translation_id, seq_region_strand FROM transcript WHERE stable_id='"+transcript_id+"';"
        cur.execute(query_transcript)
        (canonical_translation_id, strand) = cur.fetchone()

        #Get start and stop info of canonical translation
        query_translation = "SELECT seq_start, seq_end, start_exon_id, end_exon_id FROM translation WHERE translation_id='"+str(canonical_translation_id)+"';"
        cur.execute(query_translation)
        (seq_start, seq_end, start_exon_id, end_exon_id) = cur.fetchone()

        #Get starts and ends of start and end exons
        query_start_exon = "SELECT seq_region_start, seq_region_end FROM exon where exon_id='"+str(start_exon_id)+"';"
        cur.execute(query_start_exon)
        (start_exon_start, start_exon_end) = cur.fetchone()
        query_end_exon = "SELECT seq_region_start, seq_region_end FROM exon where exon_id='"+str(end_exon_id)+"';"
        cur.execute(query_end_exon)
        (end_exon_start, end_exon_end) = cur.fetchone()

        #Calculate canonical start and end
        if strand==1:
            can_start = start_exon_start + seq_start - 1
            can_stop  = end_exon_start + seq_end - 1
        elif strand==-1:
            can_start = start_exon_end - seq_start + 1
            can_stop = end_exon_end - seq_end + 1

    return (can_start, can_stop, strand)

#Get start of proteoform
def get_start(proteoform):

    #Init
    start = ""

    m = re.search('^ENST\d+?\_.+?\_(\d+?)\_', proteoform)
    if m:
        start = int(m.group(1))

    return start


#Find biotype of proteoform
def find_biotype(proteoform, ens_db):

    #Get transcript_id
    transcript_id = get_transcript_id(proteoform)

    #Connect to ensembl
    conn = sqlite.connect(ens_db)
    with conn:
        cur = conn.cursor()
        query = "select biotype from transcript where stable_id='"+transcript_id+"';"
        cur.execute(query)
        biotype = cur.fetchone()[0]

    return biotype

#Get transcript id out of proteoform accession
def get_transcript_id(proteoform):

    #Init
    transcript_id=""

    m = re.search('^(ENST\d+?)\_', proteoform)
    if m:
        transcript_id = m.group(1)

    return transcript_id

#Find the base canonical protein
def find_base_canonical(row, mapping_dict):

    #Init
    base_canonical = ""

    #Non-coding transcripts have no canonical
    if row['base_annotation']=='ntr':
        base_canonical = "non-coding"
    else:
        #Get transcript ID of base proteoform
        transcript_id=""
        m = re.search('^(ENST.+?)\_', row['base_proteoform'])
        if m:
            transcript_id = m.group(1)
        #Search in mapping dict
        if transcript_id in mapping_dict.keys():
            base_canonical = mapping_dict[transcript_id]

    return base_canonical

#Find base proteoform of each proteoform
def find_base_proteoform(row):

    #In Percolator, every protein of the protein group, contains the peptide which makes the protein a separate protein group
    #Therefore, all proteins are as important. In order to select the ideal base proteoform, certain annotations are given
    #priority. Further, the first occuring protein is taken as base proteoform.

    #Order of priority of annotations
    sorted_annotations = {'aTIS' : 0,
                          'CDS' : 1,
                          '5UTR' : 2,
                          '3UTR' : 3,
                          'ntr' : 4}
    proteins = row['protein_group'].split(';')

    #Give the priority order numbers to the proteins and select the first occuring protein with the lowest priority number.
    annotations = []
    for protein in proteins:
        m = re.search('^ENST\d+\_.+\_\d+\_(.+)\_.+db.+$', protein)
        annotations.append(sorted_annotations[m.group(1)])
    prioretized_annotation = min(annotations)
    idx_prioretized_annotation = annotations.index(prioretized_annotation)

    base_proteoform = proteins[idx_prioretized_annotation]

    #Get prioretized annotation
    base_annotation = ""
    for key in sorted_annotations.keys():
        if sorted_annotations[key]==prioretized_annotation:
            base_annotation = key

    return base_proteoform, base_annotation

#Load fasta file
def load_fasta(fasta):

    #Init
    fasta_sequences = defaultdict()

    with open(fasta, 'r') as FR:
        lines = FR.readlines()
        lines = map(lambda x: x.rstrip("\n"), lines)
        lines = map(lambda x: x.rstrip("\r"), lines)
        # Init
        accession = ""
        sequence = ""
        for i in range(0, len(lines)):
            # Check if accession -> new entry
            if ((re.search('^>', lines[i]))):
                if (i != 0):
                    # Save previous entry
                    id = get_id_from_accession(accession)
                    fasta_sequences[id] = sequence
                # Get new entry accession and empty sequence
                accession = lines[i]
                sequence = ""
            else:
                # Concat next line of sequence
                sequence += lines[i]
        # Save the last entry
        if (sequence != ""):
            id = get_id_from_accession(accession)
            fasta_sequences[id] = sequence
            accession = ""
            sequence = ""

    return fasta_sequences

#Get proteoform id out of accession
def get_id_from_accession(accession):

    #Init
    id=""

    m = re.search('^.+?\|(.+?)\|.+', accession)
    if m:
        id = m.group(1)

    return id

#Construct ID mapping dict
def construct_mapping_dict(mapping_file):

    #Init
    #mapping_dict[ensembl transcript id] = uniprot accession
    mapping_dict = defaultdict()

    with open(mapping_file, 'r') as FR:
        for line in FR:
            values = line.split('\t')
            transcript_ids = values[19].split('; ')
            for transcript_id in transcript_ids:
                mapping_dict[transcript_id] = values[0]


    return mapping_dict

#Load proteoforms from Percolator parsed csv list
def load_proteoforms(in_file):

    proteoforms = pd.read_csv(in_file, sep=",", header=0)

    return proteoforms

## Peptide in protein sequence visualizer
def peptide_in_protein_printer(protein_seq, pep_seqs):
    pep_seqs = list(pep_seqs)
    protein_seq = protein_seq.lower()
    for pep_seq in pep_seqs:
        pep_seq = pep_seq.lower()
        regex = r"^(\w*)" + re.escape(pep_seq) + r"(\w*)$"
        m = re.search(regex, protein_seq, re.IGNORECASE)
        if m:
            protein_seq = m.group(1)+pep_seq.upper()+m.group(2)
    if protein_seq.islower():
        print("Not present: "+protein_seq)
    else:
        print(protein_seq)
    return

## Data Dumper for recursively printing nested dictionaries and defaultDicts ##
def print_dict(dictionary, indent='', braces=0):
    """
    :param dictionary: dict or defaultdict to be printed
    :param indent: indentation
    :param braces:
    :return: void
    """

    for key, value in dictionary.iteritems():
        if isinstance(value, dict):
            print '%s%s%s%s' % (indent, braces * '[', key, braces * ']')
            print_dict(value, indent + '  ', braces)
        else:
            print indent + '%s = %s' % (key, value)

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