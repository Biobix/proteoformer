import traceback
import sys
import os
import time
import argparse
import re
from collections import defaultdict
import sqlite3 as sqlite
#from dict_functions import dict_funcs
import pandas as pd
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import ConnectionPatch, ArrowStyle

__author__ = 'Steven Verbruggen'
#Execute 'python analyse_proteoforms.py -h' for more information

def main():

    starttime = time.time()

    print
    print "#######################"
    print "# Combine wit UniProt #"
    print "#######################"
    print
    print "This program is part of the PROTEOFORMER pipeline"
    print

    #Parse arguments from command line
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS,
                                     description="This tool is part of the PROTEOFORMER pipeline. It analyses "
                                                 "the different proteoforms found after MS analysis of a combined database "
                                                 "of PROTEOFORMER and UniProt. It searches for the classification of all "
                                                 "proteoforms found by PROTEOFORMER, not yet in UniProt and with confiramtion"
                                                 " of MS with MaxQuant.",
                                     add_help=False)

    man_args = parser.add_argument_group("Mandatory parameters")
    man_args.add_argument("--max_protein_db", "-m", action="store", required=True, nargs="?", metavar="DB",
                          type=str, help="Max protein DB (mandatory)")
    man_args.add_argument("--ens_db", "-e", action="store", required=True, nargs="?", metavar="DB",
                          type=str, help="Ensembl database (mandatory)")
    man_args.add_argument("--mapping_file", "-M", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="Uniprot mapping file. You can get that from the Uniprot website, downloads, "
                                         "ID mapping (mandatory)")
    man_args.add_argument("--fasta_file", "-f", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="Combined fasta file of proteoformer and uniprot (mandatory)")
    man_args.add_argument("--protein_groups", "-g", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="MaxQuant output file of protein groups (mandatory)")
    man_args.add_argument("--peptides", "-p", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="MaxQuant output file of peptides(mandatory)")

    opt_args = parser.add_argument_group("Optional parameters")
    opt_args.add_argument("--help", "-h", action="help", help="Show help message and exit")
    opt_args.add_argument("--workdir", "-w", action="store", required=False, nargs="?", metavar="FOLDER", default="",
                        type=str, help="Working directory (default: CWD)")
    opt_args.add_argument("--clustalo", "-c", action="store", required=False, nargs="?", metavar="PATH", default="clustalo",
                          type=str, help="Path to the Clustal Omega executable (default: clustalo)")
    opt_args.add_argument("--plot_file", "-o", action="store", required=False, nargs="?", metavar="PATH",
                          type=str, default="proteoform_class_plot.eps", help="Plot output file (default: proteoform_class_plot.eps)")
    opt_args.add_argument("--csv_output", "-x", action="store", required=False, nargs="?", metavar="PATH",
                          type=str, default="proteoform_classifications.csv", help="CSV output file of classifications (default: proteoform_classifications.csv)")

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

    #Load proteoformer-only identifications from max protein MS results. Herein, possible new proteoforms could be included.
    identifications = load_proteoformer_identifications(args.max_protein_db)

    #Make a Ensembl transcript id - Uniprot id mapping dict
    mapping_dict = construct_mapping_dict(args.mapping_file)

    #Get all sequences from the fasta file
    fasta_sequences = load_fasta(args.fasta_file)

    #Read in protein groups
    protein_groups = read_proteingroups(args.protein_groups)

    #Read in peptides
    peptides = read_peptides(args.peptides)

    ##############
    # Test cases #
    ##############
    #N-terminal extension
    #test_id = "ENST00000373382_10_69124218_5UTR_110db1;ENST00000395098_10_69124218_5UTR_100db1;ENST00000263559_10_69124230_5UTR_100db1;O75436;S4R3Q6;ENST00000489794_10_69132940_5UTR_100db1;ENST00000395098_10_69124278_aTIS_101db1;ENST00000395098_10_69124230_5UTR_100db1;ENST00000490696_10_69124218_ntr_100db1;ENST00000263559_10_69157162_CDS_100db1"
    #test_id = "ENST00000409021_4_6112850_aTIS_100db1;A0A1B0GUE0;F2Z2K5;Q96N16"
    #N-terminal truncation
    #test_id = "ENST00000338193_12_56751067_CDS_100db1;ENST00000338193_12_56743008_CDS_100db1"
    #test_id = "ENST00000281456_4_185143475_CDS_010db2"
    #N-terminal truncation without confirming peptide
    #test_id = "ENST00000593845_19_15397845_ntr_100db1;ENST00000397410_19_15397854_CDS_100db1;ENST00000397410_19_15401290_CDS_010db2;ENST00000595465_19_15418923_aTIS_101db1;ENST00000397410_19_15418923_aTIS_101db1;Q9ULX6"
    #NTR
    #test_id = "ENST00000506634_4_170605301_ntr_100db1"
    #test_id = "ENST00000428881_1_171683767_ntr_100db1"
    #3UTR
    #test_id = "ENST00000425660_7_5528626_3UTR_100db1"
    #Exon exclusion
    #test_id = "ENST00000233596_19_1491270_aTIS_111db1;ENST00000233596_19_1491222_5UTR_100db1;ENST00000233596_19_1491213_5UTR_100db1;ENST00000233596_19_1496294_CDS_100db1;ENST00000233596_19_1495575_CDS_100db1;Q96HR9;ENST00000395479_19_1496294_CDS_100db1;ENST00000395479_19_1495575_CDS_100db1;A8MXN1;A8MWX0;ENST00000591735_19_1491222_ntr_100db1;ENST00000591735_19_1491213_ntr_100db1"
    #uORF
    #test_id = "ENST00000419065_6_3258968_5UTR_100db1;ENST00000438998_6_3258968_5UTR_100db1;ENST00000473000_6_3258968_5UTR_100db1"
    #out of frame ORF
    #test_id = "ENST00000244020_20_43460226_CDS_100db1"
    #test_id = "ENST00000412249_3_51674477_CDS_100db1;ENST00000614067_3_51674477_CDS_100db1"
    #C terminal truncation
    #test_id="ENST00000398970_17_1456117_aTIS_101db1"
    #C terminal extension
    #test_id="ENST00000372798_1_40059347_aTIS_101db1;ENST00000340450_1_40040797_5UTR_110db1;ENST00000372798_1_40040797_5UTR_100db1;ENST00000340450_1_40060135_CDS_100db1;ENST00000479759_1_40067567_ntr_100db1;ENST00000479759_1_40069737_ntr_100db1;Q5T0R9;ENST00000479759_1_40069851_ntr_100db1;ENST00000476047_3_76434856_ntr_100db1;ENST00000417455_10_43605301_ntr_100db1;ENST00000476047_3_76434832_ntr_100db1;ENST00000479759_1_40070450_ntr_100db1;Q5T0S3;Q5T0R8;ENST00000493172_6_17421556_aTIS_101db1;B7Z385;A0A087X0J3;A0A087WZ15;ENST00000465994_6_17421556_aTIS_101db1;E9PDI2;P40123"
    #N terminal splice variant
    #test_id="ENST00000319410_16_75644428_aTIS_101db1;ENST00000302445_16_75641719_CDS_100db1;ENST00000302445_16_75647609_CDS_010db2;Q15046;ENST00000566560_16_75640313_ntr_100db1;H3BVA8;ENST00000566560_16_75641719_ntr_100db1;ENST00000569298_16_75628619_ntr_100db1;J3KRL2;H3BPV7;ENST00000566772_16_75636073_aTIS_101db1;H3BMR9;ENST00000566772_16_75635974_CDS_100db1;ENST00000566772_16_75635989_CDS_100db1;ENST00000566560_16_75635974_ntr_100db1;H3BQK5;ENST00000566560_16_75635989_ntr_100db1"
    #N-terminal truncation with SAV -> multiple variations
    #test_id = "ENST00000593845_19_15397845_ntr_100db1;ENST00000397410_19_15397854_CDS_100db1;ENST00000397410_19_15401290_CDS_010db2;ENST00000595465_19_15418923_aTIS_101db1;ENST00000397410_19_15418923_aTIS_101db1;Q9ULX6-2;Q9ULX6"

    #More complex variation
    #test_id = "ENST00000345136_8_143939461_aTIS_101db1;ENST00000398774_8_143944663_aTIS_101db1;ENST00000354958_8_143953771_aTIS_101db1;ENST00000356346_8_143973472_aTIS_101db1;ENST00000354589_8_143943890_aTIS_101db1;ENST00000357649_8_143942515_aTIS_101db1;ENST00000436759_8_143975369_aTIS_101db1;ENST00000527096_8_143975369_aTIS_101db1;E9PMV1;H0YDN1;E9PKG0;E9PIA2;E9PQ28;REV__A0A0U1RR03"



    #identifications_copy = defaultdict(lambda: defaultdict())
    #identifications_copy[test_id] = identifications[test_id]
    #identifications = identifications_copy

    #dict_funcs.print_dict(identifications[test_id])

    #Analyze identifications and parse them in the right categories
    identifications = analyze_identifications(args.workdir, identifications, mapping_dict, fasta_sequences, args.clustalo, protein_groups, peptides, args.ens_db)

    for protein_group in identifications:
        print identifications[protein_group]['base_proteoform']+": "+identifications[protein_group]['classification']

    #Output identifications and counts as csv
    output_csv(identifications, args.csv_output)


    #Count classifications
    counts = count_classifications(identifications)
    print
    print
    print "Counts"
    for classification in sorted(counts.keys()):
        print classification+": "+str(counts[classification])
    print

    #Consturct distribution plot
    construct_plot(counts, args.plot_file)



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

#Output identifications
def output_csv(identifications, csv_file):
	
	id=0
	with open(csv_file, 'w') as FW:
		keys = ["base_proteoform","classification","annotations","bin_codes","gene_ids","max_proteins","protein_group"]
		header_string= "id"+","+(','.join(keys))+"\n"
		FW.write(header_string)
		for protein_group in identifications:
			value_string = str(identifications[protein_group][keys[0]])
			for key in keys[1:]:
				if type(identifications[protein_group][key]) is list:
					value_string = value_string+','+str("|".join(identifications[protein_group][key]))
				else:
					if re.search(', ', str(identifications[protein_group][key])):
						value = re.sub(', ', '|', str(identifications[protein_group][key]))
						value_string = value_string+','+value
					else:
						value_string = value_string+','+str(identifications[protein_group][key])
			FW.write(str(id)+','+value_string+"\n")
			id+=1
	
	return

#Construct distribution plot
def construct_plot(counts, output_file):

    #Present classification distribution with plot
        #Without confirming peptide: weglaten
        #Multiple variations: maybe print as different lane
        #NTR bundelen
        #Exon inclusion/exclusion bundelen als splicing isoforms
        #Found with ensembl: drop

    #Parse into the counts you want to show
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
        elif classf=='Exon substitution':
            main_counts['Splice variants'] += counts[classf]
            splice_counts[classf] = counts[classf]
        elif classf=='Found with Ensembl check':
            pass
        elif m_ntr:
            main_counts['Translation in non-coding region'] += counts[classf]
            ntr_counts[m_ntr.group(1)] = counts[classf]
        else:
            main_counts[classf] = counts[classf]

    #Remove empty groups
    if main_counts['Splice variants']==0:
        del main_counts['Splice variants']
    if main_counts['Translation in non-coding region']==0:
        del main_counts['Translation in non-coding region']
    if splice_counts['Exon inclusion']==0:
        del splice_counts['Exon inclusion']
    if splice_counts['Exon exclusion']==0:
        del splice_counts['Exon exclusion']

    #Make plot
    fig = plt.figure(figsize=(30,20))
    # set up subplot grid
    GridSpec(2, 3)

    #Big plot
    # Reorder dataframe
    df = get_ordered_dataframe(main_counts)
    labels_list = df['count'].values.tolist()
    #Plot
    ax1 = plt.subplot2grid((2, 3), (0, 0), colspan=2, rowspan=2)
    #Generate explode tuple
    explode_tuple = tuple()
    for val in df.loc[:,'classification'].values.tolist():
        if val=='Translation in non-coding region' or val=='Splice variants':
            explode_tuple = explode_tuple + (0.2,)
        else:
            explode_tuple = explode_tuple + (0,)

    #Get colors
    colors = gen_color_list(df['classification'].values.tolist())
    ax1 = plt.subplot2grid((2, 3), (0, 0), colspan=2, rowspan=2)
    #Generate explode tuple
    explode_tuple = tuple()
    for val in df.loc[:,'classification'].values.tolist():
        if val=='Translation in non-coding region' or val=='Splice variants':
            explode_tuple = explode_tuple + (0.2,)
        else:
            explode_tuple = explode_tuple + (0,)
    cmap = plt.cm.tab10
    colors = cmap(np.linspace(0., 1., len(labels_list)))
    #Resturcture colors
    b = np.copy(colors[2,:])
    c = np.copy(colors[0,:])
    colors[2,:] = np.copy(c)
    colors[0,:] = np.copy(b)
    b = np.copy(colors[3, :])
    c = np.copy(colors[-1, :])
    colors[3, :] = np.copy(c)
    colors[-1, :] = np.copy(b)
    patches, texts, autotexts = ax1.pie(df['count'], explode=explode_tuple, autopct='%.1f%%', colors=colors, textprops={'fontsize':24}, pctdistance=0.85, labels = labels_list, labeldistance=1.05)
    #Label size
    for i in range(0, len(texts)):
        texts[i].set_fontsize(18)
    ax1.set_title("Overall proteoform categories", {'fontsize': 36})
    ax1.set_ylabel("")
    box = ax1.get_position()
    ax1.set_position([box.x0 - box.width * 0.1, box.y0,
                     box.width, box.height])
    lgd_labels = df['classification'].values.tolist()
    handles, labels = ax1.get_legend_handles_labels() #Get legend information of last ax object
    leg = ax1.legend(handles, lgd_labels, fontsize=24, loc='upper center', bbox_to_anchor = (.5, 0), ncol=3)

    #Small plot 1
    df = pd.DataFrame.from_dict(splice_counts, orient="index")
    df.columns = ['count']
    df['classification'] = df.index
    labels_list = df['count'].values.tolist()
    ax2 = plt.subplot2grid((2,3), (0,2))
    colors = get_greens(df['classification'].values.tolist())
    df.plot(kind="pie", subplots="True", autopct='%.1f%%', ax=ax2, colors=colors, textprops={'fontsize':24}, fontsize=24, pctdistance=0.7, labels = labels_list, labeldistance=1.05)
    ax2.set_ylabel("")
    ax2.set_title("Splice variants", {'fontsize': 36})
    lgd_labels = df['classification'].values.tolist()
    handles, labels = ax2.get_legend_handles_labels()  # Get legend information of last ax object
    leg = ax2.legend(handles, lgd_labels, fontsize=24, loc='upper center', bbox_to_anchor=(.5, 0.075), ncol=2)

    #Small plot 2
    df = pd.DataFrame.from_dict(ntr_counts, orient="index")
    df.columns = ['count']
    df['classification'] = df.index
    labels_list = df['count'].values.tolist()
    ax3 = plt.subplot2grid((2, 3), (1, 2))
    colors = get_reds(df['classification'].values.tolist())
    df.plot(kind="pie", subplots="True", autopct='%.1f%%', ax=ax3, colors=colors, textprops={'fontsize':24}, fontsize=24, pctdistance=0.7, labels = labels_list, labeldistance=1.05)
    ax3.set_ylabel("")
    ax3.set_title("Non-coding region proteoforms", {'fontsize':36})
    lgd_labels = rewrite_labels(df['classification'].values.tolist())
    handles, labels = ax3.get_legend_handles_labels()  # Get legend information of last ax object
    leg = ax3.legend(handles, lgd_labels, fontsize=24, loc='upper center', bbox_to_anchor=(.5, 0.075), ncol=2)

    #Connecting arrows
    (x,y) = texts[0].get_position()
    xy0 = (x+0.08, y+0.04)
    xyEnd = (-1.2, 0.5)
    con = ConnectionPatch(xyA=xy0, xyB=xyEnd, coordsA="data", coordsB="data",
                      axesA=ax1, axesB=ax2, arrowstyle=ArrowStyle.CurveFilledB(head_length=2, head_width=1), color='black', linewidth=4)
    ax1.add_artist(con)

    (x, y) = texts[-1].get_position()
    xy0 = (x + 0.06, y)
    xyEnd = (-1.2, 0.7)
    con = ConnectionPatch(xyA=xy0, xyB=xyEnd, coordsA="data", coordsB="data",
                          axesA=ax1, axesB=ax3, arrowstyle=ArrowStyle.CurveFilledB(head_length=2, head_width=1),
                          color='black', linewidth=4)
    ax1.add_artist(con)

    #Output
    fig.suptitle('Classification of PROTEOFORMER-only identified proteoforms', fontsize=54)
    fig.savefig(output_file, dpi=1200)

    return



#Get ordered data frame
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

#Generate color list
def gen_color_list(class_list):

    all_labels = ['Splice variants','C-terminal extension', 'C-terminal truncation', 'Multiple variations', 'N-terminal extension', 'N-terminal truncation', 'Only amino acid substitutions', 'Out of frame ORF', 'dORF', 'uORF', 'Translation in non-coding region']

    cmap = plt.cm.tab10
    colors = cmap(np.linspace(0., 1., len(all_labels)))

    #Resturcture colors
    b = np.copy(colors[2,:])
    c = np.copy(colors[0,:])
    colors[2,:] = np.copy(c)
    colors[0,:] = np.copy(b)
    b = np.copy(colors[3, :])
    c = np.copy([0.8588,0.8588,0.5529,1])
    colors[3, :] = np.copy(c)
    colors[-1, :] = np.copy(b)

    #Remove unneeded colors
    deleted_labels = 0
    for i in range(0, len(all_labels)):
        if all_labels[i] not in class_list:
            colors = np.delete(colors, i-deleted_labels, 0)
            deleted_labels+=1

    return colors
    
#Rewrite non coding labels
def rewrite_labels(labels):

    rewritten = list()
    translate = dict({'processed_transcript':'Processed transcript', 'transcribed_processed_pseudogene':'Transcribed processed pseudogene', 'processed_pseudogene': 'Processed pseudogene', 'retained_intron':'Retained intron'})

    for label in labels:
        rewritten.append(translate[label])

    return rewritten

#Get red colors
def get_reds(class_list):

    all_labels = ['processed_transcript', 'transcribed_processed_pseudogene', 'processed_pseudogene', 'retained_intron']

    cmap = plt.cm.Reds
    colors = cmap(np.linspace(.1, .9, len(all_labels)))

    #Remove unneeded colors
    deleted_labels = 0
    for i in range(0, len(all_labels)):
        if all_labels[i] not in class_list:
            colors = np.delete(colors, i-deleted_labels, 0)
            deleted_labels+=1

    return colors

#Get green colors
def get_greens(class_list):

    all_labels = ['Exon inclusion', 'Exon exclusion', 'Exon substitution', 'C-terminal splice variant', 'N-terminal splice variant']

    cmap = plt.cm.Greens
    colors = cmap(np.linspace(.1, .9, len(all_labels)))

    #Remove unneeded colors
    deleted_labels = 0
    for i in range(0, len(all_labels)):
        if all_labels[i] not in class_list:
            colors = np.delete(colors, i-deleted_labels, 0)
            deleted_labels+=1

    return colors

#Count classifications
def count_classifications(identifications):

    #Init
    counts = defaultdict()

    for protein_group in identifications:
        classification = identifications[protein_group]['classification']
        if classification in counts.keys():
            counts[classification]+=1
        else:
            counts[classification] = 1

    return counts

#Analyze identifications
def analyze_identifications(workdir, identifications, mapping_dict, fasta_sequences, clustalo_path, protein_groups, peptides, ens_db):

    #For all protein_groups
    for protein_group in identifications:
        #Get the most important proteoform
        base_proteoform = find_base_proteoform(identifications, protein_group)
        identifications[protein_group]['base_proteoform'] = base_proteoform
        #Get the most important canonical protein
        base_canonical, canonical_found_in_ensembl = find_base_canonical(protein_group, base_proteoform, mapping_dict)
        #Get sequences
        base_proteoform_seq = fasta_sequences[base_proteoform]
        #Get annotation of base proteoform
        base_proteoform_annotation = get_annotation(base_proteoform)

        #For non-coding: capture based on beforehand based on annotation
        if base_canonical=="non_coding":
            identifications[protein_group]['classification'] = "Translation in NTR"
            #Get more info about the non-coding transcript's biotype
            biotype = find_biotype(base_proteoform, ens_db)
            identifications[protein_group]['classification'] = identifications[protein_group]['classification'] + ", " + biotype
        #For dORFs: capture also based on annotation
        elif base_proteoform_annotation=="3UTR":
            #Do additional check based on coordinates of suspective dORF and canonical ORF
            #Get transcript id of base proteoform
            transcript_id = get_transcript_id(base_proteoform)
            #Get start of base proteoform
            base_proteoform_start = get_start(base_proteoform)
            #Find start and stop of canonical ORF of that transcript
            (canonical_start, canonical_stop, strand) = get_canonical_start_stop(transcript_id, ens_db)
            if strand==1:
                if base_proteoform_start>canonical_stop:
                    identifications[protein_group]['classification'] = "dORF"
            elif strand==-1:
                if base_proteoform_start<canonical_stop:
                    identifications[protein_group]['classification'] = "dORF"
        #For extensions and truncations
        else:
            if base_canonical not in fasta_sequences:
                print "Base canonical not found: "+base_canonical+" for protein group "+protein_group
            base_canonical_seq = fasta_sequences[base_canonical]
            #Align with ClustalO
            #Prepare fasta
            tmp_fa = prepare_fasta(workdir, base_proteoform, base_proteoform_seq, base_canonical, base_canonical_seq)
            #Align
            output_fa = clustal_align(workdir, base_proteoform, tmp_fa, clustalo_path)
            #Read aligned sequences
            (aligned_base_proteoform, aligned_base_canonical) = read_clustal_output(output_fa, base_proteoform, base_canonical)
            #Remove fasta
            os.system("rm -rf "+tmp_fa)
            #Calculate mapped percentage
            perc_mapped = calc_perc_mapped(aligned_base_proteoform, aligned_base_canonical)
            #Get peptides and their position on the full protein sequence
            proteoform_peptides = get_peptides(protein_groups, peptides, protein_group)
            #Classify proteoform
            identifications[protein_group]['classification'] = classify_proteoform(proteoform_peptides,
                aligned_base_proteoform, aligned_base_canonical, perc_mapped, base_proteoform, base_canonical,
                canonical_found_in_ensembl, ens_db)
            #remove clustal output fasta
            #os.system("rm -rf "+output_fa)

        #print protein_group
        #print "Classification: "+identifications[protein_group]['classification']
        #print



        if identifications[protein_group]['classification'] == '':
            print
            print
            print "-----UNCLASSIFIED!!-------"
            print "Base proteoform: "+base_proteoform
            print "Base canonical: "+base_canonical
            print "Protein group: "+protein_group
            print "--------------------------"
            print
            print

    return identifications

#Classify proteoform
def classify_proteoform(proteoform_peptides, aligned_base_proteoform, aligned_base_canonical, perc_mapped, base_proteoform, base_canonical, canonical_found_in_ensembl, ens_db):

    #Init
    classification=""

    #Check if enough identical positions
    if perc_mapped>0.3:
        #Correct start and stop positions based on alignment
        for peptide in proteoform_peptides:
            proteoform_peptides[peptide]['corr_start'] = correct_pos(aligned_base_proteoform, proteoform_peptides[peptide]['start'])
            proteoform_peptides[peptide]['corr_end'] = correct_pos(aligned_base_proteoform, proteoform_peptides[peptide]['end'])

        #If N-terminal extension: proteoform starts with AA's and canonical with '-'s. There is also a corr peptide start in this stretch
        m1 = re.search('^[A-Za-z]+',  aligned_base_proteoform)
        m2 = re.search('^(\-+)', aligned_base_canonical)
        #If N-terminal truncation: proteoform starts with '-'s and canonical with AA's. Peptide at start or second position of truncation
        m3 = re.search('^(\-+)', aligned_base_proteoform)
        m4 = re.search('^[A-Za-z]+', aligned_base_canonical)
        #If C-terminal extension: proteoform ends with AA's and canonical with '-'s. There is a corr_peptide end in this stretch
        m5 = re.search('[A-Za-z]+$', aligned_base_proteoform)
        m6 = re.search('(\-+)$', aligned_base_canonical)
        #If C-terminal truncation: proteoform ends with '-'s and canonical with AA's. Peptide found at last or second last position of proteoform
        m7 = re.search('(\-+)$', aligned_base_proteoform)
        m8 = re.search('[A-Za-z]+$', aligned_base_canonical)
        #If exon exclusion: insert in proteoform
        m9 = re.search('^([A-Za-z]+)\-+([A-Za-z]+)$', aligned_base_proteoform)
            #multiple
        m9b = re.search('^[A-Za-z]+\-+[A-Za-z]+\-+[A-Za-z]+$', aligned_base_proteoform)
        m9c = re.search('^[A-Za-z]+\-+[A-Za-z]+\-+[A-Za-z]+\-+[A-Za-z]+$', aligned_base_proteoform)
        #If exon inclusion: insert in canonical
        m10 = re.search('^([A-Za-z]+)\-+([A-Za-z]+)$', aligned_base_canonical)
            #multiple
        m10b = re.search('^[A-Za-z]+\-+[A-Za-z]+\-+[A-Za-z]+$', aligned_base_canonical)
        m10c = re.search('^[A-Za-z]+\-+[A-Za-z]+\-+[A-Za-z]+\-+[A-Za-z]+$', aligned_base_canonical)
        #If there are no indels from start to stop in both sequence, then only SAVs
        m11 = re.search('^[A-Za-z]+$', aligned_base_proteoform)
        m12 = re.search('^[A-Za-z]+$', aligned_base_canonical)
        #Sequence absent in Uniprot
        if (aligned_base_proteoform==aligned_base_canonical):
            classification = "Sequence absent in UniProt"
        #N-terminal extension
        if (m1 and m2):
            # Get the rest of the sequences (-initiator methionine) and calculate similarity
            rest_of_proteoform = aligned_base_proteoform[(len(m2.group(1)) + 1):]
            rest_of_canonical = aligned_base_canonical[(len(m2.group(1)) + 1):]
            # Check if the rest of the sequences is identical
            if rest_of_proteoform == rest_of_canonical:
                peptide_found = False
                for peptide in proteoform_peptides:
                    if (proteoform_peptides[peptide]['corr_start'] <= len(m2.group(1))):
                        #Check if peptide is also in a canonical
                        peptide_in_canonical = False
                        for protein in proteoform_peptides[peptide]['proteins'].split(';'):
                            #if not re.search('^ENST', protein): #To be very stringent, look also for presence in other canonical proteins
                            if protein==base_canonical:
                                peptide_in_canonical = True
                        if peptide_in_canonical == False:
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
                    if rest_of_proteoform.endswith(longest_stretch) and perc>0.1:
                        classification = "N-terminal splice variant"
                    else:
                        classification = "Multiple variations"
        # C-terminal extension
        elif (m5 and m6):
            # Get the rest of the sequences and calculate similarity
            rest_of_proteoform = aligned_base_proteoform[:(len(aligned_base_proteoform)-len(m6.group(1)))]
            rest_of_canonical = aligned_base_canonical[:(len(aligned_base_canonical)-len(m6.group(1)))]
            # Check if the rest of the sequences is identical
            if rest_of_proteoform == rest_of_canonical:
                peptide_found = False
                for peptide in proteoform_peptides:
                    if (proteoform_peptides[peptide]['corr_end'] >= (len(aligned_base_proteoform) - len(m6.group(1)))):
                        # Check if peptide is also in a canonical
                        peptide_in_canonical = False
                        for protein in proteoform_peptides[peptide]['proteins'].split(';'):
                            #if not re.search('^ENST', protein):
                            if protein==base_canonical:
                                peptide_in_canonical = True
                        if peptide_in_canonical == False:
                            peptide_found = True
                if peptide_found:
                    classification = "C-terminal extension"
                else:
                    classification = 'C-terminal extension, without confirming peptide'
            else:
                if check_only_sav(rest_of_proteoform, rest_of_canonical):
                    classification = "Multiple variations"
                else:
                    #Search longest identical stretch of characters
                    longest_stretch, perc = search_longest_stretch(rest_of_proteoform, rest_of_canonical)
                    #Check if stretch can be found at the start
                    if rest_of_proteoform.startswith(longest_stretch) and perc>0.1:
                        classification = "C-terminal splice variant"
                    else:
                        classification = "Multiple variations"
        #N-terminal truncation
        elif (m3 and m4):
            # Get the rest of the sequences (-initiator methionine) and calculate similarity
            rest_of_proteoform = aligned_base_proteoform[(len(m3.group(1)) + 1):]
            rest_of_canonical = aligned_base_canonical[(len(m3.group(1)) + 1):]
            # Check if the rest of the sequences is identical
            if rest_of_proteoform == rest_of_canonical:
                #Check for peptide
                peptide_found = False
                for peptide in proteoform_peptides:
                    if (proteoform_peptides[peptide]['corr_start'] in [len(m3.group(1)) + 1, len(m3.group(1)) + 2]):
                        # Check if peptide is also in a canonical
                        peptide_in_canonical = False
                        for protein in proteoform_peptides[peptide]['proteins'].split(';'):
                            #if not re.search('^ENST', protein):
                            if protein==base_canonical:
                                peptide_in_canonical=True
                        if peptide_in_canonical == False:
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
        #C-terminal truncation
        elif (m7 and m8):
            # Get the rest of the sequences and calculate similarity
            rest_of_proteoform = aligned_base_proteoform[:(len(aligned_base_proteoform)-len(m7.group(1)))]
            rest_of_canonical = aligned_base_canonical[:(len(aligned_base_proteoform)-len(m7.group(1)))]
            # Check if the rest of the sequences is identical
            if rest_of_proteoform == rest_of_canonical:
                #Check for peptide
                peptide_found = False
                for peptide in proteoform_peptides:
                    if (proteoform_peptides[peptide]['corr_end'] in [len(aligned_base_proteoform) - len(m7.group(1)),
                                                                     len(aligned_base_proteoform) - len(
                                                                             m7.group(1)) - 1]):
                        # Check if peptide is also in a canonical
                        peptide_in_canonical = False
                        for protein in proteoform_peptides[peptide]['proteins'].split(';'):
                            #if not re.search('^ENST', protein):
                            if protein==base_canonical:
                                peptide_in_canonical = True
                        if peptide_in_canonical == False:
                            peptide_found = True
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
        #Exon exclusion
        elif (m9):
            regex = r"^"+re.escape(m9.group(1))+r"[A-Za-z]+?"+re.escape(m9.group(2))+r"$"
            m_excl = re.search(regex, aligned_base_canonical)
            if m_excl:
                classification = "exon exclusion"
            else:
                classification = "exon exclusion, more complex"
        #Exon inclusion
        elif (m10):
            regex = r"^"+re.escape(m10.group(1))+r"[A-Za-z]+?"+re.escape(m10.group(2))+r"$"
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
        #Only amino acid variations
        elif (m11 and m12):
            if aligned_base_canonical==aligned_base_proteoform:
                if canonical_found_in_ensembl:
                    classification = "Found with Ensembl check"
            else:
                #Check only sav
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




    #More complex variants (multiple possibilities): ENST00000345136_8_143939461_aTIS_101db1 (plectin)
        #Classified as a N-terminal truncation but is actually a C-terminal extension



    #Out of frame ORF?
    #SAV?





    return classification

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

#Get start of proteoform
def get_start(proteoform):

    #Init
    start = ""

    m = re.search('^ENST\d+?\_.+?\_(\d+?)\_', proteoform)
    if m:
        start = int(m.group(1))

    return start

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

#Find annotation of accession
def get_annotation(proteoform):

    #Init
    annotation = ""

    m = re.search('^ENST\d+?\_.+?\_\d+?\_(.+?)\_.+', proteoform)
    if m:
        annotation = m.group(1)

    return annotation

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

#Get peptides and their positions on the fulle protein sequence
def get_peptides(protein_groups, peptides, protein_group):

    #Init
    proteoform_peptides = defaultdict(lambda: defaultdict())

    #Get all peptides of protein group
    peptide_list = protein_groups[protein_group]['Peptide IDs'].split(';')

    #For all peptides
    for peptide in peptide_list:
        #Get start, stop position and sequence
        if peptide in peptides.keys():#Otherwise it is a contaminant or reverse peptide
            proteoform_peptides[peptide]['sequence'] = peptides[peptide]['Sequence']
            proteoform_peptides[peptide]['start'] = peptides[peptide]['Start position']
            proteoform_peptides[peptide]['end'] = peptides[peptide]['End position']
            proteoform_peptides[peptide]['proteins'] = peptides[peptide]['Proteins']

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
def clustal_align(workdir, base_proteoform, input_fa, clustalo_path):

    #Init output
    output_fa = workdir+"/"+base_proteoform+"aligned_output.fa"

    #Execute
    command = clustalo_path+" -i "+input_fa+" -o "+output_fa+" --outfmt fa --force --wrap=1000000000" #Force to overwrite, output in fasta format
    os.system(command)

    return output_fa

#Prepare tmp fasta
def prepare_fasta(workdir, base_proteoform, base_proteoform_seq, base_canonical, base_canonical_seq):

    #Define path of tmp fasta file
    tmp_fasta = workdir+"/tmp_input.fa"

    with open(tmp_fasta, 'w') as FW:
        FW.write(">"+base_proteoform+"\n")
        FW.write(base_proteoform_seq+"\n")
        FW.write(">"+base_canonical+"\n")
        FW.write(base_canonical_seq+"\n")

    return tmp_fasta

#Find the base canonical protein
def find_base_canonical(protein_group, base_proteoform, mapping_dict):

    #Init
    base_canonical = ""

    #try to find it in the protein group
    #Split by ;
    proteins = map(lambda x: str(x), protein_group.split(';'))
    for protein in proteins:
        if not re.search('^ENST', protein):
            return protein, False
    #Get the annotation of the base proteoform
    m = re.search('^ENST\d+?\_.+?\_\d+?\_(.+?)\_.+', base_proteoform)
    if m:
        #Check if the base proteoform is non-coding
        if m.group(1)=="ntr":
            return "non_coding", False
    #Check if base canonical already found in protein group
    if base_canonical=="":
        #If not found there, search the base canonical protein in ensembl db with the transcript id of the base proteoform
        base_canonical = search_base_canonical_with_mapping_info(base_proteoform, mapping_dict)

    return base_canonical, True

#Try to find the base canonical with ensembl info of the base proteoform
def search_base_canonical_with_mapping_info(base_proteoform, mapping_dict):

    #Init
    base_canonical=""

    #Get the transcript id of the base proteoform
    transcript_id = ""
    m = re.search('^(ENST.+?)\_', base_proteoform)
    if m:
        transcript_id=m.group(1)

    if transcript_id in mapping_dict.keys():
        base_canonical = mapping_dict[transcript_id]

    return base_canonical

#Find the base proteoform
def find_base_proteoform(identifications, protein_group):

    #If multiple max proteins
    if isinstance(identifications[protein_group]['max_proteins'], list):
        base_proteoform = identifications[protein_group]['max_proteins'][0] #Take for now the first proteoform as base proteoform
    elif isinstance(identifications[protein_group]['max_proteins'], str) or isinstance(identifications[protein_group]['max_proteins'], unicode):
        base_proteoform = identifications[protein_group]['max_proteins']

    #print "Base proteoform: "+base_proteoform

    return base_proteoform

#Read peptides input info from maxquant
def read_peptides(in_file):

    #Init
    peptides = defaultdict(lambda: defaultdict())

    #Open file
    with open(in_file, 'r') as FR:

        #Read lines
        lines = FR.readlines()
        lines = map(lambda x: x.rstrip("\n"), lines)
        lines = map(lambda x: x.rstrip("\r"), lines)

        #Parse column headers
        col_headers = re.split('\t', lines[0])

        #Search index of contaminants and reverses
        idx_reverses = col_headers.index("Reverse")
        idx_cont = col_headers.index("Potential contaminant")

        #Parse data
        for line_nr in range(1,len(lines)):
            data = re.split('\t', lines[line_nr])
            #Remove contaminants and reverse identifications
            if(data[idx_reverses]!='+' and data[idx_cont]!='+'):
                for val_nr in range(0, len(data)):
                    peptides[data[48]][col_headers[val_nr]] = data[val_nr] #data[48] is peptide ID

    return peptides

#Read protein group input info from maxquant
def read_proteingroups(in_file):

    #Init
    protein_groups = defaultdict(lambda: defaultdict())

    #Open file
    with open(in_file, 'r') as FR:

        #Read lines
        lines = FR.readlines()
        lines = map(lambda x: x.rstrip("\n"), lines)
        lines = map(lambda x: x.rstrip("\r"), lines)

        #Parse column headers
        col_headers = re.split('\t', lines[0])

        #Search index of contaminants and reverses
        idx_reverses = col_headers.index("Reverse")
        idx_cont = col_headers.index("Potential contaminant")

        #Parse data
        for line_nr in range(1,len(lines)):
            data = re.split('\t', lines[line_nr])
            #Remove contaminants and reverse identifications
            if(data[idx_reverses]!='+' and data[idx_cont]!='+'):
                for val_nr in range(0, len(data)):
                    protein_groups[data[0]][col_headers[val_nr]] = data[val_nr]

    return protein_groups

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

#Load proteoformer-only identifications
def load_proteoformer_identifications(max_protein_db):

    #Init
    #identifications->{protein_group}->{column} = value
    identifications = defaultdict(lambda: defaultdict())

    #DB connection
    try:
        conn = sqlite.connect(max_protein_db)
    except:
        print "Could not connect to "+max_protein_db
        sys.exit()

    with conn:
        cur = conn.cursor()

        #Get all protein groups only identified with a max protein by Proteoformer
        query = "SELECT protein_group, max_proteins, annotations, bin_codes, gene_ids FROM max_proteins WHERE sources='Proteoformer';"
        cur.execute(query)

        colnames = list(map(lambda x: x[0], cur.description))

        data = cur.fetchall()
        for i in range(0, len(data)):
            id = data[i][0]
            for j in range(0, len(colnames)):
                col = colnames[j]
                value = data[i][j]
                if re.search('\|', value):
                    identifications[id][col] = map(lambda x: str(x), value.split('|'))
                else:
                    identifications[id][col] = value

    return identifications

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