#!/usr/bin/python

# Merge pre-split SPECtre experiments into a single output. User must provide the directory of
# split SPECtre experiments to be merged, and then the columns that refer to the SPECtre_score,
# Full_score (if necessary), and FLOSS_distribution (if necessary) in the split experiment 
# files. Based on these inputs, the posterior probability of translation and/or the FLOSS score
# for each transcript will be calculated.

# Usage:
# python merge_split_files.py \
# <split-experiment-directory> \
# <cufflinks-file> \
# --SPEC <int> \
# --FULL <int> \
# --FLOSS <int> > merged_experiments.txt

# IMPORT STANDARD PYTHON LIBRARIES:
import sys
import os
import re
import math
import logging
import argparse
import itertools
import collections

from functools import partial
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool

# IMPORT THIRD-PARTY DEPENDENCIES:
import HTSeq
import numpy
import pysam
import pyfasta

from rpy2.robjects import r
from rpy2.robjects.packages import importr

#############
# UTILITIES #
#############
def hash():
	return collections.defaultdict(hash)

def flatten(d, parent_key='', sep=':'):
    items = []
    for k, v in d.items():
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, collections.MutableMapping):
            items.extend(flatten(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)

def parse_experiments(directory, fpkm_cutoff, spec_column, full_column, floss_column, orf_column):

	'''
	Takes as input a directory that contains the split experiments files to be merged.
	User-defined columns that correspond to the SPEC_score_CDS, FULL_score_CDS, and
	FLOSS_distribution_CDS are parsed to build the relevant distributions to compute
	the posteriors and FLOSS_score for each transcript.
	'''

	reference_distribution = dict(zip(range(24,37), [0.0]*len(range(24,37))))
	not_translated_spec_scores = list()
	not_translated_full_scores = list()
	not_translated_orf_scores = list()
	translated_spec_scores = list()
	translated_full_scores = list()
	translated_orf_scores = list()
	transcripts = hash()
	fields = dict()
	count = 0
	flag = False

	for infile in os.listdir(directory):
		logger.info("parse_experiments(): Parsing experiment file: " + infile + "... [STARTED]")
		if "txt" in infile:	# Change to match the file extension used to output split experiments.
			for linenum, line in enumerate(open(directory + infile)):
				if not line.startswith("#"): 
					if line.startswith("id"):	# Header line.
						fields = dict(zip([str(field) for field in line.strip().split("\t")], range(0, len(line.strip().split("\t")))))
						if flag == False:
							print line.strip()
							flag = True
						else:
							pass
					elif line.startswith("["):	# FLOSS reference distribution.
						pass
					else:
						transcript_line = line.strip().split("\t")
						if "ribo_fpkm" in fields:
							if spec_column > 0:
								if "SPEC_score_CDS" in fields:
									try:
										if float(transcript_line[fields["ribo_fpkm"]]) >= float(fpkm_cutoff):
											try:
												translated_spec_scores.append(float(transcript_line[spec_column]))
											except ValueError:
												pass
									except ValueError:
										pass
									try:
										if float(transcript_line[fields["ribo_fpkm"]]) < float(fpkm_cutoff):
											try:
												not_translated_spec_scores.append(float(transcript_line[spec_column]))
											except ValueError:
												pass
									except ValueError:
										pass
							if full_column > 0:
								if "FULL_score_CDS" in fields:
									try:
										if float(transcript_line[fields["ribo_fpkm"]]) >= float(fpkm_cutoff):
											try:
												translated_full_scores.append(float(transcript_line[full_column]))
											except ValueError:
												pass
									except ValueError:
										pass
									try:
										if float(transcript_line[fields["ribo_fpkm"]]) < float(fpkm_cutoff):
											try:
												not_translated_full_scores.append(float(transcript_line[full_column]))
											except ValueError:
												pass
									except ValueError:
										pass
							if orf_column > 0:
								if "ORF_score_CDS" in fields:
									try:
										if float(transcript_line[fields["ribo_fpkm"]]) >= float(fpkm_cutoff):
											try:
												translated_orf_scores.append(float(transcript_line[orf_column]))
											except ValueError:
												pass
									except ValueError:
										pass
									try:
										if float(transcript_line[fields["ribo_fpkm"]]) < float(fpkm_cutoff):
											try:
												not_translated_orf_scores.append(float(transcript_line[orf_column]))
											except ValueError:
												pass
									except ValueError:
										pass
						if floss_column > 0:
							if "FLOSS_distribution_CDS" in fields:
								dist = transcript_line[fields["FLOSS_distribution_CDS"]]
								if not dist == "NA":
									lengths, floss = [x for x in re.split("\[|\(|\)|\],", dist) if len(x) > 2]
									floss_dist_cds = dict(zip([int(l) for l in lengths.split(", ")], [float(f) for f in floss.split(", ")]))
									for read_length in floss_dist_cds:
										if read_length in reference_distribution:
											reference_distribution[read_length] += floss_dist_cds[read_length]
									count += 1
						transcripts[transcript_line[fields["gene_type"]]][transcript_line[fields["chr"]]][transcript_line[fields["strand"]]][transcript_line[fields["gene_id"]]][transcript_line[fields["transcript_id"]]] = line.strip()
		logger.info("parse_experiments(): Parsing experiment file: " + infile + "... [COMPLETE]")
	return transcripts, fields, not_translated_spec_scores, translated_spec_scores, not_translated_full_scores, translated_full_scores, not_translated_orf_scores, translated_orf_scores, reference_distribution, count

# TRANSCRIPT-LEVEL METRICS:
def calculate_reference_distribution(aggregate_distribution, transcript_count):
	return {k: v * (1/float(transcript_count)) for k, v in aggregate_distribution.items()}

def calculate_floss_score(reference_distribution, transcript_distribution):
	if isinstance(transcript_distribution, dict):
		if sum(transcript_distribution.values()) > 0:
			floss = 0
			for read_length in reference_distribution:
				if read_length in transcript_distribution:
					floss += abs(transcript_distribution[read_length] - reference_distribution[read_length])
				else:
					floss += abs(0.0 - float(reference_distribution[read_length]))
			return floss / 2.0
		else:
			return "NA"
	else:
		return "NA"

def calculate_posterior_probability(coding_scores, noncoding_scores, score):
	if score == "NA" or score == 0.0:
		return 0.0
	elif isinstance(score, float) or isinstance(score, int):
		if score > max(coding_scores):
			return 1.0
		else:
			prob_score_coding = len([n for n in coding_scores if n >= score])/float(len(coding_scores))
			prob_coding = len(coding_scores)/float(len(coding_scores)+len(noncoding_scores))
			prob_score = len([n for n in coding_scores+noncoding_scores if n >= score])/float(len(coding_scores)+len(noncoding_scores))
		return (prob_score_coding * prob_coding) / prob_score
	else:
		return 0.0

def posterior_probability(coding_scores, noncoding_scores, transcript_score):
	if transcript_score in ("NA", ""):
		return "NA"
	else:
		if "," in transcript_score:
			scores = [float(n) for n in transcript_score.split(",")]
			windowed = list()
			for coh in scores:
				windowed.append(calculate_posterior_probability(coding_scores, noncoding_scores, coh))
			return ",".join([str(n) for n in windowed])
		else:
			try:
				if len(transcript_score) > 0:
					score = float(transcript_score)
					return calculate_posterior_probability(coding_scores, noncoding_scores, score)
				else:
					return "NA"
			except ValueError:
				return "NA"

# EXPERIMENTAL METRICS:
def calculate_translation_threshold(translated, not_translated, fdr, spec_column):
	logger.info("ExperimentMetrics.translation_threshold(): Calculating experiment-level translation threshold [STARTED].")
	if spec_column > 0:
		r('active <- data.frame(SPEC=c(%s), biotype="translated")' %",".join(str(n) for n in translated))
		r('inactive <- data.frame(SPEC=c(%s), biotype="not_translated")' %",".join(str(n) for n in not_translated))
		r('scores <- rbind(active, inactive)')
		return str(r('quantile(scores$SPEC[scores$biotype=="not_translated"], probs=%s)' %(1-float(fdr)))[0])
	else:
		return "NA"
	logger.info("ExperimentMetrics.translation_threshold(): Calculating experiment-level translation threshold [COMPLETE].")

def calculate_spectre_auc(translated, not_translated, spec_column):
	logger.info("ExperimentMetrics.spectre_auc(): Calculating experiment-level SPECtre AUC [STARTED].")
	if spec_column > 0:
		# Instantiate necessary R packages:
		rocr = importr("ROCR")
		r('active <- data.frame(SPEC=c(%s), biotype="translated")' %",".join(str(n) for n in translated))
		r('inactive <- data.frame(SPEC=c(%s), biotype="not_translated")' %",".join(str(n) for n in not_translated))
		r('scores <- rbind(active, inactive)')
		return str(r('performance(prediction(scores$SPEC, scores$biotype), "auc")@y.values[[1]]')[0])
	else:
		return "NA"
	logger.info("ExperimentMetrics.spectre_auc(): Calculating experiment-level SPECtre AUC [COMPLETE].")

def calculate_floss_auc(translated, not_translated, floss_column):
	logger.info("ExperimentMetrics.floss_auc(): Calculating experiment-level FLOSS AUC [STARTED].")
	if floss_column > 0:
		# Instantiate necessary R packages:
		rocr = importr("ROCR")
		r('active <- data.frame(SPEC=c(%s), biotype="translated")' %",".join(str(n) for n in translated))
		r('inactive <- data.frame(SPEC=c(%s), biotype="not_translated")' %",".join(str(n) for n in not_translated))
		r('scores <- rbind(active, inactive)')
		return str(r('performance(prediction(-scores$SPEC, scores$biotype), "auc")@y.values[[1]]')[0])
	else:
		return "NA"
	logger.info("ExperimentMetrics.floss_auc(): Calculating experiment-level FLOSS AUC [COMPLETE].")

def calculate_orfscore_auc(translated, not_translated, orf_column):
	logger.info("ExperimentMetrics.orfscore_auc(): Calculating experiment-level ORFscore AUC [STARTED].")
	if orf_column > 0:
		# Instantiate necessary R packages:
		rocr = importr("ROCR")
		r('active <- data.frame(SPEC=c(%s), biotype="translated")' %",".join(str(n) for n in translated))
		r('inactive <- data.frame(SPEC=c(%s), biotype="not_translated")' %",".join(str(n) for n in not_translated))
		r('scores <- rbind(active, inactive)')
		return str(r('performance(prediction(scores$SPEC, scores$biotype), "auc")@y.values[[1]]')[0])
	else:
		return "NA"
	logger.info("ExperimentMetrics.orfscore_auc(): Calculating experiment-level ORFscore AUC [COMPLETE].")

def calculate_full_auc(translated, not_translated, full_column):
	logger.info("ExperimentMetrics.full_auc(): Calculating experiment-level Spectral Coherence AUC [STARTED].")
	if full_column > 0:
		# Instantiate necessary R packages:
		rocr = importr("ROCR")
		r('active <- data.frame(SPEC=c(%s), biotype="translated")' %",".join(str(n) for n in translated))
		r('inactive <- data.frame(SPEC=c(%s), biotype="not_translated")' %",".join(str(n) for n in not_translated))
		r('scores <- rbind(active, inactive)')
		return str(r('performance(prediction(scores$SPEC, scores$biotype), "auc")@y.values[[1]]')[0])
	else:
		return "NA"
	logger.info("ExperimentMetrics.full_auc(): Calculating experiment-level Spectral Coherence AUC [COMPLETE].")

if __name__ == "__main__":
	
	# Parse arguments from command-line:
	parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS, add_help=False)
	parser.add_argument("--help", action="help", help="show this help message and exit")
	parser.add_argument("--spec", action="store", required=True, nargs="?", metavar="INT", default="0", help="column number that corresponds to SPEC_score_CDS")
	parser.add_argument("--full", action="store", required=False, nargs="?", metavar="INT", default="0", help="column number that corresponds to FULL_score_CDS")
	parser.add_argument("--floss", action="store", required=False, nargs="?", metavar="INT", default="0", help="column number that corresponds to FLOSS_distribution_CDS")
	parser.add_argument("--orfscore", action="store", required=False, nargs="?", metavar="INT", default="0", help="column number that corresponds to ORF_score_CDS")
	spectre_args = parser.add_argument_group("parameters for SPECtre analysis:")
	spectre_args.add_argument("--min", action="store", required=False, nargs="?", default=3, metavar="FLOAT", help="minimum FPKM for active translation (default: %(default)s FPKM)")
	spectre_args.add_argument("--fdr", action="store", required=False, nargs="?", default=0.05, metavar="FLOAT", help="FDR cutoff (default: %(default)s)")
	file_args = parser.add_argument_group("input and output parameters:")
	file_args.add_argument("--dir", action="store", required=False, nargs="?", default="./", metavar="DIR", help="location of split SPECtre results files (default: %(default)s)")
	file_args.add_argument("--log", action="store", required=False, nargs="?", default="./spectre_results.log", metavar="FILE", help="track progress to (default: %(default)s)")
	args = parser.parse_args()

	# Enable event logging:
	logger = logging.getLogger(__name__)
	logger.setLevel(logging.INFO)
	handler = logging.FileHandler(args.log)
	handler.setLevel(logging.INFO)
	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	handler.setFormatter(formatter)
	logger.addHandler(handler)

	transcripts, fields, not_translated_spec_scores, translated_spec_scores, not_translated_full_scores, translated_full_scores, not_translated_orf_scores, translated_orf_scores, reference_distribution, count = parse_experiments(args.dir, float(args.min), int(args.spec), int(args.full), int(args.floss), int(args.orfscore))

	reference_floss_distribution = calculate_reference_distribution(reference_distribution, count)

	not_translated_floss_scores = list()
	translated_floss_scores = list()
	id_count = 1
	logger.info("__main__: Calculating transcript-level metrics... [STARTED]")
	for gene_type in sorted(transcripts):
		for chrom in sorted(transcripts[gene_type]):
			for strand in sorted(transcripts[gene_type][chrom]):
				for gene_id in sorted(transcripts[gene_type][chrom][strand]):
					for transcript_id in sorted(transcripts[gene_type][chrom][strand][gene_id]):
						line = transcripts[gene_type][chrom][strand][gene_id][transcript_id].strip().split("\t")
						if int(args.spec) > 0:
							spec_score_posterior_utr5 = posterior_probability(translated_spec_scores, not_translated_spec_scores, line[fields["SPEC_score_UTR5"]])
							spec_score_posterior_utr3 = posterior_probability(translated_spec_scores, not_translated_spec_scores, line[fields["SPEC_score_UTR3"]])
							spec_score_posterior_cds = posterior_probability(translated_spec_scores, not_translated_spec_scores, line[fields["SPEC_score_CDS"]])
							spec_signal_posterior_utr5 = posterior_probability(translated_spec_scores, not_translated_spec_scores, line[fields["SPEC_signal_UTR5"]])
							spec_signal_posterior_utr3 = posterior_probability(translated_spec_scores, not_translated_spec_scores, line[fields["SPEC_signal_UTR3"]])
							spec_signal_posterior_cds = posterior_probability(translated_spec_scores, not_translated_spec_scores, line[fields["SPEC_signal_CDS"]])
						if int(args.full) > 0:			
							full_score_posterior_utr5 = posterior_probability(translated_spec_scores, not_translated_spec_scores, line[fields["SPEC_score_UTR5"]])
							full_score_posterior_utr3 = posterior_probability(translated_spec_scores, not_translated_spec_scores, line[fields["SPEC_score_UTR3"]])
							full_score_posterior_cds = posterior_probability(translated_spec_scores, not_translated_spec_scores, line[fields["SPEC_score_CDS"]])
						if int(args.floss) > 0:
							if line[fields["FLOSS_distribution_UTR5"]] in ("NA", ""):
								floss_score_utr5 = "NA"
							else:
								lengths_utr5, floss_utr5 = [x for x in re.split("\[|\(|\)|\],", line[fields["FLOSS_distribution_UTR5"]]) if len(x) > 2]
								floss_distribution_utr5 = dict(zip([int(l) for l in lengths_utr5.split(", ")], [float(f) for f in floss_utr5.split(", ")]))
								floss_score_utr5 = calculate_floss_score(reference_floss_distribution, floss_distribution_utr5)			
							if line[fields["FLOSS_distribution_UTR3"]] in ("NA", ""):
								floss_score_utr3 = "NA"
							else:
								lengths_utr3, floss_utr3 = [x for x in re.split("\[|\(|\)|\],", line[fields["FLOSS_distribution_UTR3"]]) if len(x) > 2]
								floss_distribution_utr3 = dict(zip([int(l) for l in lengths_utr3.split(", ")], [float(f) for f in floss_utr3.split(", ")]))
								floss_score_utr3 = calculate_floss_score(reference_floss_distribution, floss_distribution_utr3)			
							if line[fields["FLOSS_distribution_CDS"]] in ("NA", ""):	
								floss_score_cds = "NA"
							else:
								lengths_cds, floss_cds = [x for x in re.split("\[|\(|\)|\],", line[fields["FLOSS_distribution_CDS"]]) if len(x) > 2]
								floss_distribution_cds = dict(zip([int(l) for l in lengths_cds.split(", ")], [float(f) for f in floss_cds.split(", ")]))
								floss_score_cds = calculate_floss_score(reference_floss_distribution, floss_distribution_cds)
							try:
								if float(line[fields["ribo_fpkm"]]) >= float(args.min):
									try:
										translated_floss_scores.append(floss_score_cds)
									except ValueError:
										pass
								else:
									try:
										not_translated_floss_scores.append(floss_score_cds)
									except ValueError:
										pass
							except ValueError:
								pass
						print "\t".join([str(field) for field in [id_count, line[fields["chr"]], line[fields["strand"]], line[fields["gene_id"]], line[fields["transcript_id"]], line[fields["gene_type"]], line[fields["ribo_fpkm"]], line[fields["coordinates_5UTR"]], line[fields["coordinates_CDS"]], line[fields["coordinates_3UTR"]], line[fields["SPEC_score_UTR5"]], spec_score_posterior_utr5, line[fields["SPEC_signal_UTR5"]], spec_signal_posterior_utr5, line[fields["FULL_score_UTR5"]], full_score_posterior_utr5, line[fields["FULL_signal_UTR5"]], floss_score_utr5, line[fields["FLOSS_distribution_UTR5"]], line[fields["ORF_score_UTR5"]], line[fields["ORF_reads_UTR5"]], line[fields["SPEC_score_CDS"]], spec_score_posterior_cds, line[fields["SPEC_signal_CDS"]], spec_signal_posterior_cds, line[fields["FULL_score_CDS"]], full_score_posterior_cds, line[fields["FULL_signal_CDS"]], floss_score_cds, line[fields["FLOSS_distribution_CDS"]], line[fields["ORF_score_CDS"]], line[fields["ORF_reads_CDS"]], line[fields["SPEC_score_UTR3"]], spec_score_posterior_utr3, line[fields["SPEC_signal_UTR3"]], spec_signal_posterior_utr3, line[fields["FULL_score_UTR3"]], full_score_posterior_utr3, line[fields["FULL_signal_UTR3"]], floss_score_utr3, line[fields["FLOSS_distribution_UTR3"]], line[fields["ORF_score_UTR3"]], line[fields["ORF_reads_UTR3"]]]])
						id_count += 1
	logger.info("__main__: Calculating transcript-level metrics... [COMPLETE]")

	translation_threshold = calculate_translation_threshold(translated_spec_scores, not_translated_spec_scores, float(args.fdr), int(args.spec))
	orfscore_auc = calculate_orfscore_auc(translated_orf_scores, not_translated_orf_scores, int(args.orfscore))
	spectre_auc = calculate_spectre_auc(translated_spec_scores, not_translated_spec_scores, int(args.spec))
	floss_auc = calculate_floss_auc(translated_floss_scores, not_translated_floss_scores, int(args.floss))
	full_auc = calculate_full_auc(translated_full_scores, not_translated_full_scores, int(args.full))

	logger.info("#########################################################")
	logger.info("Experimental metrics for combined split SPECtre analyses:")
	logger.info("Translational Threshold = " + str(translation_threshold))
	logger.info("SPECtre AUC = " + str(spectre_auc))
	logger.info("ORFscore AUC = " + str(orfscore_auc))
	logger.info("FLOSS AUC = " + str(floss_auc))
	logger.info("Full AUC = " + str(full_auc))
	logger.info("Reference Distribution = " + str(reference_floss_distribution))


















