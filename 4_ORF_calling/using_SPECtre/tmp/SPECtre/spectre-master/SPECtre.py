#!/usr/bin/python

# DESCRIPTION: Takes as input a BAM alignment file of ribosome profiling reads,
# gene or transcript abundance estimates from Cufflinks (mRNA-Seq or ribo-Seq)
# and a GTF annotation file (UCSC genes.gtf or Ensembl have been tested. From
# these inputs, the normalized read coverage over each transcript is computed
# and the spectral coherence relative to an enriched signal over a single
# reading frame is calculated (eg. 1:0:0, repeating). Additional implementations
# to calculate the FLOSS read distribution and FLOSS score (Ingolia, 2014), and
# ORFScore (Bazzini, 2014) are included. This version is able to handle one
# sample (BAM) at a time, however an additional script (intersect_experiments.py
# is included to merge multiple experiments into a single analysis.

# If you use SPECtre as part of your analysis, please cite the following:
#

# This script was written and tested using Python v2.7.8 on an OS X v10.9.4
# workstation and on a server running RHEL r6.4.

# DEPENDENCIES:
# HTSeq	(v0.6.1+):		https://pypi.python.org/pypi/HTSeq
# NumPy	(v1.11.0+):		https://pypi.python.org/pypi/numpy		
# pysam (v0.9.1+):		https://pypi.python.org/pypi/pysam/
# pyfasta (v0.5.2+):	https://pypi.python.org/pypi/pyfasta/
# rpy:					http://rpy.sourceforge.net/
# ROCR (in R)			https://rocr.bioinf.mpi-sb.mpg.de/

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

def add_custom_offsets(custom, default):
	custom_offsets = custom.split(",")
	for offset in custom_offsets:
		read, pos = offset.split(":")
		if "custom" in default:
			if int(read) in default["custom"]:
				pass
			else:
				default["custom"][int(read)] = int(pos)
		else:
			default["custom"] = {int(read): int(pos)}
	return default

def check_chromosomes(bam_file, gtf_file, cufflinks_file, target_chroms):
	# Checks the format of the chromosomes in each input file for compatibility.
	def get_bam_chroms(bam_file):
		header = os.popen("samtools view -H " + bam_file)
		chroms = [re.findall("SN:\w{1,5}", line)[0].split(":")[-1] for line in header if "@SQ" in line]
		return sorted(chroms)

	def get_gtf_chroms(gtf_file):
		chroms = list()
		for line in open(gtf_file):
			if not line.startswith("#"):
				fields = line.strip().split("\t")
				if fields[0] not in chroms:
					chroms.append(fields[0])
		return sorted(chroms)

	def get_cufflinks_chroms(cufflinks_file):
		pos = 0
		header, lines = open(cufflinks_file).readlines()[0], open(cufflinks_file).readlines()[1:]
		for i, j in enumerate(header.split("\t")):
			if j.upper() == "LOCUS":
				pos = i

		chroms = list()
		for line in lines:
			locus = line.strip().split("\t")[pos]
			chrom = locus.split(":")[0]
			if chrom not in chroms:
				chroms.append(chrom)
		return sorted(chroms)

	def get_target_chroms(targets):
		if targets:
			return [str(chrom) for chrom in targets.split(",")]
		else:
			return ""

	if get_target_chroms(target_chroms):
		return True if all(chrom in get_bam_chroms(bam_file) for chrom in get_gtf_chroms(gtf_file) for chrom in get_cufflinks_chroms(cufflinks_file) for chrom in get_target_chroms(target_chroms)) else False
	else:
		return True if all(chrom in get_bam_chroms(bam_file) for chrom in get_gtf_chroms(gtf_file) for chrom in get_cufflinks_chroms(cufflinks_file)) else False

def convert_cigar_to_reference_coordinates(cigar):
	coordinates = list()
	for op in cigar:
		if not op.type == "N":
			coordinates.extend(range(op.ref_iv.start, op.ref_iv.end))
	return sorted(set(coordinates))

def extract_coverage_over_interval(coverage, interval):
	interval_coverage = list()
	try:
		for i, v in coverage[interval].steps():
			if isinstance(v, set):
				interval_coverage.extend([len(v)]*len([pos for pos in i.xrange()]))
			else:
				interval_coverage.extend([v]*len([pos for pos in i.xrange()]))
	except TypeError:
		return "NA"
	return interval_coverage

def calculate_normalized_coverage(coverage):
	peak = float(max(coverage))
	return [cov/peak for cov in coverage]

def format_feature(feature):
	if feature == "exon":
		return "CDS"
	else:
		if "UTR" in feature:
			return "UTR"
	return feature

#############################
# SPECTRE SCORING FUNCTIONS #
#############################
def mean(signal):
	return math.fsum(signal) / float(len(signal))

def median(signal):
	sorted_signal = sorted(signal)
	midpoint = (len(sorted_signal)-1) // 2
	return sorted_signal[midpoint] if len(sorted_signal) % 2 else math.fsum([sorted_signal[midpoint], sorted_signal[midpoint+1]]) / 2.0

def maximum(signal):
	return sorted(signal)[-1]

def nonzero_mean(signal):
	nonzero_scores = [n for n in signal if n > 0]
	return math.fsum(nonzero_scores) / float(len(nonzero_scores))

def nonzero_median(signal):
	nonzero_scores = sorted([n for n in signal if n > 0])
	midpoint = (len(nonzero_scores)-1) // 2
	return nonzero_scores[midpoint] if len(nonzero_scores) % 2 else math.fsum([nonzero_scores[midpoint], nonzero_scores[midpoint+1]]) / 2.0

#########################
# TRANSCRIPT VALIDATION #
#########################
class Checks(object):
	# A series of QC checks for minimum length, coverage, abundance, etc. The format of the region
	# variable must be a list with raw read coverage by position. For A-site coverage, the set()
	# contents must be converted to its length.
	def __init__(self, window_length, buffers, strand, region, transcript, fpkms):
		self.window_length = window_length
		self.buffers = buffers
		self.strand = strand
		self.region = region
		self.transcript = transcript
		self.fpkms = fpkms

	@staticmethod
	def trim(sector, strand, left_trim, right_trim):
		if strand == "+":
			if left_trim == 0:
				return sector if right_trim == 0 else sector[:-right_trim]
			else:
				return sector[left_trim:] if right_trim == 0 else sector[left_trim:-right_trim]
		else:
			if left_trim == 0:
				return sector if right_trim == 0 else sector[right_trim:]
			else:
				return sector[:-left_trim] if right_trim == 0 else sector[right_trim:-left_trim]

	def coverage(self):
		# Checks for sufficient coverage over the region of interest:
		if self.region in (None, "NA", ""):
			return False
		else:
			max_left_trim, max_right_trim = self.buffers
			trimmed_region = self.trim(self.region, self.strand, max_left_trim, max_right_trim)
			if sum(trimmed_region) > 0:
				# Check for minimal coverage:
				return True
			else:
				return False

	def representation(self):
		# Checks for sufficient distribution of coverage over the region of interest:
		if self.region in (None, "NA", ""):
			return False
		else:
			max_left_trim, max_right_trim = self.buffers
			trimmed_region = self.trim(self.region, self.strand, max_left_trim, max_right_trim)
			try:
				ratio = float(len([n for n in trimmed_region if n > 0]) / len(trimmed_region))
				if ratio >= 0.0:
					# Check for distribution of coverage, increase for higher stringency:
					return True
				else:
					return False
			except ZeroDivisionError:
				return False

	def distribution(self):
		# Checks if a single position if responsible for most of the coverage over the region:
		if self.region in (None, "NA", ""):
			return False
		else:
			max_left_trim, max_right_trim = self.buffers
			trimmed_region = self.trim(self.region, self.strand, max_left_trim, max_right_trim)
			try:
				ratio = float(max([n for n in trimmed_region if n > 0]) / sum([n for n in trimmed_region if n > 0]))
				if ratio <= 0.5:
					# Check for single position over-representation, increase for higher stringency:
					return True
				else:
					return False
			except ValueError:
				return False

	def length(self):
		# Checks that the region of interest is of sufficient length:
		if self.region in (None, "NA", ""):
			return False
		else:
			max_left_trim, max_right_trim = self.buffers
			trimmed_region = self.trim(self.region, self.strand, max_left_trim, max_right_trim)
			if len(trimmed_region) >= self.window_length:
				# The region must be able to accomodate at least one window for SPECtre, otherwise
				# the Full Coherence may be optionally calculated. But for SPECtre, a region of
				# insufficient length to fit a single window is failed:
				return True
			else:
				return False

	def abundance(self):
		# Checks for minimal abundance (from Cufflinks):
		if self.transcript in self.fpkms:
			if self.fpkms[self.transcript] > 0:
				return True
			else:
				return False
		else:
			return False

##############
# ANNOTATION #
##############
def extract_fpkms(cufflinks_file, targets):
	# This function takes as input a Cufflinks isoforms.fpkm_tracking output
	# file and extracts the transcript-level expression in FPKM. Pleast note 
	# that, as written, this function should support a generic tab-delimited
	# file such that the name of the transcript is in the first column, and the
	# column holding the expression values is titled as "FPKM".
	logger.info("extract_fpkms(): Parsing transcript FPKMs from file: " + cufflinks_file + " to memory... [STARTED]")
	fpkms = dict()
	position = 0
	for linenum, line in enumerate(open(cufflinks_file)):
		if linenum == 0:
			for i, j in enumerate(line.strip().split("\t")):
				if j.upper() == "FPKM":
					position = i
		else:
			tracking_id, gene_id = line.strip().split("\t")[0], line.strip().split("\t")[3]
			if not tracking_id == gene_id:
				chrom, fpkm = line.strip().split("\t")[6].split(":")[0], float(line.strip().split("\t")[position])
				if len(targets) == 0 or chrom in targets:
					if tracking_id not in fpkms:
						fpkms[tracking_id] = fpkm
	logger.info("extract_fpkms(): Parsing transcript FPKMs from file: " + cufflinks_file + " to memory... [COMPLETE]")
	return fpkms

##############
# ANNOTATION #
##############
def parse_gtf(gtf_file, fpkms, targets):
	# This function takes as input a user-supplied GTF transcript annotation file
	# and extracts the CDS and UTR intervals of protein-coding genes, and the start
	# and end coordinates of non-coding transcripts. The coordinates are loaded
	# into HTSeq GenomicInteval() objects for more efficient access from memory.
	def partition_utrs(gtf):
		# Since Ensembl does not differentiate between 5' and 3' UTRs, they must be
		# annotated based on the position relative to the CDS start (or stop) and
		# strand of the transcript model:
		def find_cds_start(ivs):
			coordinates = list()
			strand = str()
			for iv in ivs:
				coordinates.extend([p.pos for p in iv.xrange()])
				strand = iv.strand
			return int(sorted(coordinates)[0]) if strand == "+" else int(sorted(coordinates)[-1])

		logger.info("parse_gtf/partition_utr_coordinates(): Partitioning 5' and 3' UTR entries from loaded GTF.. [STARTED]")
		flat_gtf = [(transcript, intervals) for transcript, intervals in flatten(gtf).items() if "protein_coding" in transcript if "UTR" in transcript]
		for transcript, intervals in flat_gtf:
			gene_type, chrom, strand, gene_id, transcript_id, feature = transcript.split(":")
			five_prime, three_prime = list(), list()
			try:
				cds_start = find_cds_start(gtf[gene_type][chrom][strand][gene_id][transcript_id]["CDS"])
				if strand == "+":
					[five_prime.append(iv) if iv.end <= cds_start else three_prime.append(iv) for iv in intervals]
				else:
					[five_prime.append(iv) if iv.start >= cds_start else three_prime.append(iv) for iv in intervals]
				# Remove the previous combined entry for UTRs of a transcript:
				del gtf[gene_type][chrom][strand][gene_id][transcript_id]["UTR"]
				# Replace with separate key entries for 5'UTR and 3'UTR:
				gtf[gene_type][chrom][strand][gene_id][transcript_id]["UTR5"] = five_prime
				gtf[gene_type][chrom][strand][gene_id][transcript_id]["UTR3"] = three_prime
			except KeyError:
				pass
		logger.info("parse_gtf/partition_utr_coordinates(): Partitioning 5' and 3' UTR entries from loaded GTF.. [COMPLETE]")
		return gtf

	logger.info("parse_gtf(): Parsing transcript coordinates from GTF: " + gtf_file + " into memory... [STARTED]")
	# Load the GTF file into memory:
	gtf = HTSeq.GFF_Reader(gtf_file)
	intervals = HTSeq.GenomicArrayOfSets(chroms="auto", stranded=True)
	transcripts = hash()
	try:
		for record in gtf:
			try:
				# Extract the necessary fields from the GTF record:
				chrom, strand, gene_id, transcript_id = record.iv.chrom, record.iv.strand, record.attr["gene_id"], record.attr["transcript_id"]
				if chrom in targets or len(targets) == 0:
					if transcript_id in fpkms:
						if fpkms[transcript_id] > 0.0:
							biotype_field = "gene_biotype" if "gene_biotype" in record.attr.keys() else "gene_type" if "gene_type" in record.attr.keys() else str()
							if record.attr[biotype_field] == "protein_coding":
								if record.type in ("CDS", "stop_codon"):
									if isinstance(transcripts[record.attr[biotype_field]][chrom][strand][gene_id][transcript_id]["CDS"], list):
										transcripts[record.attr[biotype_field]][chrom][strand][gene_id][transcript_id]["CDS"].append(record.iv)
									else:
										transcripts[record.attr[biotype_field]][chrom][strand][gene_id][transcript_id]["CDS"] = [record.iv]
									intervals[record.iv] += ":".join([record.attr[biotype_field], chrom, strand, gene_id, transcript_id, "CDS"])
								elif record.type == "UTR":
									if isinstance(transcripts[record.attr[biotype_field]][chrom][strand][gene_id][transcript_id]["UTR"], list):
										transcripts[record.attr[biotype_field]][chrom][strand][gene_id][transcript_id]["UTR"].append(record.iv)
									else:
										transcripts[record.attr[biotype_field]][chrom][strand][gene_id][transcript_id]["UTR"] = [record.iv]
									intervals[record.iv] += ":".join([record.attr[biotype_field], chrom, strand, gene_id, transcript_id, "UTR"])
							else:
								if record.type == "exon":
									if isinstance(transcripts[record.attr[biotype_field]][chrom][strand][gene_id][transcript_id]["exon"], list):
										transcripts[record.attr[biotype_field]][chrom][strand][gene_id][transcript_id]["exon"].append(record.iv)
									else:
										transcripts[record.attr[biotype_field]][chrom][strand][gene_id][transcript_id]["exon"] = [record.iv]
									intervals[record.iv] += ":".join([record.attr[biotype_field], chrom, strand, gene_id, transcript_id, "exon"])
			except KeyError:
				pass
	except ValueError:
		print "parse_gtf():"
		print record.iv
		print record.attr
	logger.info("parse_gtf(): Parsing transcript coordinates from GTF: " + gtf_file + " into memory... [COMPLETE]")
	# Return the parsed GFF_Reader() object, the GenomicIntervals() array, and the parsed GTF:
	return gtf, intervals, partition_utrs(transcripts)

##########################################
# READ ALIGNMENT AND COVERAGE ADJUSTMENT #
##########################################
class Coverage(object):
	# Take in as input a BAM file of read alignments, then return them as an
	# HTSeq BAM_Reader object. In addition, convert the positions BAM read
	# alignments to their A/P-site position as a GenomicArray (coverage) object. 
	def __init__(self, asite_offsets, psite_offsets, targets, bam):
		self.asite_offsets = asite_offsets
		self.psite_offsets = psite_offsets
		self.targets = targets
		self.bam = bam

	@staticmethod
	def calculate_offset_position(offsets, read_length):
		if read_length in offsets:
			return offsets[read_length]
		else:
			return read_length / 2

	def psite_coverage(self):
		# This method loads alignments from the input BAM file into a GenomicArray
		# object. This results in a hashed coverage object with the raw read coverage
		# converted to step() values.
		coverage = HTSeq.GenomicArray(chroms="auto", stranded=True, typecode="i", storage="step")
		for read in self.bam:
			# Extract the read chromosome and strand from the GenomicInterval() of the read:
			read_chrom = read.iv.chrom
			read_strand = read.iv.strand
			if read_chrom in self.targets or len(self.targets) == 0:
				# Based on the read CIGAR string, convert the aligned read position to its P-site
				# offset position:
				offset = self.calculate_offset_position(self.psite_offsets, len(read.get_sam_line().split("\t")[9]))
				if read_strand == "+":
					# Since Python indexes are 0-based:
					reference_pos = convert_cigar_to_reference_coordinates(read.cigar)[offset - 1]
				else:
					reference_pos = convert_cigar_to_reference_coordinates(read.cigar)[-offset]
				# Add the P-site offset position as a GenomicPosition in the coverage array:
				coverage[HTSeq.GenomicPosition(read_chrom, reference_pos, read_strand)] += 1
		return coverage

	def asite_coverage(self):
		# This method loads alignments from the input BAM file into a GenomicArrayOfSets
		# object. This results in a hashed coverage object of GenomicIntervals with the
		# read length and number stored as tuples in the set() array.
		coverage = HTSeq.GenomicArrayOfSets(chroms="auto", stranded=True)
		try:
			for read in self.bam:
				# Extract the read chromosome and strand from the GenomicInterval of the read:
				read_chrom = read.iv.chrom
				read_strand = read.iv.strand
				if read_chrom in self.targets or len(self.targets) == 0:
					# Based on the CIGAR string of the read, convert the aligned read position to its
					# A-Site offset position:
					name = read.get_sam_line().split("\t")[0] + "|" + str(len(read.get_sam_line().split("\t")[9]))
					offset = self.calculate_offset_position(self.asite_offsets, len(read.get_sam_line().split("\t")[9]))
					if read_strand == "+":
						reference_start = convert_cigar_to_reference_coordinates(read.cigar)[offset - 1]
						reference_end = convert_cigar_to_reference_coordinates(read.cigar)[offset]
					else:
						reference_start = convert_cigar_to_reference_coordinates(read.cigar)[-offset - 1]
						reference_end = convert_cigar_to_reference_coordinates(read.cigar)[-offset]
					# Add the A-site offset position as a 1-base length GenomicInterval in the coverage array:
					coverage[HTSeq.GenomicInterval(read_chrom, reference_start, reference_end, read_strand)] += name
		except AttributeError:
			print "Coverage/asite_coverage():"
			print read
			print read.iv
		return coverage

########################################
# TRANSCRIPT AND ORF SCORING FUNCTIONS #
########################################
class Coherence(object):
	# Take as input a transcript or ORF formatted as GenomicIntervals() and the
	# A-site adjusted coverage. Normalize the transcript or ORF coverage to the
	# highest position(s), then calculate the windowed SPECtre score, the SPECtre
	# signal, and the full length coherence score. The full transcript or ORF
	# signal is output as the normalized read coverage.
	def __init__(self, psite_coverage, psite_buffers, window_length, spectre_type, step_size, methods, fpkms, transcript_intervals):
		self.psite_coverage = psite_coverage
		self.psite_buffers = psite_buffers
		self.window_length = window_length
		self.spectre_type = spectre_type
		self.step_size = step_size
		self.methods = methods
		self.fpkms = fpkms
		self.transcript_intervals = transcript_intervals

	def annotation(self):
		transcript, intervals = self.transcript_intervals
		return transcript

	def p_coverage(self):
		# Calculate the coherence over each window in a transcript or region:
		transcript, intervals = self.transcript_intervals
		gene_type, chrom, strand, gene_id, transcript_id, feature = transcript.split(":")
		# Extract the raw read coverage from the P-site adjusted coverage input and the
		# region interval(s):
		raw_coverage = list()
		for interval in intervals:
			raw_coverage.extend(extract_coverage_over_interval(self.psite_coverage, interval))
		return raw_coverage

	def spectre_signal(self):
		# Calculate the coherence over each window in a transcript or region:
		transcript, intervals = self.transcript_intervals
		gene_type, chrom, strand, gene_id, transcript_id, feature = transcript.split(":")
		# Extract the raw read coverage from the P-site adjusted coverage input and the
		# region interval(s):
		raw_coverage = list()
		for interval in intervals:
			raw_coverage.extend(extract_coverage_over_interval(self.psite_coverage, interval))
		# Submit the raw coverage for QC:
		check = Checks(self.window_length, self.psite_buffers[format_feature(feature)], strand, raw_coverage, transcript_id, self.fpkms)
		if check.coverage() == True and check.representation() == True and check.distribution() == True and check.length() == True and check.abundance() == True:
			# Convert the raw coverage to the normalized coverage based on the highest covered position:
			normalized_coverage = calculate_normalized_coverage(raw_coverage)
			# Generate the reference coverage signal based on the length of the coverage region:
			reference_signal = ([4/6.0,1/6.0,1/6.0]*int(math.ceil(len(normalized_coverage)/3.0)))[0:len(normalized_coverage)]
			# Initialize the calculated signal variable:
			coherences = list()
			for i in range(0, len(normalized_coverage))[::int(self.step_size)]:
				j = i + self.window_length
				if (math.fsum(normalized_coverage[i:j]) == 0) or (len(normalized_coverage[i:j]) < self.window_length):
					coherences.append(0.0)
				else:
					r('window.region <- c(%s)' % ",".join(str(n) for n in normalized_coverage[i:j]))
					r('window.coding <- c(%s)' % ",".join(str(n) for n in reference_signal[i:j]))
					r('test.spec <- spec.pgram(data.frame(window.region, window.coding), spans=c(3,3), plot=FALSE)')
					coherences.append(r('test.spec$coh[which(abs(test.spec$freq-1/3)==min(abs(test.spec$freq-1/3)))]')[0])
			return coherences
		else:
			return "NA"

	def spectre_score(self):
		transcript, intervals = self.transcript_intervals
		gene_type, chrom, strand, gene_id, transcript_id, feature = transcript.split(":")
		# Extract the raw read coverage from the P-site adjusted coverage input and the
		# region interval(s):
		raw_coverage = list()
		for interval in intervals:
			raw_coverage.extend(extract_coverage_over_interval(self.psite_coverage, interval))
		# Submit the raw coverage for QC:
		check = Checks(self.window_length, self.psite_buffers[format_feature(feature)], strand, raw_coverage, transcript_id, self.fpkms)
		if check.coverage() == True and check.representation() == True and check.distribution() == True and check.length() == True and check.abundance() == True:
			if self.spectre_signal() == "NA":
				# Cannot calculate a SPECtre score:
				return "NA"
			else:
				if self.spectre_type == "mean":
					return mean(self.spectre_signal())
				elif self.spectre_type == "median":
					return median(self.spectre_signal())
				elif self.spectre_type == "maximum":
					return maximum(self.spectre_signal())
				elif self.spectre_type == "nonzero_mean":
					return nonzero_mean(self.spectre_signal())
				elif self.spectre_type == "nonzero_median":
					return nonzero_median(self.spectre_signal())
				else:
					# Default to calculate the mean of the spectre_signal():
					return mean(self.spectre_signal())
		else:
			return "NA"

	def spectre_phase(self):
		transcript, intervals = self.transcript_intervals
		gene_type, chrom, strand, gene_id, transcript_id, feature = transcript.split(":")
		# Extract the raw read coverage from the P-site adjusted coverage input and the
		# region interval(s):
		raw_coverage = list()
		for interval in intervals:
			raw_coverage.extend(extract_coverage_over_interval(self.psite_coverage, interval))
		# Submit the raw coverage for QC:
		check = Checks(self.window_length, self.psite_buffers[format_feature(feature)], strand, raw_coverage, transcript_id, self.fpkms)
		if check.coverage() == True and check.representation() == True and check.distribution() == True and check.length() == True and check.abundance() == True:
			if self.spectre_signal() in (None, "NA", ""):
				return "NA"
			else:
				# Convert the raw coverage to the normalized coverage based on the highest covered position:
				normalized_coverage = calculate_normalized_coverage(raw_coverage)
				# Generate the reference coverage signal based on the length of the normalized region:
				reference_signal = ([4/6.0,1/6.0,1/6.0]*int(math.ceil(len(normalized_coverage)/3.0)))[0:len(normalized_coverage)]
				r('test.region <- c(%s)' % ",".join(str(n) for n in normalized_coverage))
				r('test.coding <- c(%s)' % ",".join(str(n) for n in reference_signal))
				r('spec.coding <- spec.pgram(data.frame(test.region, test.coding), spans=c(3,3), plot=FALSE)')
				frame = r('round(test.spec$phase[which(abs(test.spec$freq-1/3)==min(abs(test.spec$freq-1/3)))])')[0]
				return "+0" if frame == 0 else "+1" if frame == 2 else "+2" if frame == -2 else "NA"
		else:
			return "NA"

	def coherence_signal(self):
		transcript, intervals = self.transcript_intervals
		gene_type, chrom, strand, gene_id, transcript_id, feature = transcript.split(":")
		# Extract the raw read coverage from the P-site adjusted coverage input and the
		# region interval(s):
		raw_coverage = list()
		for interval in intervals:
			raw_coverage.extend(extract_coverage_over_interval(self.psite_coverage, interval))
		# Submit the raw coverage for QC:
		check = Checks(self.window_length, self.psite_buffers[format_feature(feature)], strand, raw_coverage, transcript_id, self.fpkms)
		if check.coverage() == True and check.representation() == True and check.distribution() == True and check.length() == True and check.abundance() == True:
			normalized_coverage = calculate_normalized_coverage(raw_coverage)
			return normalized_coverage if "Full" in self.methods else "NA"
		else:
			return "NA"

	def coherence_score(self):
		transcript, intervals = self.transcript_intervals
		gene_type, chrom, strand, gene_id, transcript_id, feature = transcript.split(":")
		# Extract the raw read coverage from the P-site adjusted coverage input and the
		# region interval(s):
		raw_coverage = list()
		for interval in intervals:
			raw_coverage.extend(extract_coverage_over_interval(self.psite_coverage, interval))
		# Submit the raw coverage for QC:
		check = Checks(self.window_length, self.psite_buffers[format_feature(feature)], strand, raw_coverage, transcript_id, self.fpkms)
		if check.coverage() == True and check.representation() == True and check.distribution() == True and check.length() == True and check.abundance() == True:
		# As long as there is sufficient read coverage, the full coherence should be calculated for transcripts
		# of all lengths:
			if self.coherence_signal() in (None, "NA", ""):
				return "NA"
			else:
				if "Full" in self.methods:
					normalized_coverage = calculate_normalized_coverage(raw_coverage)
					# Load the normalized and reference signals into R and calculate the spectral coherence
					# over the full length of the normalized region:
					reference_signal = ([4/6.0, 1/6.0, 1/6.0]*int(math.ceil(len(normalized_coverage)/3.0)))[0:len(normalized_coverage)]
					r('test.region <- c(%s)' %",".join(str(n) for n in normalized_coverage))
					r('test.coding <- c(%s)' %",".join(str(n) for n in reference_signal))
					r('spec.coding <- spec.pgram(data.frame(test.region, test.coding), spans=c(3,3), plot=FALSE)')
					return r('spec.coding$coh[which(abs(spec.coding$freq-1/3)==min(abs(spec.coding$freq-1/3)))]')[0]
				else:
					return "NA"
		else:
			return "NA"

class FLOSS(object):
	# Calculate the FLOSS distribution for the input transcript or ORF based on the intervals provided:
	def __init__(self, fpkms, methods, window_length, asite_buffers, asite_coverage, transcript_intervals):
		self.fpkms = fpkms
		self.methods = methods
		self.window_length = window_length
		self.asite_buffers = asite_buffers
		self.asite_coverage = asite_coverage
		self.transcript_intervals = transcript_intervals

	def annotation(self):
		transcript, intervals = self.transcript_intervals
		return transcript

	def a_coverage(self):
		transcript, intervals = self.transcript_intervals
		gene_type, chrom, strand, gene_id, transcript_id, feature = transcript.split(":")
		# Extract the raw read coverage from the A-site adjusted coverage input and the
		# region interval(s):
		raw_coverage = list()
		for interval in intervals:
			raw_coverage.extend(extract_coverage_over_interval(self.asite_coverage, interval))
		return raw_coverage

	def distribution(self):
		def calculate_read_distribution(asite_coverage, ivs):
			distribution = dict(zip(range(24,37), [0]*len(range(24,37))))
			try:
				for interval in ivs:
					for iv, reads in list(asite_coverage[interval].steps()):
						if len(reads) > 0:
							for read in reads:
								read_name, read_length = read.split("|")
								if int(read_length) <= 24:
									distribution[24] += 1
								elif int(read_length) >= 36:
									distribution[36] += 1
								else:
									distribution[int(read_length)] += 1
				return {k: v * (1/float(sum(distribution.values()))) for k, v in distribution.items()}
			except ValueError:
				return "NA"

		if "FLOSS" in self.methods:
			transcript, intervals = self.transcript_intervals
			gene_type, chrom, strand, gene_id, transcript_id, feature = transcript.split(":")
			# Extract the raw coverage over the region:
			raw_coverage = self.a_coverage()
			# Submit the raw coverage for QC:
			check = Checks(self.window_length, self.asite_buffers[format_feature(feature)], strand, raw_coverage, transcript_id, self.fpkms)
			if check.coverage() == True and check.length() == True and check.abundance() == True:
				read_distribution = calculate_read_distribution(self.asite_coverage, intervals)
				return read_distribution
			else:
				return "NA"
		else:
			return "NA"

class ORF(object):
	# Calculate the ORFscore based on the enrichment of reads in the first frame of a transcript or ORF relative to the
	# number of reads in the other two frames:
	def __init__(self, fpkms, methods, orf_buffers, window_length, psite_coverage, transcript_intervals ):
		self.fpkms = fpkms
		self.methods = methods
		self.orf_buffers = orf_buffers
		self.window_length = window_length
		self.psite_coverage = psite_coverage
		self.transcript_intervals = transcript_intervals

	def annotation(self):
		transcript, intervals = self.transcript_intervals
		return transcript

	def reads(self):
		frame_coverage = list()
		transcript, intervals = self.transcript_intervals
		gene_type, chrom, strand, gene_id, transcript_id, feature = transcript.split(":")
		if "ORFscore" in self.methods:
			# Extract the raw coverage over the region:
			raw_coverage = list()
			for interval in intervals:
				raw_coverage.extend(extract_coverage_over_interval(self.psite_coverage, interval))
			# Submit the raw coverage for QC:
			check = Checks(self.window_length, (self.orf_buffers[format_feature(feature)]), strand, raw_coverage, transcript_id, self.fpkms)
			if check.coverage() == True and check.representation() == True and check.distribution() == True and check.abundance() == True:
				if strand == "+":
					frame_coverage = raw_coverage[self.orf_buffers[format_feature(feature)][0]:] if self.orf_buffers[format_feature(feature)][-1] == 0 else raw_coverage[self.orf_buffers[format_feature(feature)][0]:-self.orf_buffers[format_feature(feature)][-1]]
				else:
					frame_coverage = raw_coverage[self.orf_buffers[format_feature(feature)][-1]:] if self.orf_buffers[format_feature(feature)][0] == 0 else raw_coverage[self.orf_buffers[format_feature(feature)][-1]:-self.orf_buffers[format_feature(feature)][0]]
				if sum(frame_coverage) == 0:
					return "NA"
				else:
					masked_region = list()
					for x in xrange(len(frame_coverage)):
						masked_region.append(0) if frame_coverage[x]/math.fsum(frame_coverage) >= 0.7 else masked_region.append(frame_coverage[x])
					return [sum(masked_region[i::3]) for i in (0,1,2)] if strand == "+" else [sum(masked_region[::-1][i::3]) for i in (0,1,2)]	
			else:
				return "NA"
		else:
			return "NA"

	def score(self):
		def frame_score(reads, mean_reads):
			return math.pow((reads - mean_reads), 2) / mean_reads

		transcript, intervals = self.transcript_intervals
		gene_type, chrom, strand, gene_id, transcript_id, feature = transcript.split(":")
		if "ORFscore" in self.methods:
			psite_reads = self.reads()
			if psite_reads == "NA":
				return "NA"
			else:
				frames_mean = sum(psite_reads) / float(len(psite_reads))
				score = math.log(math.fsum([frame_score(psite_reads[0], frames_mean), frame_score(psite_reads[1], frames_mean), frame_score(psite_reads[2], frames_mean), 1]), 2)
				return score if (psite_reads[0] > psite_reads[1]) and (psite_reads[0] > psite_reads[2]) else -score
		else:
			return "NA"

#################################################
# TRANSCRIPT SCORE CALCULATIONS AND AGGREGATION #
#################################################
def calculate_transcript_scores(gtf, fpkms, fpkm_cutoff, asite_buffers, psite_buffers, orfscore_buffers, bam_file, window_length, step_size, spectre_type, methods, offsets, target_chroms, threads):

	'''
	For each transcript, this function is to calculate (if so designated) its SPECtre metrics
	(score, and windowed coherence signal), FLOSS metrics (score, and read distribution), and
	ORFscore metrics (score, and distribution of reads over each frame). If directed, the
	SPECtre metrics over a given window length(s) will also be calculated.
	'''

	def parse_offsets(offsets):
		if "custom" in offsets:
			return offsets["ingolia"], offsets["custom"]
		else:
			return offsets["ingolia"], offsets["bazzini"]

	def calculate_reference_distribution(distributions):
		reference_distribution = {26: 0, 27: 0, 28: 0, 29: 0, 30: 0, 31: 0, 32: 0, 33: 0, 34: 0}
		reference_transcripts = 0
		for dist in distributions:
			if dist == "NA":
				pass
			elif sum(dist.values()) > 0:
				for read_length in dist:
					if read_length in reference_distribution:
						reference_distribution[read_length] += dist[read_length]
				reference_transcripts += 1
			else:
				pass
		return {k: v * (1/float(reference_transcripts)) for k, v in reference_distribution.items()}

	def calculate_floss_score(reference_distribution, transcript_score):
		annotation, distribution = transcript_score
		gene_type, chrom, strand, gene_id, transcript_id, feature = annotation.split(":")
		if isinstance(distribution, dict):
			if sum(distribution.values()) > 0:
				floss = 0
				for read_length in reference_distribution:
					if read_length in distribution:
						floss += abs(distribution[read_length] - reference_distribution[read_length])
					else:
						floss += abs(0.0 - float(reference_distribution[read_length]))
				return floss / 2.0
			else:
				return "NA"
		else:
			return "NA"

	def build_translation_distributions(fpkms, fpkm_cutoff, transcript_scores):
		translated, not_translated = list(), list()
		for annotation, score in transcript_scores:
			gene_type, chrom, strand, gene_id, transcript_id, feature = annotation.split(":")
			if transcript_id in fpkms and not score == "NA":
				translated.append(score) if float(fpkms[transcript_id]) >= float(fpkm_cutoff) else not_translated.append(score)
		return translated, not_translated

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
		annotation, score = transcript_score
		gene_type, chrom, strand, gene_id, transcript_id, feature = annotation.split(":")
		if isinstance(score, int) or isinstance(score, float):
			return annotation, calculate_posterior_probability(coding_scores, noncoding_scores, score)
		elif isinstance(score, list):
			windowed = list()
			for coh in score:
				windowed.append(calculate_posterior_probability(coding_scores, noncoding_scores, coh))
			return annotation, windowed
		else:
			return annotation, "NA"

	def populate_hash(transcripts, metric, transcript_metrics):
		for annotation, score in transcript_metrics:
			gene_type, chrom, strand, gene_id, transcript_id, feature = annotation.split(":")
			transcripts[gene_type][chrom][strand][gene_id][transcript_id][feature][metric] = score
		return transcripts

	logger.info("calculate_transcript_scores(): Calculating transcript-level metrics [STARTED].")
	# Instantiate the Pool objects for multiprocessing:
	pool = ThreadPool(threads)
	# Parse the transcript annotations and coordinates into separate lists:
	transcripts, intervals = zip(*flatten(gtf).iteritems())

	# Partition the transcript offsets:
	asite_offsets, psite_offsets = parse_offsets(offsets)

	########################################
	# EXTRACT COVERAGE FROM INPUT BAM FILE #
	######################################## 
	logger.info("calculate_transcript_scores(): Calculating A/P-site read distribution and coverage... [STARTED].")
	# Load reads from the input BAM file into an HTSeq Reader() object:
	logger.info("calculate_transcript_scores(): Loading alignments from: " + bam_file + " into BAM_Reader()... [STARTED].")
	bam = HTSeq.BAM_Reader(bam_file)
	logger.info("calculate_transcript_scores(): Loading alignments from: " + bam_file + " into BAM_Reader()... [COMPLETE].")
	# Extract the A/P-site adjusted coverage based on the input BAM file:
	array = Coverage(asite_offsets, psite_offsets, target_chroms, bam)
	# Partition coverages to their respective containers:
	asite_coverage = array.asite_coverage()
	psite_coverage = array.psite_coverage()
	logger.info("calculate_transcript_scores(): Calculating A/P-site read distribution and coverage... [COMPLETE].")

	###################################################################################################
	# CALCULATE THE SPECTRE, FULL COHERENCE, FLOSS READ DISTRIBUTION AND ORFSCORE FOR EACH TRANSCRIPT #
	###################################################################################################
	# Calculate the SPECtre and Coherence scores (if required) for each transcript:
	logger.info("calculate_transcript_scores(): Calculating SPECtre/Coherence scores... [STARTED].")
	Coherence_func = partial(Coherence, psite_coverage, psite_buffers, window_length, spectre_type, step_size, methods, fpkms)
	coherences = pool.map(Coherence_func, flatten(gtf).iteritems())
	# Partition coherence scores to their respective containers:
	transcript_coherence_signals = dict(zip([scores.annotation() for scores in coherences], [scores.coherence_signal() for scores in coherences]))
	transcript_coherence_scores = dict(zip([scores.annotation() for scores in coherences], [scores.coherence_score() for scores in coherences]))
	transcript_spectre_psite_coverages = dict(zip([scores.annotation() for scores in coherences], [scores.p_coverage() for scores in coherences]))
	transcript_spectre_signals = dict(zip([scores.annotation() for scores in coherences], [scores.spectre_signal() for scores in coherences]))
	transcript_spectre_scores = dict(zip([scores.annotation() for scores in coherences], [scores.spectre_score() for scores in coherences]))
	transcript_spectre_phases = dict(zip([scores.annotation() for scores in coherences], [scores.spectre_phase() for scores in coherences]))
	logger.info("calculate_transcript_scores(): Calculating SPECtre/Coherence scores... [COMPLETE].")

	if "FLOSS" in methods:
		# Calculate FLOSS read distributions for each transcript:
		logger.info("calculate_transcript_scores(): Calculating FLOSS read length distributions... [STARTED].")
		FLOSS_func = partial(FLOSS, fpkms, methods, window_length, asite_buffers, asite_coverage)
		floss_distributions = pool.map(FLOSS_func, flatten(gtf).iteritems())	
		# Partition FLOSS distributions to their respective containers:
		transcript_floss_asite_coverages = dict(zip([floss.annotation() for floss in floss_distributions], [floss.a_coverage() for floss in floss_distributions]))
		transcript_floss_distributions = dict(zip([floss.annotation() for floss in floss_distributions], [floss.distribution() for floss in floss_distributions]))
		logger.info("calculate_transcript_scores(): Calculating FLOSS read length distributions... [COMPLETE].")
		# Build FLOSS reference distribution:
		logger.info("calculate_transcript_scores(): Calculating FLOSS reference read distribution... [STARTED].")
		protein_coding_distributions = [distribution for transcript, distribution in transcript_floss_distributions.items() if "protein_coding" in transcript.lower() if "cds" in transcript.lower()]
		reference_distribution = calculate_reference_distribution(protein_coding_distributions)
		logger.info("calculate_transcript_scores(): Calculating FLOSS reference read distribution... [COMPLETE].")
		# Calculate FLOSS scores for each transcript using reference distribution:
		logger.info("calculate_transcript_scores(): Calculating FLOSS from read length distributions... [STARTED].")
		FLOSS_metric = partial(calculate_floss_score, reference_distribution)
		floss_scores = pool.map(FLOSS_metric, transcript_floss_distributions.iteritems())
		transcript_floss_scores = dict(zip(transcripts, floss_scores))
		logger.info("calculate_transcript_scores(): Calculating FLOSS from read length distributions... [COMPLETE].")

	if "ORFscore" in methods:
		# Calculate ORFscore for each transcript:
		logger.info("calculate_transcript_scores(): Calculating ORFscore from frame distributions... [STARTED].")
		ORF_func = partial(ORF, fpkms, methods, orfscore_buffers, window_length, psite_coverage)
		orf_scores = pool.map(ORF_func, flatten(gtf).iteritems())
		# Parition ORF scores to their respective containers:
		transcript_orf_scores = dict(zip([orf.annotation() for orf in orf_scores], [orf.score() for orf in orf_scores]))
		transcript_orf_reads = dict(zip([orf.annotation() for orf in orf_scores], [orf.reads() for orf in orf_scores]))
		logger.info("calculate_transcript_scores(): Calculating ORFscore from frame distributions... [COMPLETE].")

	#####################################################################################
	# BUILD DISTRIBUTIONS BASED ON TRANSLATIONAL STATUS FOR PROTEIN-CODING TRANSCRIPTS: #
	#####################################################################################
	if len(target_chroms) > 0:
		logger.info("calculate_transcript_scores(): Populating NAs for split analysis posteriors... [STARTED].")		
	else:
		logger.info("calculate_transcript_scores(): Building distributions for posteriors calculation... [STARTED].")
		# For transcripts with SPECtre scores:
		logger.info("calculate_transcript_scores(): Buliding distributions for SPECtre posteriors... [STARTED].")
		protein_coding_spectre_scores = [(transcript, score) for transcript, score in transcript_spectre_scores.items() if "protein_coding" in transcript.lower() if "cds" in transcript.lower()]
		translated_spectre_scores, untranslated_spectre_scores = build_translation_distributions(fpkms, fpkm_cutoff, protein_coding_spectre_scores)
		logger.info("calculate_transcript_scores(): Buliding distributions for SPECtre posteriors... [COMPLETE].")

		# For transcripts with Coherence scores:
		logger.info("calculate_transcript_scores(): Building distributions for Coherence posteriors... [STARTED].")
		protein_coding_coherence_scores = [(transcript, score) for transcript, score in transcript_coherence_scores.items() if "protein_coding" in transcript.lower() if "cds" in transcript.lower()]
		translated_coherence_scores, untranslated_coherence_scores = build_translation_distributions(fpkms, fpkm_cutoff, protein_coding_coherence_scores)
		logger.info("calculate_transcript_scores(): Building distributions for Coherence posteriors... [STARTED].")
		logger.info("calculate_transcript_scores(): Building distributions for posteriors calculation... [COMPLETE].")

	##############################################################################
	# CALCULATE THE SPECTRE AND COHERENCE SCORE POSTERIORS FOR EACH TRANSCRIPT : #
	##############################################################################
	if len(target_chroms) > 0:
		transcript_coherence_posteriors = dict(zip([transcript for transcript, coherence in transcript_coherence_scores.iteritems()], ["NA" for transcript, coherence in transcript_coherence_scores.iteritems()]))
		transcript_windowed_posteriors = dict(zip([transcript for transcript, signal in transcript_spectre_signals.iteritems()], ["NA" for transcript, signal in transcript_spectre_signals.iteritems()]))
		transcript_spectre_posteriors = dict(zip([transcript for transcript, score in transcript_spectre_scores.iteritems()], ["NA" for transcript, score in transcript_spectre_scores.iteritems()]))
		logger.info("calculate_transcript_scores(): Populating NAs for split analysis posteriors... [COMPLETE].")		
	else:
		logger.info("calculate_transcript_scores(): Calculating transcript posteriors... [STARTED].")
		# Build partial functions:
		Coherence_post = partial(posterior_probability, translated_coherence_scores, untranslated_coherence_scores)
		Windowed_post = partial(posterior_probability, translated_spectre_scores, untranslated_spectre_scores)
		SPECtre_post = partial(posterior_probability, translated_spectre_scores, untranslated_spectre_scores)
		# Calculate the posterior probabilities for each transcript:
		logger.info("calculate_transcript_scores(): Calculating SPECtre/Coherence transcript posteriors... [STARTED].")
		coherence_posteriors = pool.map(Coherence_post, transcript_coherence_scores.iteritems())
		windowed_posteriors = pool.map(Windowed_post, transcript_spectre_signals.iteritems())
		spectre_posteriors = pool.map(SPECtre_post, transcript_spectre_scores.iteritems())
		# Parition SPECtre/Coherence posteriors to their respective containers:
		transcript_coherence_posteriors = dict(zip([transcript for transcript, posterior in coherence_posteriors], [posterior for transcript, posterior in coherence_posteriors]))
		transcript_windowed_posteriors = dict(zip([transcript for transcript, posterior in windowed_posteriors], [posterior for transcript, posterior in windowed_posteriors]))
		transcript_spectre_posteriors = dict(zip([transcript for transcript, posterior in spectre_posteriors], [posterior for transcript, posterior in spectre_posteriors]))
		logger.info("calculate_transcript_scores(): Calculating SPECtre/Coherence transcript posteriors... [COMPLETE].")
		logger.info("calculate_transcript_scores(): Calculating transcript-level metrics [COMPLETE].")

	#####################################################################################
	# OUTPUT THE TRANSCRIPT COVERAGES, SCORES AND POSTERIORS TO A COMPOSITE DICTIONARY: #
	#####################################################################################
	logger.info("calculate_transcript_scores(): Output transcript metrics to hash()... [STARTED].")
	# Instantiate the composite transcript hash:
	metrics = hash()

	# Populate the hash() with the A-site and P-site read coverages:
	logger.info("calculate_transcript_scores(): Output transcript A/P-site read coverages... [STARTED].")
	if "FLOSS" in methods:
		metrics = populate_hash(metrics, "A_coverage", transcript_floss_asite_coverages.items())
	metrics = populate_hash(metrics, "P_coverage", transcript_spectre_psite_coverages.items())
	logger.info("calculate_transcript_scores(): Output transcript A/P-site read coverages... [COMPLETE].")

	# Populate the hash() with the SPECtre, Coherence, FLOSS, and ORFscore metrics for each transcript:
	logger.info("calculate_transcript_scores(): Output SPECtre scores, signals and posteriors... [STARTED].")
	metrics = populate_hash(metrics, "SPEC_score", transcript_spectre_scores.items())
	metrics = populate_hash(metrics, "SPEC_signal", transcript_spectre_signals.items())
	metrics = populate_hash(metrics, "SPEC_score_posterior", transcript_spectre_posteriors.items())
	metrics = populate_hash(metrics, "SPEC_signal_posterior", transcript_windowed_posteriors.items())
	metrics = populate_hash(metrics, "SPEC_phase", transcript_spectre_phases.items())
	logger.info("calculate_transcript_scores(): Output SPECtre scores, signals and posteriors... [COMPLETE].")

	if "Full" in methods:
		logger.info("calculate_transcript_scores(): Output Coherence scores, signals and posteriors... [STARTED].")
		metrics = populate_hash(metrics, "FULL_score", transcript_coherence_scores.items())
		metrics = populate_hash(metrics, "FULL_signal", transcript_coherence_signals.items())
		metrics = populate_hash(metrics, "FULL_score_posterior", transcript_coherence_posteriors.items())
		logger.info("calculate_transcript_scores(): Output Coherence scores, signals and posteriors... [COMPLETE].")
	if "FLOSS" in methods:
		logger.info("calculate_transcript_scores(): Output FLOSS scores and distributions... [STARTED].")
		metrics = populate_hash(metrics, "FLOSS_score", transcript_floss_scores.items())
		metrics = populate_hash(metrics, "FLOSS_distribution", transcript_floss_distributions.items())
		logger.info("calculate_transcript_scores(): Output FLOSS scores and distributions... [COMPLETE].")
	if "ORFscore" in methods:
		logger.info("calculate_transcript_scores(): Output ORFscores... [STARTED].")
		metrics = populate_hash(metrics, "ORF_score", transcript_orf_scores.items())
		metrics = populate_hash(metrics, "ORF_reads", transcript_orf_reads.items())
		logger.info("calculate_transcript_scores(): Output ORFscores... [COMPLETE].")
	logger.info("calculate_transcript_scores(): Output transcript metrics to hash()... [COMPLETE].")
	if "FLOSS" in methods:
		return metrics, reference_distribution
	else:
		return metrics, "NA"

##### STILL NEED EDITS #####
class ExperimentMetrics(object):

	'''
	Once all transcripts have been scored according to the metrics requested, calculate
	experiment-level metrics including thresholds for protein-coding versus non-coding
	SPECtre scores, translated versus untranslated SPECtre scores, and analytical AUCs
	using the ROCR package in R.
	'''

	def __init__(self, stats, fpkms, analyses, cutoff, fdr):
		self.stats = stats
		self.fpkms = fpkms
		self.analyses = analyses
		self.cutoff = cutoff
		self.fdr = fdr

	@staticmethod
	def build_score_distributions(transcript_stats, transcript_fpkms, fpkm_cutoff, method):
		# Further refinement required after checking input format...
		translated, not_translated = list(), list()
		# Get the scores by transcript and method:
		metric = "_".join([method.upper(), "score"])
		for gene_type in transcript_stats:
			if gene_type == "protein_coding":
				for chrom in transcript_stats[gene_type]:
					for strand in transcript_stats[gene_type][chrom]:
						for gene in transcript_stats[gene_type][chrom][strand]:
							for transcript in transcript_stats[gene_type][chrom][strand][gene]:
								for feature in transcript_stats[gene_type][chrom][strand][gene][transcript]:
									if feature == "CDS":
										if metric in transcript_stats[gene_type][chrom][strand][gene][transcript][feature]:
											if not transcript_stats[gene_type][chrom][strand][gene][transcript][feature][metric] == None:
												if not transcript_stats[gene_type][chrom][strand][gene][transcript][feature][metric] == "NA":
													score = transcript_stats[gene_type][chrom][strand][gene][transcript][feature][metric]
													if transcript in transcript_fpkms:
														if transcript_fpkms[transcript] >= fpkm_cutoff:
															translated.append(score)
														else:
															not_translated.append(score)
		return translated, not_translated

	def translation_threshold(self):
		logger.info("ExperimentMetrics.translation_threshold(): Calculating experiment-level translation threshold [STARTED].")
		translated, untranslated = self.build_score_distributions(self.stats, self.fpkms, self.cutoff, "SPEC")
		r('active <- data.frame(SPEC=c(%s), biotype="translated")' %",".join(str(n) for n in translated))
		r('inactive <- data.frame(SPEC=c(%s), biotype="not_translated")' %",".join(str(n) for n in untranslated))
		r('scores <- rbind(active, inactive)')
		logger.info("ExperimentMetrics.translation_threshold(): Calculating experiment-level translation threshold [COMPLETE].")
		return str(r('quantile(scores$SPEC[scores$biotype=="not_translated"], probs=%s)' %(1-float(self.fdr)))[0])

	def spectre_auc(self):
		logger.info("ExperimentMetrics.spectre_auc(): Calculating experiment-level SPECtre AUC [STARTED].")
		if "SPECtre" not in self.analyses:
			return "NA"
		else:
			# Instantiate necessary R packages:
			rocr = importr("ROCR")
			translated, untranslated = self.build_score_distributions(self.stats, self.fpkms, self.cutoff, "SPEC")
			r('active <- data.frame(SPEC=c(%s), biotype="translated")' %",".join(str(n) for n in translated))
			r('inactive <- data.frame(SPEC=c(%s), biotype="not_translated")' %",".join(str(n) for n in untranslated))
			r('scores <- rbind(active, inactive)')
			logger.info("ExperimentMetrics.spectre_auc(): Calculating experiment-level SPECtre AUC [COMPLETE].")
			return str(r('performance(prediction(scores$SPEC, scores$biotype), "auc")@y.values[[1]]')[0])

	def floss_auc(self):
		logger.info("ExperimentMetrics.floss_auc(): Calculating experiment-level FLOSS AUC [STARTED].")
		if "FLOSS" not in self.analyses:
			return "NA"
		else:
			# Instantiate necessary R packages:
			rocr = importr("ROCR")
			translated, untranslated = self.build_score_distributions(self.stats, self.fpkms, self.cutoff, "FLOSS")
			r('active <- data.frame(SPEC=c(%s), biotype="translated")' %",".join(str(n) for n in translated))
			r('inactive <- data.frame(SPEC=c(%s), biotype="not_translated")' %",".join(str(n) for n in untranslated))
			r('scores <- rbind(active, inactive)')
			logger.info("ExperimentMetrics.floss_auc(): Calculating experiment-level FLOSS AUC [COMPLETE].")
			return str(r('performance(prediction(-scores$SPEC, scores$biotype), "auc")@y.values[[1]]')[0])

	def orfscore_auc(self):
		logger.info("ExperimentMetrics.orfscore_auc(): Calculating experiment-level ORFscore AUC [STARTED].")
		if "ORFscore" not in self.analyses:
			return "NA"
		else:
			# Instantiate necessary R packages:
			rocr = importr("ROCR")
			translated, untranslated = self.build_score_distributions(self.stats, self.fpkms, self.cutoff, "ORF")
			r('active <- data.frame(SPEC=c(%s), biotype="translated")' %",".join(str(n) for n in translated))
			r('inactive <- data.frame(SPEC=c(%s), biotype="not_translated")' %",".join(str(n) for n in untranslated))
			r('scores <- rbind(active, inactive)')
			logger.info("ExperimentMetrics.orfscore_auc(): Calculating experiment-level ORFscore AUC [COMPLETE].")
			return str(r('performance(prediction(scores$SPEC, scores$biotype), "auc")@y.values[[1]]')[0])

	def full_auc(self):
		logger.info("ExperimentMetrics.full_auc(): Calculating experiment-level Spectral Coherence AUC [STARTED].")
		if "Full" not in self.analyses:
			return "NA"
		else:
			# Instantiate necessary R packages:
			rocr = importr("ROCR")
			translated, untranslated = self.build_score_distributions(self.stats, self.fpkms, self.cutoff, "Full")
			r('active <- data.frame(SPEC=c(%s), biotype="translated")' %",".join(str(n) for n in translated))
			r('inactive <- data.frame(SPEC=c(%s), biotype="not_translated")' %",".join(str(n) for n in untranslated))
			r('scores <- rbind(active, inactive)')
			logger.info("ExperimentMetrics.full_auc(): Calculating experiment-level Spectral Coherence AUC [COMPLETE].")
			return str(r('performance(prediction(scores$SPEC, scores$biotype), "auc")@y.values[[1]]')[0])

def print_metrics(output_file, transcript_stats, experiment_stats, reference_distribution, gtf, fpkms, analyses, split, parameters):

	def format_coordinates(coords):
		if len(coords) == 0:
			return "NA"
		else:
			return ",".join(sorted([str(iv.start) + "-" + str(iv.end) for iv in coords]))

	def return_coordinates(gtf, gene_type, chrom, strand, gene, transcript, feature):
		try:
			return format_coordinates(gtf[gene_type][chrom][strand][gene][transcript][feature])
		except KeyError:
			return "NA"

	def return_metric(stats, gene_type, chrom, strand, gene, transcript, feature, metric):
		if gene_type == "protein_coding":
			try:
				if isinstance(stats[gene_type][chrom][strand][gene][transcript][feature][metric], list):
					return ",".join([str(x) for x in stats[gene_type][chrom][strand][gene][transcript][feature][metric]])
				elif isinstance(stats[gene_type][chrom][strand][gene][transcript][feature][metric], dict):
					if len(stats[gene_type][chrom][strand][gene][transcript][feature][metric]) == 0:
						return "NA"
					else:
						return zip(*stats[gene_type][chrom][strand][gene][transcript][feature][metric].items())
				else:
					return stats[gene_type][chrom][strand][gene][transcript][feature][metric]
			except KeyError:
				return "NA"
		else:
			if feature == "CDS":
				try:
					if isinstance(stats[gene_type][chrom][strand][gene][transcript]["exon"][metric], list):
						return ",".join([str(x) for x in stats[gene_type][chrom][strand][gene][transcript]["exon"][metric]])
					elif isinstance(stats[gene_type][chrom][strand][gene][transcript]["exon"][metric], dict):
						if len(stats[gene_type][chrom][strand][gene][transcript]["exon"][metric]) == 0:
							return "NA"
						else:	
							return zip(*stats[gene_type][chrom][strand][gene][transcript]["exon"][metric].items())
					else:
						return stats[gene_type][chrom][strand][gene][transcript]["exon"][metric]
				except KeyError:
					return "NA"
			else:
				try:
					if isinstance(stats[gene_type][chrom][strand][gene][transcript][feature][metric], list):
						return ",".join([str(x) for x in stats[gene_type][chrom][strand][gene][transcript][feature][metric]])
					elif isinstance(stats[gene_type][chrom][strand][gene][transcript][feature][metric], dict):
						if len(stats[gene_type][chrom][strand][gene][transcript][feature][metric]) == 0:
							return "NA"
						else:
							return zip(*stats[gene_type][chrom][strand][gene][transcript][feature][metric].items())
					else:
						return stats[gene_type][chrom][strand][gene][transcript][feature][metric]
				except KeyError:
					return "NA"

	def format_fpkm(transcript_fpkms, transcript):
		if transcript in transcript_fpkms:
			return str(transcript_fpkms[transcript])
		else:
			return "NA"

	def write_parameters(parameters, analyses):
		parameters_string = "\n# Input = " + str(parameters.input) + " \n# Output = " + str(parameters.output) + "\n# FPKM = " + str(parameters.fpkm) + "\n# GTF = " + str(parameters.gtf)
		parameters_string += "\n# PARAMETERS\n# Analyses = " + str(analyses) + "\n# Window Length = " + str(parameters.len) + "\n# FPKM Cutoff = " + str(parameters.min) + "\n# FDR = " + str(parameters.fdr) + "\n# Test = " + str(parameters.type)
		return parameters_string

	def return_experiment_metric(experiment_stats, analysis):
		try:
			if analysis == "Full":
				return str(experiment_stats.full_auc())
			elif analysis == "FLOSS":
				return str(experiment_stats.floss_auc())
			elif analysis == "ORFscore":
				return str(experiment_stats.orfscore_auc())
		except:
			return "NA"

	def write_experiment_metrics(experiment_stats):
		metric_string = "\n# Translation Threshold = " + str(experiment_stats.translation_threshold()) + "\n# SPECtre AUC = " + str(experiment_stats.spectre_auc()) + "\n# Full AUC = " + return_experiment_metric(experiment_stats, "Full") + "\n# FLOSS AUC = " + return_experiment_metric(experiment_stats, "FLOSS") + "\n# ORFscore AUC = " + return_experiment_metric(experiment_stats, "ORFscore")
		return metric_string

	# Instantiate necessary R packages:
	rocr = importr("ROCR")
	graphics = importr("grDevices")

	logger.info("print_metrics/main(): Output results to file... [STARTED]")
	# Initialize the default header text:
	header = "\nid\tchr\tstrand\tgene_id\ttranscript_id\tgene_type\tribo_fpkm\tcoordinates_5UTR\tcoordinates_CDS\tcoordinates_3UTR"
	for feature in ("UTR5", "CDS", "UTR3"):
		for analysis in analyses:
			if analysis == "SPECtre":
				header += "\t" + "\t".join(["SPEC_score_" + feature, "SPEC_score_posterior_" + feature, "SPEC_signal_" + feature, "SPEC_signal_posterior_" + feature])
			if analysis == "Full":
				header += "\t" + "\t".join(["FULL_score_" + feature, "FULL_score_posterior_" + feature, "FULL_signal_" + feature])
			if analysis == "FLOSS":
				header += "\t" + "\t".join(["FLOSS_score_" + feature, "FLOSS_distribution_" + feature])
			if analysis == "ORFscore":
				header += "\t" + "\t".join(["ORF_score_" + feature, "ORF_reads_" + feature])

	if len(target_chroms) > 0:
		output_file.write("# FILE I/O:" + write_parameters(parameters, analyses))
	else:
		output_file.write("# FILE I/O:" + write_parameters(parameters, analyses))
		output_file.write("\n# METRICS:" + write_experiment_metrics(experiment_stats))
		if isinstance(reference_distribution, dict):
			output_file.write("\n" + str(zip(*reference_distribution.items())))
	output_file.write(header)

	count = 1
	for gene_type in sorted(transcript_stats):
		for chrom in sorted(transcript_stats[gene_type]):
			for strand in transcript_stats[gene_type][chrom]:
				for gene in sorted(transcript_stats[gene_type][chrom][strand]):
					for transcript in sorted(transcript_stats[gene_type][chrom][strand][gene]):
						if gene_type == "protein_coding":
							line = "\n" + "\t".join(str(field) for field in [count, chrom, strand, gene, transcript, gene_type, format_fpkm(fpkms, transcript),
									return_coordinates(gtf, gene_type, chrom, strand, gene, transcript, "UTR5"),
									return_coordinates(gtf, gene_type, chrom, strand, gene, transcript, "CDS"),
									return_coordinates(gtf, gene_type, chrom, strand, gene, transcript, "UTR3")])
						else:
							line = "\n" + "\t".join(str(field) for field in [count, chrom, strand, gene, transcript, gene_type, format_fpkm(fpkms, transcript),
									return_coordinates(gtf, gene_type, chrom, strand, gene, transcript, "UTR5"),
									return_coordinates(gtf, gene_type, chrom, strand, gene, transcript, "exon"),
									return_coordinates(gtf, gene_type, chrom, strand, gene, transcript, "UTR3")])							
						for feature in ("UTR5", "CDS", "UTR3"):
							for analysis in analyses:
								if analysis == "SPECtre":
									line += "\t" + "\t".join([str(return_metric(transcript_stats, gene_type, chrom, strand, gene, transcript, feature, "SPEC_score")),
											str(return_metric(transcript_stats, gene_type, chrom, strand, gene, transcript, feature, "SPEC_score_posterior")),
											str(return_metric(transcript_stats, gene_type, chrom, strand, gene, transcript, feature, "SPEC_signal")),
											str(return_metric(transcript_stats, gene_type, chrom, strand, gene, transcript, feature, "SPEC_signal_posterior"))])
								if analysis == "Full":
									line += "\t" + "\t".join([str(return_metric(transcript_stats, gene_type, chrom, strand, gene, transcript, feature, "FULL_score")),
											str(return_metric(transcript_stats, gene_type, chrom, strand, gene, transcript, feature, "FULL_score_posterior")),
											str(return_metric(transcript_stats, gene_type, chrom, strand, gene, transcript, feature, "FULL_signal"))])
								if analysis == "FLOSS":
									line += "\t" + "\t".join([str(return_metric(transcript_stats, gene_type, chrom, strand, gene, transcript, feature, "FLOSS_score")),
											str(return_metric(transcript_stats, gene_type, chrom, strand, gene, transcript, feature, "FLOSS_distribution"))])
								if analysis == "ORFscore":
									line += "\t" + "\t".join([str(return_metric(transcript_stats, gene_type, chrom, strand, gene, transcript, feature, "ORF_score")),
											str(return_metric(transcript_stats, gene_type, chrom, strand, gene, transcript, feature, "ORF_reads"))])
						output_file.write(line)
						count += 1
	logger.info("print_metrics/main(): Output results to file... [COMPLETE]")

if __name__ == "__main__":
	
	# Parse arguments from command-line:
	parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS, add_help=False)
	parser.add_argument("--help", action="help", help="show this help message and exit")
	parser.add_argument("--full", action="store_true", default="store_false", help="calculate un-windowed spectral coherence")
	parser.add_argument("--floss", action="store_true", default="store_false", help="calculate FLOSS and distribution")
	parser.add_argument("--orfscore", action="store_true", default="store_false", help="calculate ORFScore and distribution")
	spectre_args = parser.add_argument_group("parameters for SPECtre analysis:")
	spectre_args.add_argument("--nt", action="store", required=False, nargs="?", default=1, metavar="INT", help="number of threads for multi-processing (default: %(default)s)")
	spectre_args.add_argument("--len", action="store", required=False, nargs="?", default=30, metavar="INT", help="length of sliding window (default: %(default)s)")
	spectre_args.add_argument("--min", action="store", required=False, nargs="?", default=3, metavar="FLOAT", help="minimum FPKM for active translation (default: %(default)s FPKM)")
	spectre_args.add_argument("--fdr", action="store", required=False, nargs="?", default=0.05, metavar="FLOAT", help="FDR cutoff (default: %(default)s)")
	spectre_args.add_argument("--step", action="store", required=False, nargs="?", default=3, metavar="INT", help="distance between sliding windows (default: %(default)s)")
	spectre_args.add_argument("--type", action="store", required=False, nargs="?", default="median", metavar="TYPE", choices=["mean","median","max","nonzero_mean","nonzero_median"], help="metric for SPECtre analysis (choices: mean,[median],max,nonzero_mean,nonzero_median)")
	spectre_args.add_argument("--target", action="store", required=False, nargs="?", default="", metavar="LIST", help="specify a single chromosome ('X') or comma-delimited list of chromosomes to enable split chromosome analysis (faster)")
	spectre_args.add_argument("--offsets", action="store", required=False, nargs="?", default="", metavar="LIST", help="comma-delimited user-defined list of read_length:offset_position definitions (eg. 28:12,29:14,30:15,31:15,32:15)")
	file_args = parser.add_argument_group("input and output parameters:")
	file_args.add_argument("--input", action="store", required=True, nargs="?", metavar="BAM", type=str, help="location of BAM alignment file")
	file_args.add_argument("--output", action="store", required=True, nargs="?", default="./spectre_results.txt", metavar="FILE", help="write results to (default: %(default)s)")
	file_args.add_argument("--fpkm", action="store", required=True, nargs="?", metavar="FILE", type=str, help="location of Cufflinks isoforms.fpkm_tracking file")
	file_args.add_argument("--gtf", action="store", required=True, nargs="?", metavar="FILE", type=str, help="location of GTF annotation file")
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

	# Transcript annotation types were pulled from Ensembl under the entry
	# "Biotype" at http://useast.ensembl.org/Help/Glossary. Mitochondrial and
	# those related to rRNA not included:
	coding = ["protein_coding"]
	noncoding = ["miRNA", "misc_RNA", "ncRNA", "scRNA", "snlRNA", "snoRNA", "snRNA"]
	lncrna = ["lincRNA", "lncRNA"]
	pseudogene = ["process_pseudogene", "pseudogene", "transcribed_processed_pseudogene",
					"transcribed_unprocessed_pseudogene", "translated_processed_pseudogene",
					"translated_unprocessed_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene"]

	# Default transcript boundary buffers:
	psite_buffers = {"CDS": (0,0), "UTR": (0,0)}
	asite_buffers = {"CDS": (45,15), "UTR": (15,15)}
	orfscore_buffers = {"CDS": (3,3), "UTR": (3,3)}

	# Offset positions for A/P sites:
	offsets = {"ingolia": {26: 15, 27: 14, 28: 14, 30: 15, 31: 15}, "bazzini": {26: 12, 27: 12, 28: 12, 29: 12, 30: 12, 31: 12}}
	if len(args.offsets) > 0:
		offsets = add_custom_offsets(args.offsets, offsets)

	# Initial check of chromosome format in BAM, GTF and Cufflinks input:
	logger.info("main(): Checking validity of input files... [STARTED].")
	if check_chromosomes(args.input, args.gtf, args.fpkm, args.target) == True:
		logger.info("main(): Checking validity of input files... [COMPLETE].")

		# Initialize the chromosome(s) to be targeted for analysis:
		target_chroms = [str(chrom) for chrom in args.target.split(",") if len(chrom) > 0]
		# Extract transcripts, transcript intervals and expression using the provided GTF:
		transcript_fpkms = extract_fpkms(args.fpkm, target_chroms)
		transcript_array, transcript_intervals, transcript_gtf = parse_gtf(args.gtf, transcript_fpkms, target_chroms)

		# Initialize the types of analyses to be conducted (default: SPECtre):
		analyses = ["SPECtre"]
		if args.full == True:
			analyses.append("Full")
		if args.floss == True:
			analyses.append("FLOSS")
		if args.orfscore == True:
			analyses.append("ORFscore")


		# Calculate the designated transcript-level scores based on the analyses to be conducted:
		transcript_metrics, reference_read_distribution = calculate_transcript_scores(transcript_gtf, transcript_fpkms, float(args.min), asite_buffers, psite_buffers, orfscore_buffers, args.input, int(args.len), int(args.step), args.type, analyses, offsets, target_chroms, int(args.nt))
		
		# Perform a second-pass global analysis based on the transcript-level metrics, such as:
		# ROC analyses, posterior probability as a function of empirical FDR, and codon window
		# signal plots (based on windowed spectral coherence).
		experiment_metrics = ExperimentMetrics(transcript_metrics, transcript_fpkms, analyses, float(args.min), float(args.fdr))
		# Print the results table to the output file:
		print_metrics(open(args.output,"w"), transcript_metrics, experiment_metrics, reference_read_distribution, transcript_gtf, transcript_fpkms, analyses, target_chroms, args)
	else:
		print "[ERROR]: Please check that your BAM, GTF and Cufflinks input files are formatted correctly..."
	sys.exit()
