#!/usr/bin/python

import sys
import os
import re
import math
import itertools
import collections

gtf = sys.argv[1]
bam_file = sys.argv[2]

def hash():
	return collections.defaultdict(hash)

def decode(flag):
	if int(flag) & 16:
		return "-"
	else:
		return "+"

def mode(read, positions):
	# Offset should be approximately half the read length, therefore offsets of both shorter
	# and longer length should be ignored:
	try:
		return max([p for p in positions if p>=0.4*read and p<=0.6*read], key=positions.count)
	except ValueError:
		# Default to half the read length:
		return read / 2 

def transcript_coordinates(coords):
	coordinates = list()
	if "," in coords:
		coords = coords.split(",")
	if isinstance(coords, str):
		start, end = coords.split("-")
		region = range(int(start), int(end)+1)
		for pos in region:
			coordinates.append(pos)
	elif isinstance(coords, list):
		for coord in coords:
			if isinstance(coord, str):
				start, end = coord.split("-")
			else:
				start, end = coord
			region = range(int(start), int(end)+1)
			for pos in region:
				coordinates.append(pos)
	else:
		start, end = coords
		region = range(int(start), int(end)+1)
		for pos in region:
			coordinates.append(pos)
	return sorted(coordinates)

def extract_coordinates_from_cigar(pos, cigar):
	# Return a list of expanded coordinates based on the read position and its CIGAR string:
	coordinates = list()
	ops = re.findall("[0-9]*[DHIMNPSX=]{1}", cigar)
	for op in ops:
		increment, modifier = int(op[:-1]), op[-1]
		if not modifier == "N":
			increments = range(increment)
			for x in xrange(len(increments)):
				coordinates.append(pos + increments[x])
		else:
			pos = coordinates[-1] + increment
	return coordinates

class SAM(object):

	# Parse the default fields from a BAM/SAM alignment record:
	def __init__(self, sam):
		self.sam = sam
	def qname(self):
		return self.sam.strip().split("\t")[0]
	def flag(self):
		return self.sam.strip().split("\t")[1]
	def rname(self):
		return self.sam.strip().split("\t")[2]
	def pos(self):
		return self.sam.strip().split("\t")[3]
	def mapq(self):
		return self.sam.strip().split("\t")[4]
	def cigar(self):
		return self.sam.strip().split("\t")[5]
	def rnext(self):
		return self.sam.strip().split("\t")[6]
	def pnext(self):
		return self.sam.strip().split("\t")[7]
	def tlen(self):
		return self.sam.strip().split("\t")[8]
	def seq(self):
		return self.sam.strip().split("\t")[9]
	def qual(self):
		return self.sam.strip().split("\t")[10]
	def optional(self):
		return self.sam.strip().split("t")[11:]

def parse_gtf(gtf_file):
	transcripts = hash()
	for line in open(gtf_file):
		if not line.startswith("#"):
			seq_name, source, feature, start, end, score, strand, frame, attributes = line.strip().split("\t")
			gene_type = re.findall('gene_[a-z]{0,3}type\s\"[\w\.]+\"', attributes)[0].split('"')[1]
			# Parse annotated protein-coding CDS into the GTF dictionary:
			if gene_type == "protein_coding" and feature == "start_codon":
				gene, transcript = re.findall('gene_id\s\"[\w\_\:\:\.]+\"', attributes)[0].split('"')[1], re.findall('transcript_id\s\"[\w\W]+\"', attributes)[0].split('"')[1]
				transcripts[gene_type][seq_name][strand][gene][transcript][feature] = (int(start), int(end))
	# Filter transcripts GTF based on minimum length, expression, and sanitize of overlapping transcripts, as requested:
	print "# GTF parsed..."
	return transcripts

start_codons = parse_gtf(gtf)
start_reads = hash()
for gene_type in start_codons:
	if gene_type == "protein_coding":
		for chrom in start_codons[gene_type]:
			for strand in start_codons[gene_type][chrom]:
				for gene in start_codons[gene_type][chrom][strand]:
					for transcript in start_codons[gene_type][chrom][strand][gene]:
						if "start_codon" in start_codons[gene_type][chrom][strand][gene][transcript]:
							start, end = start_codons[gene_type][chrom][strand][gene][transcript]["start_codon"]
							reads = os.popen("samtools view " + bam_file + " " + chrom + ":" + str(start) + "-" + str(end))
							for read in reads:
								bam = SAM(read)
								read_strand = decode(bam.flag())
								if read_strand == strand:
									if read_strand == "+":
										read_start = int(bam.pos())
										read_end = int(bam.pos()) + len(bam.seq()) - 1
										read_coordinates = extract_coordinates_from_cigar(read_start, bam.cigar())
										try:
											atg_index = read_coordinates.index(int(start))
											if len(bam.seq()) in start_reads:
												if not atg_index == 0:
													start_reads[len(bam.seq())].append(atg_index+1)
											else:
												if not atg_index == 0:
													start_reads[len(bam.seq())] = [atg_index+1]
										except ValueError:
											pass
									else:
										read_start = int(bam.pos())
										read_end = int(bam.pos()) + len(bam.seq()) - 1
										read_coordinates = extract_coordinates_from_cigar(read_start, bam.cigar())[::-1]
										try:
											atg_index = read_coordinates.index(int(start))
											if len(bam.seq()) in start_reads:
												if not atg_index == 0:
													start_reads[len(bam.seq())].append(atg_index+1)
											else:
												if not atg_index == 0:
													start_reads[len(bam.seq())] = [atg_index+1]
										except ValueError:
											pass
print "# Reads parsed..."

print "# Printing offsets..."
for read_length in sorted(start_reads):
	print str(read_length) + "\t" + str(mode(read_length, start_reads[read_length]))