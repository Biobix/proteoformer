#!/usr/bin/python

import sys
import re
import logging
import argparse

# DESCRIPTION
# This script takes as input: 1) a GTF derived from the UCSC knownGene database, 2) an Ensembl GTF of the
# same genome version, and 3) a tab-delimited file that maps UCSC to Ensembl IDs. Instructions to generate
#both files are detailed below. Instructions as listed below are for human references.

# Derive a GTF annotation file from the UCSC knownGene database (from the Linux/Unix command-line):
# 1) Download the genePredToGtf binary from http://hgdownload.cse.ucsc.edu/admin/exe/
# 2) Navigate to the folder that contains the genePredToGtf binary and execute the following:
#	chmod 600 genePredToGtf
# 3) In your $HOME directory, create and save the following 4-line configuration file (shown using VI):
# 	cd
#	vi .hg.conf
#	db.host=genome-mysql.cse.ucsc.edu
#	db.user=genomep
#	db.password=password
#	central.db=hgcentral
#	:wq  
# 4) Set permissions for the configuration file:
#	chmod 600 .hg.conf
# 5) Navigate to the folder that contains your knownGene file and execute the following:
#	genePredToGtf -utr -honorCdsStat -source=hg38 hg38 knownGene knownGene.gtf
#
# NOTES:
# Replace the two references to "hg38" in genePredToGtf executable to your desired genome version. This
# genome version should match the Ensembl/GENCODE version used in the mapping file (see below).

# Generate a UCSC to Ensembl mapping file (from your WWW browser):
# 1) Navigate to the UCSC Table Browser at http://genome.ucsc.edu/cgi-bin/hgTables
# 2) For this example, set the following fields to:
#	clade=Mammal
#	genome=Human
#	assembl=hg38
#	group=Genes and Gene Predictions
#	track=GENCODE v22
#	table=knownToEnsembl
#	output_format=all fields from selected table
#	output_file=ensGene.txt
#
# This will generate a two-column tab-delimited file with UCSC transcript IDs to Ensembl transcript IDs.

# USAGE: python map_ucsc_to_ensembl.py --ucsc <path/to/knownGene.gtf> --ensembl <path/to/ensembl.gtf> --mapfile <path/to/ensGene.txt> > <path/to/mapped.gtf>

class Inputs(object):
	
	'''
	Take as in put a GTF derived from the UCSC knownGene database and map Ensembl identifiers, and other annotation variables
	using a mapping file generated from the ensGene database.
	'''

	def __init__(self, ensembl_file, map_file):
		self.ensembl_file = ensembl_file
		self.map_file = map_file

	def mapped_ids(self):
		logger.info("Mapping UCSC to Ensembl transcript IDs from: " + self.map_file + " [STARTED]")
		mapped = dict()
		for line in open(self.map_file):
			if not line.startswith("#"):
				ucsc_id, ensembl_id = line.strip().split("\t")
				if not ucsc_id in mapped:
					mapped[ucsc_id] = ensembl_id.split(".")[0]	# Transcript version number not necessary.
		logger.info("Mapping UCSC to Ensembl transcript IDs from: " + self.map_file + " [COMPLETE]")
		return mapped

	def ensembl_gtf(self):
		logger.info("Parsing Ensembl GTF: " + self.ensembl_file + " into memory [STARTED].")
		ensembl = dict()
		for line in open(self.ensembl_file):
			if not line.startswith("#"):
				seqname, source, feature, start, end, score, strand, frame, attributes = line.strip().split("\t")
				if feature in ("CDS", "exon", "UTR", "5UTR", "3UTR"):
					gene_id = re.findall("ENS[A-Z]*G[0-9]*", attributes)[0]
					gene_type = re.findall("gene_biotype .*;{1}", attributes)[0].split('"')[1]
					transcript_id = re.findall("ENS[A-Z]*T[[0-9]*", attributes)[0]
					transcript_type = re.findall("transcript_biotype .*;{1}", attributes)[0].split('"')[1]
					if transcript_id not in ensembl:
						ensembl[transcript_id] = {"gene": gene_id, "gene_biotype": gene_type, "transcript_biotype": transcript_type}
		logger.info("Parsing Ensembl GTF: " + self.ensembl_file + " into memory [COMPLETE].")
		return ensembl
			
if __name__ == "__main__":
	
	# Parse input and output arguments from command-line:
	parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS, add_help=False)
	parser.add_argument("--help", action="help", help="show this help message and exit")
	file_args = parser.add_argument_group("required input and output parameters:")
	file_args.add_argument("--ucsc", action="store", required=True, nargs="?", metavar="GTF", type=str, help="location of knownGene.gtf file")
	file_args.add_argument("--ensembl", action="store", required=True, nargs="?", metavar="GTF", type=str, help="location of Ensembl GTF file")
	file_args.add_argument("--mapfile", action="store", required=True, nargs="?", metavar="TXT", type=str, help="location of ensGene.txt file")
	file_args.add_argument("--log", action="store", nargs="?", metavar="TXT", type=str, default="./ensGene.log", help="location of log file (default: %(default)s)")
	args = parser.parse_args()

	# Enable event logging:
	logger = logging.getLogger(__name__)
	logger.setLevel(logging.INFO)
	handler = logging.FileHandler(args.log)
	handler.setLevel(logging.INFO)
	formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
	handler.setFormatter(formatter)
	logger.addHandler(handler)

	# Load the Ensembl GTF and mapping file into memory:
	parsed = Inputs(args.ensembl, args.mapfile)
	mapping = parsed.mapped_ids()
	ensembl = parsed.ensembl_gtf()

	# Get basic starting statistics:
	ensembl_count = len(ensembl)
	ucsc_count = len(mapping)
	
	# Initialize basic counters:
	skipped_ensembl = list()
	skipped_ucsc = list()

	for line in open(args.ucsc):
		if not line.startswith("#"):
			seqname, source, feature, start, end, score, strand, frame, attributes = line.strip().split("\t")
			seqname = seqname.replace("chr", "")
			gene_id = re.findall("gene_id .*;{1}", attributes)[0].split('"')[1]
			transcript_id = re.findall("transcript_id .*;{1}", attributes)[0].split('"')[1]
			if transcript_id in mapping:
				ensembl_id = mapping[transcript_id]
				if ensembl_id in ensembl:
					attributes += " ensembl_gene_name " + ensembl[ensembl_id]["gene"] + "; ensembl_transcript_name " + ensembl_id + "; gene_biotype " + ensembl[ensembl_id]["gene_biotype"] + "; transcript_biotype " + ensembl[ensembl_id]["transcript_biotype"] + ";"
					print seqname + "\t" + source + "\t" + feature + "\t" + start + "\t" + end + "\t" + score + "\t" + strand + "\t" + frame + "\t" + attributes
				else:
					if ensembl_id not in skipped_ensembl:
						skipped_ensembl.append(ensembl_id)
			else:
				if transcript_id not in skipped_ucsc:
						skipped_ucsc.append(transcript_id)

	logger.info("Ensembl Transcripts = " + str(ensembl_count))
	logger.info("UCSC Transcripts = " + str(ucsc_count))
	logger.info("Ensembl Transcripts Skipped = " + str(len(skipped_ensembl)))
	logger.info("UCSC Transcripts Skipped = " + str(len(skipped_ucsc)))
	

