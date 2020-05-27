# SPECtre

## Description
This software is designed to identify regions of active translation from ribosome profiling sequence data. This analytical pipeline scores the translational status of each annotated region (5'UTR, CDS, exon, 3'UTR) as a function of its spectral coherence over user-defined N nucleotide sliding windows to an idealized reference coding signal. It is used by calling *SPECtre.py* with the listed parameters.

## Required Resources:
```
R:		https://www.r-project.org/
Rpy:		http://rpy.sourceforge.net
ROCR:		https://rocr.bioinf.mpi-sb.mpg.de/
Python:		https://www.python.org
HTSeq:		https://pypi.python.org/pypi/HTSeq
NumPy:		https://pypi.python.org/pypi/numpy              
pysam:		https://pypi.python.org/pypi/pysam/
pyfasta:	https://pypi.python.org/pypi/pyfasta/
```

## Quick Start
Download and Install:
```
git clone git@github.com:mills-lab/spectre.git
cd spectre
chmod +x SPECtre.py
```

Index Alignments:
```
samtools index <in.bam>
```

Run SPECtre with default parameters:
```
python SPECtre.py \
	--input <in.bam> \
	--output <spectre_results.txt> \
	--log <spectre_results.log> \
	--fpkm <isoforms.fpkm_tracking> \
	--gtf <ensembl.gtf>
```

## Supporting Files
Sample BAM alignment file, Cufflinks output, and Ensembl-formatted GTF are available for testing purposes in the folder *test*.

## Output
*SPECtre* outputs transcript-level and experiment-level translational metrics in tab-delimited text format. Example output is shown in the folder *test*.

## Usage
```
python SPECtre.py [parameters]
```

### Parameters:

#### Required File Arguments:
```
	--input, alignment file in BAM format
	--output, file to output results
	--log, file to track progress
	--fpkm, location of isoforms.fpkm_tracking file from Cufflinks output
	--gtf, location of annotation file in GTF format (only Ensembl supported currently)
```

#### User-defined Analytical Arguments:
```
	--nt <INTEGER>, number of threads for multi-processing (default: 1)
	--len <INTEGER>, length in nucleotides of sliding window for SPECtre analysis (default: 30 nt)
	--min <FLOAT>, minimum FPKM or reads required for classification as active translation (default: 5 FPKM)
	--fdr <FLOAT>, FDR cutoff to use for calculation of posterior probabilities (default: 0.05)
	--step <INT>, distance between sliding windows (default: 3)
	--target <LIST>, specify a single chromosome ("X") or comma-delimited list of chromosomes, to enable splic analysis (faster)
	--offsets <LIST>, comma-delimited list of read_length:offset_position definitions (eg. "28:12,29:14,30:15,31:15,32:15")
	--type <STRING>, summary statistic to use for SPECtre score (default: median)

```

#### Optional Arguments:
```
	--full, enables calculation of un-windowed spectral coherence over full length of transcript
	--floss, enables calculation of FLOSS metric (Ingolia, 2014) for each transcript
	--orfscore, enables calculation of ORFscore (Bazzini, 2014) for each transcript
```

## Test Data:
Test data is derived from ribosome profiling of human SH-SY5Y cells and limited to a single chromosome (3) for space allocation purposes. Similarly, the *.fpkm_tracking input file and test GTF are limited to a single choromosome. To run test analysis with default parameters:
```
python SPECtre.py \
	--input test/test.bam \
	--output test/spectre_test.txt \
	--log test/spectre_test.log \
	--fpkm test/isoforms-test.fpkm_tracking \
	--gtf test/Homo_sapiens.GRCh38.78.test.gtf \
	--full \
	--floss \
	--orfscore
```

##E xample Analytical Pipeline:
Example SPECtre analysis of mESC (Ingolia, 2014) ribosome profiling data.

### Download Required Sequence and Annotation Files
#### Transcript GTF:
```
	ftp://ftp.ensembl.org/pub/release-78/gtf/mus_musculus/
	
	Download the archived Mus_musculus.GRCm38.78.gtf file.
```
#### Reference FASTAs:
```
	ftp://ftp.ensembl.org/pub/release-78/fasta/mus_musculus/dna/
	
	Download the unmasked genomic FASTA files (1-19, MT, X and Y), files should be named according to the format:
	Mus_musculus.GRCm38.dna.chromosome.N.fa.gz (where N = chromosome ID, see above)
	
	Individual chromosome FASTA files may be concatenated for ease-of-use by downstream applications (see: Index Genomic FASTAs)
```
#### rRNA Contaminant FASTA:
```
	wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz

	Extract the appropriate musRibosomal.fa file from the contaminants folder.
```
Note the directory location of these files for future steps.

### Index Genomic and Ribosomal Sequences
#### Index Genomic FASTAs:
```
	bowtie-build <reference_fasta_files> <genomic_index_name>
	
	If TopHat2 is required for alignment, bowtie2 indexes must also be built.
```
#### Index rRNA Contaminants:
```
	bowtie-build <rRNA_fasta_files> <rRNA_index_name>
```

### Library Pre-processing and Alignment
#### Remove Adapters and Trim Library Sequences:
```
	fastx_clipper -Q33 -a CTGTAGGCACCATCAAT -l 24 -c -n –v -i <fastq_file> 2> logs/<log_file> > <library_clipped.fastq>
	fastx_trimmer -Q33 -f 2 -m 24 -i <library_clipped.fastq> > <library_trimmed.fastq>
	
	fastx_clipper and fastx_trimmer are part of the FASTX Toolkit (http://hannonlab.cshl.edu/fastx_toolkit/)
```
#### Align Trimmed Sequences to rRNA Contaminant Database:
```
	bowtie -p 8 -l 23 --solexa-quals --maqerr=60 -S --un <unmapped_reads.fastq> /path/to/<rRNA_index_name> /path/to/<library_trimmed.fastq> 2>> path/to/<log_file> > <rRNA_reads.sam>
```
#### Align to Genome and Transcript GTF:
```
	tophat -p 8 --bowtie1 --solexa-quals --GTF /path/to/Mus_musculus.GRCm38.78.gtf --no-novel-juncs --library-type fr-unstranded -o /path/to/output /path/to/<genomic_index_name> /path/to/<unmapped_reads.fastq>
```

### Alignment QC and Abundance Estimation
#### Filter Alignments Based on Quality:
```
	samtools view -b -q 10 /path/to/tophat/<accepted_hits.bam> > <filtered_hits.bam>
```
#### Estimated RPF Abundance Over Transcripts:
```
	cufflinks -p 4 -o /path/to/output -g /path/to/Mus_musculus.GRCm38.78.gtf /path/to/<filtered_hits.bam>
	
	Alternatively, the -G/--GTF argument may be used instead of -g to estimate abundance only against transcripts annotated in the GTF.
```

### SPECtre Analysis
```
	# CALCULATE CUSTOM P-SITE OFFSETS (IF NECESSARY):
	python calculate_psite_offsets.py /path/to/Mus_musculus.GRCm38.78.gtf /path/to/<filtered_hits.bam>

	# SPECTRE COMMAND LINE INPUT:
	python SPECtre.py \
		--input /path/to/<filtered_hits.bam> \
		--output /path/to/<spectre_output.txt> \
		--log /path/to/<spectre_results.log> \
		--gtf /path/to/Mus_musculus.GRCm38.78.gtf \
		--fpkm /path/to/cufflinks/<isoforms.fpkm_tracking>
	
	Default parameters are as follows:
		--len 30		# Specifies the length of the sliding window used.
		--min 5.0		# Enforces the designated FPKM cutoff for translational status.
		--fdr 0.05		# FDR cutoff for translation threshold.
		--step 3		# Number of nucleotides between sliding windows.
		--type median	# Metric for SPECtre analysis (median, mean, maximum, etc.).

	Alternate P-site offsets can be supplied with the following flag:
		--offsets 28:12,29:12,30:12,31:14,32:15
	
	Additional analyses may be specified by inclusion of the following flags:
		--floss
		--orfscore
		--full
```
#### Expected Output:
```
	Results will be output to the designated file with the general format defined by the number and types of additional analyses requested, for example:

	A summart of parameters used, and experiment-level metrics are output as comments at the beginning of the output file.

	In tab-delimited format, the columns output by SPECtre should be:
		1) unique ID
		2) chromosome
		3) strand
  		4) Ensembl gene ID
  		5) Ensembl transcript ID
	  	6) gene biotype
	  	7) ribosome profiling FPKM
	  	9) spectre_metric [5'UTR]
	   10) spectre_posterior_probability [5'UTR]
	   11) spectre_metric [CDS]
	   12) spectre_posterior_probability [CDS]
	   13) spectre_metric [3'UTR]
	   14) spectre_posterior_probability [3'UTR]
	   ---
	   Full Transcript Spectral Coherence, FLOSS, and ORFscore will follow the same generalized output format as that of SPECtre (above).
