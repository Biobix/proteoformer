README Mapping/Parsing/General (GerbenM)

1° General prerequisites
------------------------

- Perl v.5.10
- SQLite
- Getopt::Long package
- DBI package
- fastx toolkit (http://hannonlab.cshl.edu/fastx_toolkit/download.html) if using Bowtie mapper
- STAR, TopHat2, Bowtie(2) aligner


2° Mapping of the FastQ files
-----------------------------

Prerequisites:

 - FastQ files of all the sequenced lanes (with Illumina GAIIx or HiSEQ). The FastQ files ($file) need to be placed in the a fastq folder under the $workDir. The actual looping over the different sequencing fastq input files is handled in a separate shell script (e.g. mapping_SM.sh and mapping_CL.sh, different if run on single machine or cluster environment).

 - Igenome of organism (http://tophat.cbcb.umd.edu/igenomes.html) 

 - The Igenome folder structure contain all necessary information for the creation of the mapping indexes (if using transcriptome mappers: STAR and/or TopHat2). In a layered approach (using Bowtie1/Bowtie2, one still needs to download the complete transcriptome/cDNA from e.g. Ensembl and needs to create bowtie index files (manually created with bowtie_build based on download from eg. http://feb2012.archive.ensembl.org/info/data/ftp/index.html. Needs to be placed under the iGenome subdirectory holding the BOWTIE indexes with the name eg. Mus_musculus.NCBIM37.66.cdna.all.fa).

 - The STAR indexes are automatically created if they're not found:
	- for the complete genome
	- for the rRNA sequences (one still needs to supply an fasta file containing the rRNA sequences (places within the Igenome folder structure under ~/Sequence/AbundantSequences/hsa_rRNA_seqs.fasta)
############ TODO Create rRNA indexes from SQLite entries
	- for the phix control (this fasta file is already present in the Igenome folder structure under ~/Sequence/AbundantSequences)
	- the STAR indexes are created under ~/Sequence/STARIndex/

- The Bowtie1 and Bowtie2 genome indexes are already present in the Igenome folder structure under ~/Sequence/BowtieIndex and ~/Sequence/Bowtie2Index
  
The input variables of the mapping script are:

GetOptions(
"dir=s"=>\$work_dir,            # Path to the working directory,                                    mandatory argument
"file=s"=>\$seqFileName,        # Base name of the fastq file, e.g. "lane",                         mandatory argument
"name=s"=>\$run_name,           # Name of the run,                                                  mandatory argument
"species=s"=>\$species,         # Species, eg mouse/human/fruitfly,                                 mandatory argument
"ensembl=i"=>\$ensemblversion,  # Ensembl annotation version, eg 70 (Feb2012),                      mandatory argument
"cores=i"=>\$cores,             # Number of cores to use for Bowtie Mapping,                        mandatory argument
"readtype=s"=>\$readtype,       # The readtype (ribo, PE_polyA, SE_polyA, PE_total, SE_total)       mandatory  argument
"mapper:s"=>\$mapper,           # The mapper used for alignment (Bowtie,Bowtie2,STAR,TopHat2)       optional argument   (default = STAR)
"readlength:i"=>\$readlength,   # The readlength (if RiboSeq take 50 bases),                        optional  argument  (default = 50)
"adaptor:s"=>\$adaptorSeq,      # The adaptor sequence that needs to be clipped with fastx_clipper, optional  argument  (default = CTGTAGGCACCATCAATAGATCGGA) => Ingolia paper
"unique=s" =>\$unique,          # Retain the uniquely (and multiple) mapping reads (Y or N),        mandatory  argument
"tmp:s" =>\$tmpfolder         # Folder where temporary files are stored,                           optional  argument (default = $TMP env setting)
);

Usage example: 
./1_mapping.pl --dir /data/RIBO_runs/RIBO_HCT116_GerbenM/ --file lane1 --name HCT116 --species human --ensembl 70 --cores 24 --readtype ribo --unique N ( --mapper STAR --readlength 50 --tmp /data/RIBO_runs/RIBO_HCT116_GerbenM/tmp/ --adaptor CTGTAGGCACCATCAATAGATCGGA )


3° Parsing the profile locations
--------------------------------

Prerequisites:

- Perl package Parallel::ForkManager to parallellize process over multiple chromosomes
- an SQLite database containing necessary information on the annotation bundle (e.g. Ensembl 70, human), called e.g. ENS_HSA_70.db


4° Output
---------

- The mapping output, placed under $work_dir/$mapper/$file
- A BED and BEDGRAPH file holding the genome-wide ribosome-profile densities (only the first base of the A-site is used), placed under $work_dir/output
- An SQLite database holding this genome-wide ribosome-profile information, next to some mapping statistics, placed under $work_dir/SQLite




