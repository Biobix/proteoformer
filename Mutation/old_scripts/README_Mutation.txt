README Mutation (AlexanderK)


	Ribosomal Sequencing

	-- looking for SNPs in ribosomal sequencing reads



There are two options: use TopHat (and bowtie) to map the ribo-seq reads or use the STAR aligner. The script that uses bowtie takes about 14 hours to finish, whereas the STAR script needs only 3 hours.

- TO DO
Cufflinks doesn't seem to work when using STAR. Maybe something about the file output format of STAR? With cufflinks, the STAR script would of course take more time to finish.
There's a difference in number of SNPs discovered between bowtie (27434) and STAR (21530).


#######
# STAR
#######

what you need:
--------------
STAR
STAR genome annotation/index directory (optional)
cufflinks (optional)
perl
perl DBI package
samtools
fasta reference sequence of the genome
SQLite

how to use it:
--------------
perl findRiboSeqSNPsSTAR.pl --riboreads path/to/sequencing/reads --genome path/to/genome/sequence --stargenome path/to/STAR/genome/annotation --cufflinks

command line options:
-r --riboreads: specify the path to the ribo-seq reads
-g --genome: specify the path to the fasta reference genome sequence
-s --stargenome: specify the path to the STAR genome annotation/index (optional)
-c --cufflinks: flag to indicate if cufflinks should be used to build the transcripts

output:
-------
STAR output files
Samtools output files (snpRiboSeq.vcf = SNP file)

SQLite database:
--> snps table
id | chromosome | position | reference base | alternative base


#########
# Tophat
#########

what you need:
--------------
bowtie index of the genome
fasta reference sequence of the genome
(and some patience)
perl
perl DBI package
SQLite
samtools
tophat
bowtie


how to use it:
--------------
perl findRiboSeqSNPs.pl path/to/bowtie/index path/to/sequencing/reads
- this script takes ~14h to finish -

output:
-------
/tophat_out
/cufflinks_out
snpRiboSeq.raw.bcf
snpRiboSeq.vcf = SNP file

SQLite database:
riboseqSNP.db
--> snps table
id | chromosome | position | reference base | alternative base


	-- done!






