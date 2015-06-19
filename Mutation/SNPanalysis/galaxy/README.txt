README SNP calling (alexander)

1¬∞ General prerequisites
------------------------

- bash (awk, date, sort, unique, touch, rep, printf)
- perl 5.12+
- SQLite3
- picard (http://sourceforge.net/projects/picard/files/picard-tools/) --> specify the location of the picard-tools directory at the top of the snp_calling bash script!
- samtools (http://sourceforge.net/projects/samtools/files/samtools/)
- igenome of organism (http://tophat.cbcb.umd.edu/igenomes.html --> /home/galaxy/data/igenomes/Mus_musculus/)
- Ensembl SQLite database (Ensembl Mus musculus 72 --> /home/galaxy/app/tools/ribo_prof_pipeline/tool-data/ENS_mmu_72.db)
- SNPdb data in VCF file format
	• mouse: ftp://ftp.ncbi.nih.gov/snp/organisms/mouse_10090/VCF/genotype/SC_MOUSE_GENOMES.genotype.vcf.gz --> /home/galaxy/app/tools/ribo_prof_pipeline/tool-data/snpdb_mouse_GRCm38.vcf
	• human: ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/common_all.vcf.gz --> /home/galaxy/app/tools/ribo_prof_pipeline/tool-data/snpdb_human_GRCm37.vcf
	• fruitfly: --> /home/galaxy/app/tools/ribo_prof_pipeline/tool-data/snpdb_fruitfly_BDGPr5.vcf
		the SNP data for drosophila is not available in one file, the following lines of code show you how to download the separate files and how to merge them in a single vcf file
		wget ftp://ftp.ncbi.nih.gov/snp/organisms/fruitfly_7227/VCF/vcf_chr_2L.vcf.gz
		wget ftp://ftp.ncbi.nih.gov/snp/organisms/fruitfly_7227/VCF/vcf_chr_2R.vcf.gz
		wget ftp://ftp.ncbi.nih.gov/snp/organisms/fruitfly_7227/VCF/vcf_chr_3L.vcf.gz
		wget ftp://ftp.ncbi.nih.gov/snp/organisms/fruitfly_7227/VCF/vcf_chr_3R.vcf.gz
		wget ftp://ftp.ncbi.nih.gov/snp/organisms/fruitfly_7227/VCF/vcf_chr_X.vcf.gz
		gunzip vcf_chr_*.vcf.gz
		mv vcf_chr_2L.vcf snpdb_fruitfly_BDGPr5.vcf
		grep -v '^#' vcf_chr_2R.vcf >> snpdb_fruitfly_BDGPr5.vcf
		grep -v '^#' vcf_chr_3L.vcf >> snpdb_fruitfly_BDGPr5.vcf
		grep -v '^#' vcf_chr_3R.vcf >> snpdb_fruitfly_BDGPr5.vcf
		grep -v '^#' vcf_chr_X.vcf >> snpdb_fruitfly_BDGPr5.vcf
		rm vcf_chr_2R.vcf vcf_chr_3L.vcf vcf_chr_3R.vcf vcf_chr_X.vcf

2¬∞ SNP calling on mapped reads
------------------------------

Prerequisites:

- the following perl scripts:
	• filterSAMfile.pl (/home/galaxy/app/tools/ribo_prof_pipeline/filterSAMfile.pl)
		--> extracts variants from mapped reads in a SAM file
	• snpIndexBuilder.pl (/home/galaxy/app/tools/ribo_prof_pipeline/snpIndexBuilder.pl)
		--> creates a unique index for every variant
	• splitVCFaltRecords.pl (/home/galaxy/app/tools/ribo_prof_pipeline/splitVCFaltRecords.pl)
		--> splits variants with multiple alternatives into separate records

- sam file with the mapped reads
- SQLite database to store the results, should also contain the table "arguments" from which the organism name can be obtained


The input variables of the SNP calling script are:

required variables:
-s/--sqlitedb = path to the SQLite database where the results need to be stored
-e/--ensembldb = path to the Ensembl database
--removeduplicates = "y" or "n", indicates whether or not duplicate reads have to be removed (uses picard)
--snpdb = "true" or "false", indicates whether or not the mapped reads will be searched for known SNPs data (Sanger institute)
-r/--reads = path to the sam file that contains the mapped reads

optional variables:
-t/--threads = number of threads to run the script on, the SNPs are called for each chromosome individually, so the for best performance this value has to be set to the number of chromosomes of the analyzed organism (default = 1)
--mincoverage = samtools parameter, the minimal number of reads that need to map at a location so that a SNP can be called there (default = 3)
--maxcoverage = samtools parameter, the maximal number of reads that can map at a location for a SNP to be called there (default = 100)
--high_af, --lower_af & --upper_af = high, lower and upper allelic frequency, input variables for mergeVCFfiles.pl, select SNPs & INDELS when their allelic frequency is between lower_af and upper_af or higher than high_af (default = 0.95, 0.3, 0.7)

Usage:

snp_calling_samtools --sqlitedb path/to/results/database.db --ensembldb path/to/ensembl/database.db --removeduplicates [y|n] --snpdb [true|false] --reads path/to/mapped/reads.sam --threads 22 --mincoverage 3 --maxcoverage 100 --high_af 0.95 --lower_af 0.3 --upper_af 0.7

