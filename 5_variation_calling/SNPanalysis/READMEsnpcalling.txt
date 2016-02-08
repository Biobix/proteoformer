README SNP calling (alexander, adaptations 2015 steven)

1¬∞ General prerequisites
------------------------

- bash (awk, date, sort, unique, touch, rep, printf)
- perl 5.12+
- SQLite3
- picard (http://sourceforge.net/projects/picard/files/picard-tools/) --> the location of picard tools need to be specified in the toolsdir argument
- samtools (http://sourceforge.net/projects/samtools/files/samtools/)
- vcfutils.pl (https://github.com/lh3/samtools/blob/master/bcftools/vcfutils.pl)
- splitVCFaltRecords.pl, snpIndexBuilder,pl and filterSAMfile.pl
- igenome of organism (get_igenomes.py to download all necessary files)
- Ensembl SQLite database (ENS_db.py to download the right database)
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
	• filterSAMfile.pl
		--> extracts variants from mapped reads in a SAM file
	• snpIndexBuilder.pl
		--> creates a unique index for every variant
	• splitVCFaltRecords.pl
		--> splits variants with multiple alternatives into separate records
	These scripts need to be available in one folder, mentioned in the ——tooldir argument
- sam file with the mapped reads
- SQLite database to store the results, should also contain the table "arguments" from which the organism name can be obtained


The input variables of the SNP calling script are:

required variables:
-s/--sqlitein = path to the SQLite database where previous results of the pipeline were stored
--sqliteout = path to the SQLite database where the results need to be stored
-e/--ensembldb = path to the Ensembl database
--removeduplicates = "y" or "n", indicates whether or not duplicate reads have to be removed (uses picard)
--snpdbselected = “y” or “n”, indicates whether or not the mapped reads will be searched for known SNPs data (Sanger institute)
-r/--reads = path to the sam file that contains the mapped reads
—-toolsdir = path to the folder where filterSAMfile.pl, snpIndexBuilder.pl and splitVCFaltRecords.pl are

optional variables:
—-galaxydir = galaxy specific parameter
—-picardpath = path to map where picard jar files are (mandatory if you want to remove duplicates)
—-snpdb = path to the snpDB (mandatory if you want to search in known SNP data)
--mincoverage = samtools parameter, the minimal number of reads that need to map at a location so that a SNP can be called there (default = 3)
--maxcoverage = samtools parameter, the maximal number of reads that can map at a location for a SNP to be called there (default = 100)
--high_af, --lower_af & --upper_af = high, lower and upper allelic frequency, input variables for mergeVCFfiles.pl, select SNPs & INDELS when their allelic frequency is between lower_af and upper_af or higher than high_af (default = 0.95, 0.3, 0.7)

Usage:

bash snp_calling --sqlitein path/to/results/database.db —-sqliteout path/to/output/database.db --removeduplicates [y|n] ——picardpath /path/to/picardmap —-snpdbselected [y|n] --snpdb --path/to/snpdb tooldir /path/to/tooldir --reads path/to/mapped/reads.sam --mincoverage 3 --maxcoverage 100 --high_af 0.95 --lower_af 0.3 --upper_af 0.7

