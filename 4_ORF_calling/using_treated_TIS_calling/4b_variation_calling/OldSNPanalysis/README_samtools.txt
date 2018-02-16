README SNP calling samtools (alexander)

1¬∞ General prerequisites
------------------------

- bash
- SQLite
- samtools (http://sourceforge.net/projects/samtools/files/samtools/)
- mergeVCFfiles.pl

2¬∞ SNP calling on mapped reads
------------------------------

Prerequisites:

- sam file with the mapped reads
- igenome of organism (http://tophat.cbcb.umd.edu/igenomes.html)

The input variables of the SNP calling script are:

required variables:
-e/--experiment = name of the experiment, is used to select the SQLite database to store the results
-r/--reads = path to the sam file that contains the mapped reads
-c/--chromosomes = path to the igenomes folder with the chromosome sequence files (fasta format)
-o/--organism = name of the organism: mouse, human or fruitfly

optional variables:
-t/--threads = number of threads to run the script on, the SNPs are called for each chromosome individually, so the for best performance this value has to be set to the number of chromosomes of the analyzed organism (default = 1)
--mincoverage = samtools parameter, the minimal number of reads that need to map at a location so that a SNP can be called there (default = 3)
--maxcoverage = samtools parameter, the maximal number of reads that can map at a location for a SNP to be called there (default = 100)
--high_af, --lower_af & --upper_af = high, lower and upper allelic frequency, input variables for mergeVCFfiles.pl, select SNPs & INDELS when their allelic frequency is between lower_af and upper_af or higher than high_af (default = 0.95, 0.3, 0.7)

Usage:

snp_calling_samtools -e experimentName -r pathToMappedReads -c pathToChromosomeSequencesFolder -o organism -t 1 --mincoverage 3 --maxcoverage 100 --high_af 0.95 --lower_af 0.3 --upper_af 0.7