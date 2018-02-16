README SNP calling gatk & samtools (alexander)

1¬∞ General prerequisites
------------------------

- bash
- SQLite
- GATK (http://www.broadinstitute.org/gatk/download)
- picard (http://sourceforge.net/projects/picard/files/picard-tools/)
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
-g/--gatk_dir = path to the GATK directory, example: /home/user/GenomeAnalysisTK-2.7-2
-p/--picard_dir = path to the picard toolkit directory, example: /home/user/picard-tools-1.97
-s/--sequence = path to the genome reference genome sequence (fasta file)
-c/--chromosomes = path to the igenomes folder with the chromosome sequence files (fasta format)
-o/--organism = name of the organism: mouse, human or fruitfly

optional variables:
-t/--threads = number of threads to run the script on, the SNPs are called for each chromosome individually, so the for best performance this value has to be set to the number of chromosomes of the analyzed organism (default = 1)
--mincoverage = samtools parameter, the minimal number of reads that need to map at a location so that a SNP can be called there (default = 3)
--maxcoverage = samtools parameter, the maximal number of reads that can map at a location for a SNP to be called there (default = 100)
--high_af, --lower_af & --upper_af = high, lower and upper allelic frequency, input variables for mergeVCFfiles.pl, select SNPs & INDELS when their allelic frequency is between lower_af and upper_af or higher than high_af (default = 0.95, 0.3, 0.7)

Usage:

snp_calling_gatk_samtools -g pathToGATKdirectory -p pathToPicardDirectory -e experimentName -r pathToMappedReads -s pathToReferenceSequence -c pathToChromosomeSequencesFolder -o organism -t 1 --mincoverage 3 --maxcoverage 100 --high_af 0.95 --lower_af 0.3 --upper_af 0.7

Difference between this script and snp_calling_gatk & snp_calling_samtools:

This script combines picard, GATK and samtools for the SNP calling in mapped reads. It starts off in the same way as snp_calling_gatk by sorting the reads, removing the duplicates and indexing the reads with picard. Next, it uses the GATK toolkit to realign the mapped reads around INDELS. These realigned reads are then fed to samtools for the SNP calling.