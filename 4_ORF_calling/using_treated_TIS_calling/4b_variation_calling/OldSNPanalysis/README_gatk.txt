README SNP calling gatk (alexander)

1¬∞ General prerequisites
------------------------

- bash
- SQLite
- GATK (http://www.broadinstitute.org/gatk/download)
- picard (http://sourceforge.net/projects/picard/files/picard-tools/)
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

optional variables:
--high_af, --lower_af & --upper_af = high, lower and upper allelic frequency, input variables for mergeVCFfiles.pl, select SNPs & INDELS when their allelic frequency is between lower_af and upper_af or higher than high_af (default = 0.95, 0.3, 0.7)
--dmq = GATK parameter, default mapping quality, used to reset the mapping quality scores, which is necessary to avoid filtering out all ribo-seq reads (default = 60)
--stand_call_conf = GATK parameter, minimum phred-scaled confidence threshold at which variants should be called (default = 30.0)
--stand_emit_conf = GATK parameter, minimum phred-scaled confidence threshold at which variants should be emitted (default = 10.0)
--dcov = GATK parameter, downsamples the coverage at any location to the value provided, e.g. if dcov = 200 and the amount of reads that mapped at a certain location = 400, GATK will randomly select 200 reads from these 400 for the SNP calling (default = 500)
--min_mapping_quality = GATK parameter, mapping quality threshold to use a read for the SNP calling (default = 20)

Usage:

snp_calling_gatk -e experimentName -r pathToMappedReads -g pathToGATKdirectory -p pathToPicardDirectory -s pathToReferenceSequence --high_af 0.95 --lower_af 0.3 --upper_af 0.7 --dmq 60 --stand_call_conf 30.0 --stand_emit_conf 10.0 --dcov 500 --min_mapping_quality 20