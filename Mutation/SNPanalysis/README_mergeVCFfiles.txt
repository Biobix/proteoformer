README merge vcf files (alexander)

1¬∞ General prerequisites
------------------------

- perl 5+

2¬∞ Merge all vcf files in a directory
-------------------------------------

Prerequisites:

- ./tmp/ directory with one or more .vcf files

The input variables of the script are:

my $cutoff1 = $ARGV[0]; # this is the high_af value
my $cutoff2 = $ARGV[1]; # this is the lower_af value
my $cutoff3 = $ARGV[2]; # this is the upper_af value
my $resultsFile = "./tmp/allSNPsID.txt";

Usage:

This script is called by both snp_calling_samtools and snp_calling_gatk to combine the vcf files they produce in one file.
mergeVCFfiles.pl high_af lower_af upper_af