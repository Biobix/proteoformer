#!/usr/bin/perl -w

use strict;

my @input = @ARGV;
my $cutoff1 = shift(@input);
my $cutoff2 = shift(@input);
my $cutoff3 = shift(@input);
# the remaining input variables (1 or more) are the VCF files
my @snpFiles = @input;
my $resultsFile = "allSNPsID.txt";

# open the first VCF file
# remove the annotation lines (those that start with #),
# select the appropriate columns and select the SNPs where the allelic frequency is between cutoff2 and cutoff3 or greater than cutoff1
# write the result to the output text file
# move to the next VCF file
open (MERGED, ">>$resultsFile") or die "couldn't open $resultsFile\n$!\n\n";
my $id = 1;
foreach my $snpF (@snpFiles){
    open (VCF, $snpF) or die "couldn't open $snpF\n$!\n\n";
    while (<VCF>){
        if ($_ !~ /^#/){ # skip the information lines
            chomp;
            my @line = split("\t",$_);
            my $snpChr = $line[0];
            my $snpPos = $line[1];
	    my $snpRef = $line[3];
	    my $snpAlt = $line[4];
	    my $info = $line[7];
	    if ($info =~ /AF1?=(.*?);/){
		my $allelicFrequency = $1;
		$allelicFrequency =~ s/,.+$//;
		if ($allelicFrequency > $cutoff1 || ($allelicFrequency > $cutoff2 && $allelicFrequency < $cutoff3)){
		    print MERGED "$id;$snpChr;$snpPos;$snpRef;$snpAlt;$allelicFrequency\n";
		    $id++;
		}
	    }
        }
    }
    close VCF;
}
close MERGED;
