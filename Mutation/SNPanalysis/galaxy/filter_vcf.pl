#!/usr/bin/perl -w

use strict;

my $high_af = $ARGV[0];
my $lower_af = $ARGV[1];
my $upper_af = $ARGV[2];
my $snpFile = $ARGV[3];
my $resultsFile = $ARGV[4];

# open the  VCF file
# remove the annotation lines (those that start with #),
# select the appropriate columns and select the SNPs where the allelic frequency is between cutoff2 and cutoff3 or greater than cutoff1
# write the result to the output text file
open (FILTERED, ">>$resultsFile") or die "couldn't open $resultsFile\n$!\n\n";
my $id = 1;

open (VCF, $snpFile) or die "couldn't open $snpFile\n$!\n\n";
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
	    if ($allelicFrequency > $high_af || ($allelicFrequency > $lower_af && $allelicFrequency < $upper_af)){
		print FILTERED "$id;$snpChr;$snpPos;$snpRef;$snpAlt;$allelicFrequency\n";
	        $id++;
	    }
	}
    }
}
close VCF;
close FILTERED;
