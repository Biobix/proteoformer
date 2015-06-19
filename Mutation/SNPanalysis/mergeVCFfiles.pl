#!/usr/bin/perl -w

use strict;

my $cutoff1 = $ARGV[0];
my $cutoff2 = $ARGV[1];
my $cutoff3 = $ARGV[2];
my $resultsFile = "./tmp/allSNPsID.txt";

# start by merging all the files together, remove the annotation lines (start with #),
# select the appropriate columns and select the SNPs where the allelic frequency is between cutoff2 and cutoff3 or greater than cutoff1
my @snpFiles = glob("./tmp/*.vcf");


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
