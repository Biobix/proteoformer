#!/usr/bin/perl -w
#
# Split the 'alt' column in processed SNPdb files or SNP calling result files whenever there is more than one alternative base.
# When there are multiple alternative bases, they are separated by a comma so we'll just look for that.
# Input = text file with 4 or 6 columns: chromosome name, position, reference base, alternative base (+ read depth & allelic frequency)
#

use strict;

my $usage = "\nUsage: $0 <infile.txt>\n\n";
my $input = shift or die $usage;

open(IN, "<$input") or die $!;
open(OUT, ">>split.txt") or die $!;

while (<IN>){
    chomp;
    my @line = split(/;/, $_);
    if ($line[3] =~ /,/){ # check whether the 'alt' column contains a comma
        my @alt = split(/,/, $line[3]);
        my $lineLength = @line;
        if ($lineLength == 4){
            foreach my $altBase (@alt){
                print OUT $line[0].";".$line[1].";".$line[2].";".$altBase."\n";
            }
        } elsif ($lineLength == 6){ # there could be an extra column for the read depth and for the allelic frequency
            foreach my $altBase (@alt){
                print OUT $line[0].";".$line[1].";".$line[2].";".$altBase.";".$line[4].";".$line[5]."\n";
            }
        }
    } else {
        print OUT $_."\n";
    }
}

close(IN);
close(OUT);