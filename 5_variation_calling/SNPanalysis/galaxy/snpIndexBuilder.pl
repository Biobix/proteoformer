#!/usr/bin/perl -w
#
# Create a unique index for a SNP, using the chromosome, location, reference base and alternative base.
# Input = text file with at least 4 columns: chromosome name, position, reference base, alternative base
# Extra columns will not influence the script, but the first four need to be in this order.
# The index is built as follows:
#   - convert the chromosome name to its numeric unicode value (chr) and add 90 to those values smaller than 10,
#     this will ensure that the indexes for eg chr1 and chr10 can't be the same by accident
#   - the position value can be used as it is (pos)
#   - convert every character in the 'alt' column to its numeric unicode value and sum them all up (altSum)
# Combine these three numeric values by concatenating them: chr.pos.altSum
# This will create a unique integer for every SNP or INDEL.
#

use strict;
use utf8;

my $usage = "\nUsage: $0 <input file>\n\n";
my $input = shift or die $usage;
my $output = $input;
$output =~ s/\.txt$/_indexed\.txt/;

open IN, "<$input" or die $!."\n";
open OUT, ">>$output" or die $!."\n";

#my $maxIndex = 0;

while (<IN>){
    
    chomp;
    
    if ($_ =~ /^$/){
        next;
    }
    
    my @columns = split(/;/, $_);
    my $pos = $columns[1];
    if (!$pos){ # skip the line if the position is undefined
        next;
    }
    if ($pos =~ /[a-zA-Z]/){
        next; # skip the line if the position column contains one or more characters instead of only numbers
    }
    my $posLength = length($pos);
    my $altNumber;
    my $altString = $columns[3];
    
    if (!$altString){ # skip the line if the alternative base is undefined
        next;
    }
    
    while (length($altString) > 0){
        
        if ($altString =~ /^(A+)/){
            if (length($1) > 4){
                $altNumber .= length($1)."1";
            } else {
                $altNumber .= "1"x(length($1));
            }
            $altString =~ s/^$1//;
        } elsif ($altString =~ /^(T+)/){
            if (length($1) > 4){
                $altNumber .= length($1)."2";
            } else {
                $altNumber .= "2"x(length($1));
            }
            $altString =~ s/^$1//;
        } elsif ($altString =~ /^(C+)/){
            if (length($1) > 4){
                $altNumber .= length($1)."3";
            } else {
                $altNumber .= "3"x(length($1));
            }
            $altString =~ s/^$1//;
        } elsif ($altString =~ /^(G+)/){
            if (length($1) > 4){
                $altNumber .= length($1)."4";
            } else {
                $altNumber .= "4"x(length($1));
            }
            $altString =~ s/^$1//;
        } elsif ($altString =~ /^([^ATCG]+)/) {
            if (length($1) > 4){
                $altNumber .= length($1)."0";
            } else {
                $altNumber .= "0"x(length($1));
            }
            $altString =~ s/^$1//;
        }
        
    }
    
    # transform the chromosome names to numeric values
    my $chr = $columns[0];
    if ($chr eq "X"){
        $chr = 80;
    } elsif ($chr eq "Y"){
        $chr = 81;
    } elsif ($chr eq "MT"){
        $chr = 82;
    } elsif ($chr =~ /(\d)([RL])/i){
        if ($2 eq "R"){
            $chr = $1 + 10;
        } elsif ($2 eq "L"){
            $chr = $1 + 20;
        }
    } elsif ($chr < 10){ # avoid potential overlap between eg chr1 and chr10 by adding 90 to those chr < 10
        $chr += 90;
    }
    # combine the different numbers to create a unique index
    my $numLeft = 9 - $posLength;
    my $index = $chr."0"x$numLeft.$pos.$altNumber;
    print OUT $index;
    foreach my $c (@columns){
        print OUT ";".$c;
    }
    print OUT "\n";
    
    #if ($index > $maxIndex){
    #    $maxIndex = $index;
    #}
    
}

#print "maximal index = $maxIndex\n";
#print "maximal SQLite INTEGER size = 9223372036854775808\n";
#print "index - max size = ".(9223372036854775808 - $maxIndex)."\n";

close IN;
close OUT;