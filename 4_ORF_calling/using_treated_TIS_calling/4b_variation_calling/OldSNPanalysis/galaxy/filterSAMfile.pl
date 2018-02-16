#!/usr/bin/perl -w
#
# Extract mismatches (SNPs and INDELs), aka alternative bases, and their corresponding reference bases from mapped reads in a SAM file.
# Mismatches are found based on the CIGAR string and MD tag (for more info: http://samtools.sourceforge.net/SAMv1.pdf):
#   - does the MD tag contain an A, T, C or G?
#       yes: there are one or more mismatches
#       no: no mismatches, go to the next read
#   - find the position of the mismatch in the read based on the numbers of matching bases that precede the mismatch
#     using the number in the MD tag right before the A, T, C or G
#   - check the CIGAR string for clipped bases (indicated with an S) and calculate the genomic position of the mismatch
# Input = SAM file
# Output = mismatches.txt file with the following columns:
#   chromosome name | position | reference base | alternative base(s)
#
# Remark: the output will contain duplicates, so don't forget to remove those (eg sort mismatches.txt | uniq > unique_mismatches.txt).
#

use strict;

my $startTime = time;

my $usage = "\nUsage: $0 <infile.sam>\n\n";
my $samFile = shift or die $usage;

open(SAM, "<$samFile") or die "Couldn't open the sam file $samFile\n".$!."\n";
open(SNP, ">>mismatches.txt") or die "Couldn't open mismatch file\n".$!."\n";

print "Reading sam file...\n";

while (<SAM>){
    
    chomp;
    my $samFileLine = $_;
    
    if ($samFileLine !~ /^@/){
        
        my @lineArray = split('\t', $samFileLine);
        # the columns of interest are:
        # column 3: chromosome
        # column 4: position
        # column 6: CIGAR string
        # column 10: read sequence
        # column 19: MD tag
        my $chr = $lineArray[2];
        my $pos = $lineArray[3];
        my $cigar = $lineArray[5];
        my $read = $lineArray[9];
        my $md = $lineArray[-1]; # the number of columns depends on whether duplicate reads were removed with picard or not, but the MD tag will be the last column either way
        $md =~ s/MD:Z://;
        
        my $validCigar = 1;
        if ($md =~ /[ATCG]/ and $cigar !~ /[ID]/){ # is there a SNP in the MD tag? INDELS are also not considered atm
            
            # assign a genomic reference position to each base of the mapped read
            my @genomicPositions;
            my $readPos = $pos;
            my $cigarCopy = $cigar;
            while ($cigarCopy !~ /^$/){
                if ($cigarCopy =~ /^([0-9]+[MNIDS])/){
                    my $cigar_part = $1;
                    if ($cigar_part =~ /(\d+)M/){
                        for (my $t = 0; $t < $1; $t++){
                            push(@genomicPositions, $readPos);
                            $readPos++;
                        }
                    # this part will have to deal with INDELS:
                    #} elsif ($cigar_part =~ /(\d+)I/){
                    #    for (my $t = 0; $t < $1; $t++){
                    #        push(@genomicPositions, "I");
                    #    }
                    #} elsif ($cigar_part =~ /(\d+)D/){
                    #    $readPos += $1;
                    } elsif ($cigar_part =~ /(\d+)S/){
                        for (my $t = 0; $t < $1; $t++){
                            push(@genomicPositions, "S");
                        }
                    } elsif ($cigar_part =~ /(\d+)N/){
                        $readPos += $1;
                    }
                    $cigarCopy =~ s/$cigar_part//;
                } else {
                    $validCigar = 0;
                    last;
                }
            }
            
            if ($validCigar){
                
                # go through the MD tag information character by character
                # numbers indicate the number of matching bases (beware: a number can consist of more than one digit!)
                # a letter indicates a SNP (the letter in the MD tag = the reference base, the corresponding base in the read sequence = the alternate base)
                # letters that follow a caret sign ("^") indicate deleted bases
                my $nrMatches = "";
                my $deletion = "";
                my $readPosition;
                # clipped bases are not represented in the MD tag
                # so if we don't account for possible clipped bases at the start of the read, the positions will be off
                if ($cigar =~ /^([0-9]+)S/){
                    $readPosition = $1;
                } else {
                    $readPosition = 0;
                }
                
                for (my $t = 0; $t < length($md); $t++){
                    
                    my $character = substr $md, $t, 1;
                    
                    if ($character =~ /[0-9]/){
                        $nrMatches .= $character;
                        $deletion = "";
                    } elsif ($character =~ /[ACTG]/){
                        if ($deletion eq ""){
                            $readPosition += $nrMatches;
                            $nrMatches = "";
                            my $snpPos = $genomicPositions[$readPosition];
                            my $alt = substr $read, $readPosition, 1;
                            if ($snpPos ne "S"){
                                while ($snpPos eq "I"){
                                    $readPosition++;
                                    $snpPos = $genomicPositions[$readPosition];
                                }
                                print SNP "$chr;$snpPos;$character;$alt\n";
                            }
                            $readPosition++; # this is to account for the SNP position itself
                        } else {
                            $deletion .= $character;
                        }
                    } elsif ($character eq "^"){ # deletion from the reference chromosome
                        $deletion .= "^";
                        $readPosition += $nrMatches;
                        $nrMatches = "";
                    } elsif ($character eq "N"){ # ignore SNPs where the alternative base = "N"
                        if ($nrMatches ne ""){
                            $readPosition += $nrMatches;
                        }
                        $readPosition++;
                        $nrMatches = "";
                    }
                    
                }
            }
            
        }
        
    }
    
}

close SAM;
close SNP;

my $runTime = time - $startTime;
printf("runtime: %02d:%02d:%02d", $runTime/3600, ($runTime % 3600)/60, ($runTime % 3600) % 60);
print "\n";