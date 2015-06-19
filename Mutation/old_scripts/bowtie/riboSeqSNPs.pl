#!/usr/bin/perl -w

#######################################################################################################
#
#	this script runs TopHat, Cufflinks and Samtools on ribosomal sequencing reads to identify SNPs
#       the results are stored in an SQLite table
#
#	usage:
#	$ perl findRiboSeqSNPs.pl --genome path/to/genome/sequence --index path/to/bowtie/index --riboreads path/to/sequencing/reads --cufflinks
#
#       result:
#       in the working directory it creates the following output:
#        - tophat_out/ (folder)
#        - cufflinks_out/ (folder)
#        - Samtools output (snpRiboSeq.raw.bcf & snpRiboSeq.vcf)
#        - SQLite table (in database riboseqSNP.db) with the SNPs:
#          snps
#          id | chromosome | position | reference base | alternative base
#
#
#
#######################################################################################################


use strict;
use Getopt::Long;
use DBI;

my $startTime = time;

print "\n-----------------------------";
print "\n Ribo-SEQ SNP analysis\n";
print " finding the SNPs\n";
print "-----------------------------\n\n\n";

# get the command line arguments
my ($riboSeqReads, $genomeSequence, $bowtieIndex, $cufflinks);
GetOptions(
    "riboreads=s"=>\$riboSeqReads, # path to the ribo-seq results, mandatory argument
    "genome=s"=>\$genomeSequence, # path to the fasta genome sequence, mandatory argument
    "index=s"=>\$bowtieIndex, # path to the bowtie index, mandatory argument
    "cufflinks"=>\$cufflinks # cufflinks flag to indicate wheather or not the cufflinks analysis should be performed
);

if ($riboSeqReads){
    print "Ribo-Seq data: $riboSeqReads\n";
} else {
    die "\nDon't forget to pass the location of the Ribo-Seq fastq file using the --riboreads or -r argument!\n\n";
}
if ($genomeSequence){
    print "Genome fasta file: $genomeSequence\n";
} else {
    die "\nDon't forget to pass the location of the genome fasta file using the --genome or -g argument!\n\n";
}
if ($bowtieIndex){
    print "Bowtie index: $bowtieIndex\n";
} else {
    die "\nDon't forget to pass the location of the bowtie index files using the --index or -i argument!\n\n";
}


##################
#
# step 1: TopHat
# 
# aligns the RiboSeq reads to the genome using bowtie
# and analyzes the mapping results to find splice junctions between exons
#
#
##################

my $tophatStart = time;
print "running TopHat ...\n\n";
# tophat -N 3 -o tophat/tophat_out -p 2 --segment-mismatches 3 /data/bowtie/indexes/Mus_musculus/Ensembl/NCBIM37/Sequence/BowtieIndex/genome SRR315061_norrna.fq
system("tophat -N 3 -o tophat_out -p 2 $bowtieIndex $riboSeqReads"); # --segment-mismatches 3
systemError("TopHat",$?,$!); # $? = exit status (-1 if failed) and $! = error reason
my $tophatEnd = time - $tophatStart;
print "tophat is finished\n";
printf("runtime: %02d:%02d:%02d\n\n",int($tophatEnd/3600), int(($tophatEnd % 3600)/60), int($tophatEnd % 60));


#####################
#
# step 2: Cufflinks
#
# assembles the aligned RiboSeq reads into transcripts
#
#
#####################

my $cuffStart;
my $cuffEnd;
if ($cufflinks){
    $cuffStart = time;
    print "\nrunning Cufflinks ...\n\n";
    system("mkdir cufflinks_out");
    system("cufflinks --quiet -o cufflinks_out tophat_out/accepted_hits.bam");
    systemError("Cufflinks",$?,$!);
    $cuffEnd = time - $cuffStart;
    print "cufflinks is finished\n";
    printf("runtime: %02d:%02d:%02d\n\n",int($cuffEnd/3600), int(($cuffEnd % 3600)/60), int($cuffEnd % 60));
}


#####################
#
# step 3: Samtools
#
# sort: sorts the alignments by leftmost coordinates
# mpileup: calculate genotype likelihoods
# bcftools: SNP calling based on these likelihoods & converting bcf format (binary) to vcf (plain text)
# vcfutils.pl: filter out SNPs with read depth > 100
#
#
#####################

my $samStart = time;
print "running Samtools ...\n";
print " - sorting TopHat hits...\n";
system("samtools sort tophat_out/accepted_hits.bam accepted_hits.sorted");
systemError("Samtools sort",$?,$!);
print (" - running mpileup...\n");
system("samtools mpileup -uD -f $genomeSequence accepted_hits.sorted.bam | bcftools view -bvcg - > snpRiboSeq.raw.bcf");
systemError("Samtools mpileup",$?,$!);
print(" - running vcfutils.pl...\n");
my $snpFile = "snpRiboSeq.vcf";
system("bcftools view snpRiboSeq.raw.bcf | vcfutils.pl varFilter -D100 > $snpFile");
systemError("vcfutils.pl",$?,$!);
my $samEnd = time - $samStart;
print "Samtools is finished\n";
printf("runtime: %02d:%02d:%02d\n\n",int($samEnd/3600), int(($samEnd % 3600)/60), int($samEnd % 60));


#####################
#
# step 4: create SQLite table
#
#
#####################

print "\n-----------------------------";
print "\n Ribo-SEQ SNP analysis\n";
print " SNP SQLite table creation\n";
print "-----------------------------\n\n\n";

# create the database
my $sqliteStart = time;
print "creating SNP database...";
system("sqlite3 riboseqSNP .quit");
# connect to the database
my $dbh = DBI -> connect(
    "dbi:SQLite:dbname=riboseqSNP.db",
    "",
    "",
    { RaiseError => 1 },
) or die $DBI::errstr;
# create a table to save the SNP data
$dbh -> do("DROP TABLE IF EXISTS snps");
$dbh -> do("CREATE TABLE snps(id INT PRIMARY KEY, chr VARCHAR, pos INT, ref CHAR, alt CHAR)");
print "done!\n\nloading SNP data from $snpFile ...";

open(SNP,$snpFile) || die "couldn't open SNP file $snpFile\n$!\n\n";
my $snpCount = 1;
while(<SNP>){
    if ($_ !~ /^#/){ # skip the information lines
        chomp;
        my @line = split("\t",$_);
	my $snpChr = $line[0];
        my $snpPos = $line[1];
	my $snpRef = $line[3];
	my $snpAlt = $line[4];
	$dbh -> do("INSERT INTO snps VALUES ('".$snpCount."','".$snpChr."','".$snpPos."','".$snpRef."','".$snpAlt."')");
	$snpCount++;
    }
}
close SNP;
$dbh -> disconnect();
my $sqliteEnd = time - $sqliteStart;
print "SQLite is finished\n";
printf("total runtime: %02d:%02d:%02d\n\n",int($sqliteEnd/3600), int(($sqliteEnd % 3600)/60), int($sqliteEnd % 60));


print "\n---------------------------------------\n";
$snpCount--;
print "$snpCount SNPs were detected & loaded in the SQLite database riboseqSNP.db\n";
my $endTime = time - $startTime;
printf("runtime: %02d:%02d:%02d\n",int($endTime/3600), int(($endTime % 3600)/60), int($endTime % 60));
printf(" --> tophat: %02d:%02d:%02d\n",int($tophatEnd/3600), int(($tophatEnd % 3600)/60), int($tophatEnd % 60));
if ($cufflinks) {printf(" --> cufflinks: %02d:%02d:%02d\n",int($cuffEnd/3600), int(($cuffEnd % 3600)/60), int($cuffEnd % 60));}
printf(" --> samtools: %02d:%02d:%02d\n",int($samEnd/3600), int(($samEnd % 3600)/60), int($samEnd % 60));
printf(" --> SQLite: %02d:%02d:%02d\n",int($sqliteEnd/3600), int(($sqliteEnd % 3600)/60), int($sqliteEnd % 60));
print "---------------------------------------\n\n";


################################################################################
################################################################################
#
#   subroutines
#

sub systemError
{
    my ($command,$returnValue,$errorMessage) = @_;
    if ($returnValue == -1){
        die "$command failed!\n$errorMessage\n\n";
    }
}

