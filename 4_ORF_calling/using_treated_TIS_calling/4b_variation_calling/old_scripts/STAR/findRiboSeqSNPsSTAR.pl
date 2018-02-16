#!/usr/bin/perl -w

#######################################################################################################
#
#	this script runs STAR, Cufflinks and Samtools on ribosomal sequencing reads to identify SNPs
#       the results are stored in an SQLite table
#
#	usage:
#	$ perl findRiboSeqSNPsSTAR.pl --riboreads path/to/sequencing/reads --genome path/to/genome/sequence --stargenome path/to/STAR/genome/annotation --cufflinks
#
#	example:
#	nohup perl ~/perlScripts/findRiboSeqSNPsSTAR.pl --riboreads ../SRR315061_norrna.fq --genome ../bowtie/musMusculus.fa --stargenome GenomeDir/ &
#
#       results:
#       in the working directory it creates the following output:
#        - GenomeDir/ (if it doesn't exist yet, contains genome annotation info)
#        - STAR output files
#        - cufflinks_out/ (folder)
#        - Samtools output (snpRiboSeq.raw.bcf & snpRiboSeq.vcf)
#        - SQLite database riboseqSNP.db with the SNPs table:
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
my ($riboSeqReads, $genomeSequence, $starGenomeDir, $cufflinks);
GetOptions(
    "riboreads=s"=>\$riboSeqReads, # path to the ribo-seq results, mandatory argument
    "genome=s"=>\$genomeSequence, # path to the fasta genome sequence, mandatory argument
    "stargenome:s"=>\$starGenomeDir, # path to the STAR genome annotation, optional argument (if the annotation doesn't exist, it's created by this script)
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
if ($starGenomeDir){
    print "STAR genome annotation directory: $starGenomeDir\n\n";
} else {
    # choose the default name for the STAR genome annotation directory
    $starGenomeDir = "GenomeDir/";
}


##################
#
# step 1: generate STAR genome
#
# check if a STAR reference genome folder exists, if it doesn't create it
#
#
##################

my $STARstart = time;
print "checking for STAR genome folder...\n";
if (-d $starGenomeDir){
    print "ok, folder found\n\n";
} else {
    print "no STAR genome directory found\ncreating STAR genome annotation...\n";
    system("mkdir GenomeDir/");
    system("STAR --runMode genomeGenerate --genomeDir GenomeDir/ --genomeFastaFiles $genomeSequence --runThreadN 12");
    systemError("STAR genome",$?,$!);
    print "done!\n\n";
}


##################
#
# step 2: run STAR alignment
#
#
##################

print "running STAR...\n";
# alignment
system("STAR --genomeDir $starGenomeDir --readFilesIn $riboSeqReads --runThreadN 12 --outFilterMismatchNmax 3");
systemError("STAR",$?,$!);
# convert SAM output file to BAM file
print "converting SAM output to BAM...\n";
system("samtools view -bS -o Aligned.out.bam Aligned.out.sam");
systemError("Samtools view",$?,$!);
# sort the BAM file
print "sorting STAR hits...\n";
system("samtools sort Aligned.out.bam aligned.sorted");
systemError("Samtools sort",$?,$!);
my $STARend = time - $STARstart;
print "STAR is finished\n";
printf("runtime: %02d:%02d:%02d\n\n",int($STARend/3600), int(($STARend % 3600)/60), int($STARend % 60));


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
# check if the cufflinks flag was passed to the script
if ($cufflinks){
    $cuffStart = time;
    print "running Cufflinks...\n";
    system("mkdir cufflinks_out");
    system("cufflinks --quiet --library-type fr-firststrand -o cufflinks_out aligned.sorted.bam");
    # the --library-type option is used to indicate that the RNA-seq data is stranded, without it cufflinks can't deal with the STAR output
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
print (" - running mpileup...\n");
system("samtools mpileup -uD -f $genomeSequence aligned.sorted.bam | bcftools view -bvcg - > snpRiboSeq.raw.bcf");
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
# step 4: create SQLite database
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
printf("total runtime: %02d:%02d:%02d\n",int($endTime/3600), int(($endTime % 3600)/60), int($endTime % 60));
printf(" --> STAR: %02d:%02d:%02d\n",int($STARend/3600), int(($STARend % 3600)/60), int($STARend % 60));
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

