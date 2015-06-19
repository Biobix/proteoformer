#!/usr/bin/perl -w
use strict;
use warnings;
use DBI;
use DBD::SQLite;
use Data::Dumper;
use Getopt::Long;
use v5.10;
use Parallel::ForkManager;
use Storable 'dclone';
use Cwd;


my $sorf_table = "TIS_sORFs_2_transcripts";
my $sqlite_db = "SQLite/results_tophat.db";
my $sorfs = {};
my $IGENOMES_ROOT = "/storage/igenomes/Drosophila_melanogaster/Ensembl/BDGP5/";
my $BIN_chrom_dir = $IGENOMES_ROOT."/Sequence/Chromosomes_BIN";
my $strand;
my $web_logo_begin;
my $web_logo_end;
my $web_logo_seq;

# DB settings
# Sqlite Riboseq
my $db_results  = $sqlite_db;
my $dsn_results = "DBI:SQLite:dbname=$db_results";
my $us_results  = "";
my $pw_results  = "";

# Get transcripts
#my $query = "SELECT id||'_'||sorf_begin as sorf_id,chr,strand,sorf_begin, sorf_end from ".$sorf_table." where annotation != 'aTIS' and annotation != 'exonic' ";
my $query = "SELECT id||'_'||sorf_begin as sorf_id,chr,strand,sorf_begin, sorf_end from ".$sorf_table." ";
my $dbh = dbh($dsn_results,$us_results,$pw_results);
my $sth = $dbh->prepare($query);
$sth->execute();
$sorfs = $sth->fetchall_hashref('sorf_id');

### Output fasta
open FASTA, "+>> sorf_weblogo.fa " or die $!;

foreach my $sorf (keys %{$sorfs}){
    
    $strand = $sorfs->{$sorf}{'strand'};
    $web_logo_begin = ($strand eq "1") ? $sorfs->{$sorf}{'sorf_begin'} - 3 : $sorfs->{$sorf}{'sorf_end'} - 3 ;
    $web_logo_end = ($strand eq "1") ? $sorfs->{$sorf}{'sorf_begin'} + 3 : $sorfs->{$sorf}{'sorf_end'} + 3 ;
    
    $web_logo_seq = ($strand eq "1") ? get_sequence($sorfs->{$sorf}{'chr'},$web_logo_begin,$web_logo_end) : revdnacomp(get_sequence($sorfs->{$sorf}{'chr'},$web_logo_begin,$web_logo_end));
    
    print FASTA ">".$sorfs->{$sorf}{'sorf_id'}."\n";
    print FASTA $web_logo_seq."\n";
}

close (FASTA);

### get reverse complement sequence ###

sub revdnacomp {
    my $dna = shift;
    my $revcomp = reverse($dna);
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}


### get sequence from binary chromosomes ###

sub get_sequence {
    
    #Catch
    my $chr   =     $_[0];
    my $start =     $_[1];
    my $end   =     $_[2];
    
    open(IN, "< ".$BIN_chrom_dir."/".$chr.".fa");
    binmode(IN);
    
    # Always forward strand sequence
    my $buffer; # Buffer to get binary data
    my $length = $end - $start + 1; # Length of string to read
    my $offset = $start-1; # Offset where to start reading
    
    seek(IN,$offset,0);
    read(IN,$buffer,$length);
    close(IN);
    
    # Return
    return($buffer);
}

### DBH ###

sub dbh {
    
    # Catch
    my $db  =   $_[0];
    my $us	=   $_[1];
    my $pw	=   $_[2];
    
    # Init DB
    my $dbh = DBI->connect($db,$us,$pw,{ RaiseError => 1 },) || die "Cannot connect: " . $DBI::errstr;
    
    #Return
    return($dbh);
}
