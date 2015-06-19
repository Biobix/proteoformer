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


my $table = "TIS_1_transcripts";
my $sqlite_db = "SQLite/results_STARclip_26-35_1-14MM.db";
my $transcripts = {};
my $IGENOMES_ROOT = "/storage/igenomes/Mus_musculus/Ensembl/GRCm38/";
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

# Get R values
#my $query = "SELECT tr_stable_id||'_'||start as id_start,chr,strand,start from ".$table." ";
my $query = "SELECT tr_stable_id||'_'||start as id_start,Rltm_min_Rchx from ".$table." ";
my $dbh = dbh($dsn_results,$us_results,$pw_results);
my $sth = $dbh->prepare($query);
$sth->execute();
$transcripts = $sth->fetchall_hashref('id_start');
my $R;

### Output fasta
open CSV, "+>> R.csv " or die $!;

foreach my $tr (keys %{$transcripts}){
    
    if ($transcripts->{$tr}{'Rltm_min_Rchx'} eq 'NA'){
        next;
    }else{
        $R = sprintf("%.2f", $transcripts->{$tr}{'Rltm_min_Rchx'});
    print CSV $R."\n";
    }
}

close (CSV);

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
