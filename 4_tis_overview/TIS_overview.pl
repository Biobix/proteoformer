#!/usr/bin/perl -w
use strict;
use warnings;
use DBI;
use DBD::SQLite;
use Getopt::Long;
use Cwd;

#################
# Perl commands #
#################

## Command line ##

# ./TIS_overview.pl --sqlite /home/galaxy/data/SQLite/results_STAR8.db

# Galaxy
#--out_tis_overview tis_overview --dir getcwd

# get the command line arguments
my ($work_dir,$tis_overview,$sqlite_db);
GetOptions(
"dir=s"=>\$work_dir,
"sqlite_db:s" =>\$sqlite_db,
"out_tis_overview:s" =>\$tis_overview
);

my $CWD             = getcwd;
my $TMP             =  "$CWD/tmp" ;

#Check if tmpfolder exists, if not create it...
if (!-d "$TMP") {
    system ("mkdir ". $TMP);
}

print "The following tmpfolder is used                          : $TMP\n";
print "The following results db is used                         : $sqlite_db\n";

# comment on these
if ($work_dir){
    print "The following working directory is used                  : $work_dir\n";
} else {
    $work_dir = $CWD;
    print "The following working directory is used                  : $CWD\n";
    #die "\nDon't forget to pass the working directory using the --dir or -d argument!\n\n";
}

if (!defined($tis_overview))     {$tis_overview        = $work_dir."/SQLite/overview.tis";}

# DB settings
# Sqlite Riboseq
my $db_results  = $sqlite_db;
my $dsn_results = "DBI:SQLite:dbname=$db_results";
my $us_results  = "";
my $pw_results  = "";

############################################################
## Creating TAB separated file based on TIS_overview table #
##                                                         #
############################################################

print "\n Create TAB seperated file based on TIS_overview table...";
## Create TIS_overview file
create_tab_file($dsn_results,$us_results,$pw_results,$tis_overview);
print "\nDone";

############
# THE SUBS #
############

### Get Analysis ID

sub create_tab_file{
    
    # Catch
    my $dsn             =   $_[0];
    my $us              =   $_[1];
    my $pw              =   $_[2];
    my $tis_overview    =   $_[3];
    
    # Init
    my $dbh = dbh($dsn,$us,$pw);
    my $temp_tis_overview = $TMP."/overview.tis";
    open TempTIS, "+>> ".$temp_tis_overview or die $!;
    
    # Create table 
    my $query = "select * from TIS_overview";
    my $sth = $dbh->prepare($query);
	$sth->execute();
    
    #Print to TIS_overview
    while (my $row = $sth->fetchrow_hashref) {
        print TempTIS $row->{"ID"}."\t".
        $row->{"local_max"}."\t".
        $row->{"min_count_aTIS"}."\t".
        $row->{"R_aTis"}."\t".
        $row->{"min_count_5UTR"}."\t".
        $row->{"R_5UTR"}."\t".
        $row->{"min_count_CDS"}."\t".
        $row->{"R_CDS"}."\t".
        $row->{"min_count_3UTR"}."\t".
        $row->{"R_3UTR"}."\t".
        $row->{"min_count_ntr"}."\t".
        $row->{"R_ntr"}."\t".
        $row->{"SNP"}."\t".
        $row->{"filter"}."\n";
    }
    system("mv ".$temp_tis_overview." ".$tis_overview);
    #system("rm -rf ".$temp_tis_overview);
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
