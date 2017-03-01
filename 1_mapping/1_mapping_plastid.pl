#!/usr/bin/perl -w

$|=1;

use strict;
use warnings;
use DBI;
use Data::Dumper;
use Getopt::Long;
use v5.10;
use Cwd;

##############
##Command-line
##############
# ./1_mapping_plastid.pl --out_sqlite SQLite/results.db --treated untreated

# get the command line arguments
my ($work_dir,$tmpfolder,$out_sqlite,$treated);

GetOptions(
"tmp:s" =>\$tmpfolder,                  	# Folder where temporary files are stored,                          			optional  argument (default = $TMP or $CWD/tmp env setting)
"work_dir:s" =>\$work_dir,              	# Working directory ,                                               			optional  argument (default = $CWD env setting)
"out_sqlite:s" =>\$out_sqlite,          	# sqlite DB output file,                                             			optional  argument (default = results.db)
"treated:s" =>\$treated                     # Untreated (no treat, CHX,...) or treated (LTM, HARR,...)                      optional  argument (default = 'untreated')
);

###########################################################################
#Check all input variable and/or get default values and set extra variables
###########################################################################

my $CWD             = getcwd;
my $HOME            = $ENV{'HOME'};
my $TMP             = ($ENV{'TMP'}) ? $ENV{'TMP'} : ($tmpfolder) ? $tmpfolder : "$CWD/tmp"; # (1) get the TMP environment variable, (2) get the $tmpfolder variable, (3) get current_working_dir/tmp
print "The following tmpfolder is used                          : $TMP\n";

#Check if tmpfolder exists, if not create it...
if (!-d "$TMP") {
    system("mkdir ". $TMP);
}

# comment on input
if ($work_dir){
    print "Working directory                                        : $work_dir\n";
} else {
    $work_dir = $CWD;
    print "Working directory                                        : $work_dir\n";
}
if (!defined($out_sqlite))     		{$out_sqlite        = $work_dir."/SQLite/results.db";}
print "SQLite database                                          : $out_sqlite\n";

if (!defined($treated)){
    $treated = "untreated";
} elsif ($treated ne "untreated" && $treated ne "treated"){
    print "ERROR: treated parameter should be 'untreated' or 'treated'!";
    die;
}
print "Sample treatment                                              : $treated\n";


#Get all other necessary arguments out of SQLite DB
my $db_sqlite_results  = $out_sqlite;
# Sqlite results
my $dsn_sqlite_results = "DBI:SQLite:dbname=$db_sqlite_results";
my $us_sqlite_results  = "";
my $pw_sqlite_results  = "";

# Get arguments vars
my ($species,$ensemblversion,$IGENOMES_ROOT,$run_name,$bam_untr,$bam_tr,$min_length,$max_length) = get_ARG_vars($dsn_sqlite_results,$us_sqlite_results,$pw_sqlite_results);
# Get executables
my $sqlite_loc = "sqlite3";

#Select right bam file based on treated or untreated
my $bam;
if ($treated eq "untreated"){
    $bam = $bam_untr;
} elsif ($treated eq "treated"){
    $bam = $bam_tr;
}

#Comment on input variables from argument table
print "The following species is used                    : $species\n";
print "The following Ensembl version is used            : $ensemblversion\n";
print "The following igenomes folder is used			: $IGENOMES_ROOT\n";
print "The following bam file is used                   : $bam\n";

#Conversion for species terminology
my $spec = (uc($species) eq "MOUSE") ? "Mus_musculus" : (uc($species) eq "HUMAN") ? "Homo_sapiens" : (uc($species) eq "ARABIDOPSIS") ? "Arabidopsis_thaliana" : (uc($species) eq "FRUITFLY") ? "Drosophila_melanogaster" : "";
my $spec_short = (uc($species) eq "MOUSE") ? "mmu" : (uc($species) eq "HUMAN") ? "hsa" : (uc($species) eq "ARABIDOPSIS") ? "ath" : (uc($species) eq "FRUITFLY") ? "dme" : "";
#Old mouse assembly = NCBIM37, new one is GRCm38. Old human assembly = GRCh37, the new one is GRCh38
my $assembly = (uc($species) eq "MOUSE" && $ensemblversion >= 70 ) ? "GRCm38"
: (uc($species) eq "MOUSE" && $ensemblversion < 70 ) ? "NCBIM37"
: (uc($species) eq "HUMAN" && $ensemblversion >= 76) ? "GRCh38"
: (uc($species) eq "HUMAN" && $ensemblversion < 76) ? "GRCh37"
: (uc($species) eq "ARABIDOPSIS") ? "TAIR10"
: (uc($species) eq "FRUITFLY" && $ensemblversion < 79) ? "BDGP5"
: (uc($species) eq "FRUITFLY" && $ensemblversion >= 79) ? "BDGP6" : "";

# Define genes gtf file path
my $genes_gtf = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Annotation/Genes/genes_".$ensemblversion.".gtf";

############
## PLASTID #
############

my $start = time;

generate_metagene($genes_gtf,$run_name);

index_bam($bam);

calculate_offset($bam,$run_name,$min_length,$max_length);

dump_offsets_in_sqlite($dsn_sqlite_results,$us_sqlite_results,$pw_sqlite_results,$run_name, $treated);

#Resturcture plastid output files
if (!-d "plastid"){
    system("mkdir plastid");
}
if (!-d "plastid/".$treated){
    system("mkdir plastid/".$treated);
}
system("mv ".$run_name."_rois.bed plastid/".$treated."/".$run_name."_".$treated."_rois.bed");
system("mv ".$run_name."_rois.txt plastid/".$treated."/".$run_name."_".$treated."_rois.txt");
system("mv ".$run_name."_metagene_profiles.txt plastid/".$treated."/".$run_name."_".$treated."_metagene_profiles.txt");
system("mv ".$run_name."_p_offsets.txt plastid/".$treated."/".$run_name."_".$treated."_p_offsets.txt");
system("mv ".$run_name."_p_offsets.png plastid/".$treated."/".$run_name."_".$treated."_p_offsets.png");

my $end = time - $start;
printf("runtime plastid: %02d:%02d:%02d\n\n",int($end/3600), int(($end % 3600)/60), int($end % 60));

############
# THE SUBS #
############

## Generate metagene
sub generate_metagene{
    #Catch
    my $genes_gtf = $_[0];
    my $run_name = $_[1];
    
    #Build command
    my $command = "metagene generate -q ".$run_name." --landmark cds_start --annotation_files ".$genes_gtf." 2> /dev/null";
    
    #Execute command
    print "Generate metagene\n".$command."\n\n";
    system($command);
    
    return;
}

## Index bamfile
sub index_bam{
    #Catch
    my $bam = $_[0];
    
    #Build command
    my $command = "samtools index ".$bam;
    
    #Execute command
    print "Index bam\n".$command."\n\n";
    system($command);
    
    return
}


## Run the psite module
sub calculate_offset{
    #catch
    my $bam = $_[0];
    my $run_name = $_[1];
    my $min_l = $_[2];
    my $max_l =$_[3];
    
    #Build command
    my $command = "psite -q ".$run_name."_rois.txt ".$run_name." --min_length ".$min_l." --max_length ".$max_l." --require_upstream --count_files ".$bam." 2> /dev/null";
    
    #Execute command
    print "Calculate psite\n".$command."\n\n";
    system($command);
    
    return;
}

## Dump offsets in SQLite table
sub dump_offsets_in_sqlite{
    # Catch
    my $dsn = $_[0];
    my $user = $_[1];
    my $pw = $_[2];
    my $run_name = $_[3];
    my $treated = $_[4];
    
    #Init
    my %offset_hash;
    
    #Connect to db
    my $dbh_results = dbh($dsn, $user, $pw);
    
    #Create table
    my $create_query = "CREATE TABLE IF NOT EXISTS `p_offsets_".$treated."` (
    `length` int(11) default NULL,
    `offset` int(11) default NULL
    );";
    $dbh_results->do($create_query);
    
    #Read in file
    my $offset_file = $run_name."_p_offsets.txt";
    open(my $FR, $offset_file) or die "Could not open $offset_file: $!";
    while(my $line = <$FR>){
        if($line =~ /^(\d+)\s+(\d+)$/){
            my $length = $1;
            my $offset = $2;
            $offset_hash{$length} = $offset;
        }
    }
    close $FR;
    
    #print Dumper (\%offset_hash);
    #print "\n\n";
    
    my ($last) = sort {$a<=>$b} values %offset_hash; #To determine offsets which were given value 50 in plastid, start for the first with minimal offset value
    for my $key (sort {$a<=>$b} keys %offset_hash) {
        my $length = $key;
        my $offset = $offset_hash{$key};
        if ($offset==50){
            $offset = $last #Take offset of last RPF length to correct for values given value 50
        }
        $last = $offset;
        #print $length."\t".$offset."\n";
        
        #Insert into SQLite table
        my $insert_query = "INSERT INTO p_offsets_".$treated." (length, offset) VALUES (".$length.", ".$offset.");";
        $dbh_results->do($insert_query);
    }
    
    return;
}


##Get paramters from arguments table

sub get_ARG_vars{
    # Catch
    my $dsn = $_[0];
    my $user = $_[1];
    my $pw = $_[2];
    
    my ($query,$sth);
    
    # Connect to db
    my $dbh_results = dbh($dsn, $user, $pw);
    
    # Get ARG variables
    $query = "select value from arguments where variable = \'species\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $species = $sth->fetch()->[0];
    $sth->finish();
    
    $query = "select value from arguments where variable = \'ensembl_version\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $version = $sth->fetch()->[0];
    $sth->finish();
    
    $query = "select value from arguments where variable = \'igenomes_root\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $IGENOMES_ROOT = $sth->fetch()->[0];
    $sth->finish();
    
    $query = "select value from arguments where variable = \'run_name\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $run_name = $sth->fetch()->[0];
    $sth->finish();
    
    $query = "select value from arguments where variable = \'out_bam_untr\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $out_bam_untr = $sth->fetch()->[0];
    $sth->finish();
    
    $query = "select value from arguments where variable = \'out_bam_tr\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $out_bam_tr = $sth->fetch()->[0];
    $sth->finish();
    
    $query = "select value from arguments where variable = \'min_l_plastid\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $min_l_plastid = $sth->fetch()->[0];
    $sth->finish();
    
    $query = "select value from arguments where variable = \'max_l_plastid\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $max_l_plastid = $sth->fetch()->[0];
    $sth->finish();
    
    $dbh_results -> disconnect();
    
    # Return ARG variables
    return($species,$version,$IGENOMES_ROOT,$run_name,$out_bam_untr,$out_bam_tr,$min_l_plastid,$max_l_plastid);
}

### DBH ###
sub dbh {
    
    # Catch
    my $db  = $_[0];
    my $us	= $_[1];
    my $pw	= $_[2];
    
    # Init DB
    my $dbh = DBI->connect($db,$us,$pw,{ RaiseError => 1 },) || die "Cannot connect: " . $DBI::errstr;
    
    return($dbh);
}
