#!/usr/bin/perl -w
use strict;
use warnings;
use DBI;
use DBD::SQLite;
use Data::Dumper;
use Getopt::Long;
use v5.10;
use Parallel::ForkManager;
use Cwd;

##############
##Command-line
# ./Translation_assembly.pl --sqliteRES SQLite/results_sORFs_LTM.db --cores 23 --tis_ids 1

# get the command line arguments
my ($sqlite_db,$cores,$maf_root,$tis_ids,$out_sqlite,$tmpfolder,$work_dir);

GetOptions(
"sqliteRES=s"=>\$sqlite_db,                 # The sqlite DB holding all RIBO-pipeline results,                     mandatory argument
"cores=i"=>\$cores,                         # Number of cores to use for Bowtie Mapping,                            mandatory argument
"tis_ids=s"  =>\$tis_ids,                   # list of analysis ids                                                  mandatory argument
"out_sqlite=s"=>\$out_sqlite                # The sqlite DB holding all the output                                  mandatory argument
);

my $CWD             = getcwd;
my $HOME            = $ENV{'HOME'};
my $TMP             = ($ENV{'TMP'}) ? $ENV{'TMP'} : ($tmpfolder) ? $tmpfolder : "$CWD/tmp"; # (1) get the TMP environment variable, (2) get the $tmpfolder variable, (3) get

print "The following tmpfolder is used                          : $TMP\n";
print "The following results db folder is used                  : $CWD/$sqlite_db\n";

#Check if tmpfolder exists, if not create it...
if (!-d "$TMP") {
    system ("mkdir ". $TMP);
}

#comment on these
if ($out_sqlite){
    print "Output SqliteDB used is                                  : $CWD/$out_sqlite\n";
} else {
    $out_sqlite = $sqlite_db;
    print "Output SqliteDB used is                                  : $out_sqlite\n";
}
if ($work_dir){
    print "Working directory                                        : $work_dir\n";
} else {
    #Choose default value
    $work_dir = $CWD;
    print "Working directory                                        : $CWD\n";
}
if ($cores){
    print "Number of cores to use for assembling                    : $cores\n";
} else {
    die "\nDon't forget to pass number of cores to use for mapping using the --cores or -c argument!\n\n";
}

# Create output files for command line script
if (!defined($out_sqlite))       {$out_sqlite          = $work_dir."/".$sqlite_db;}

# DB settings
# Sqlite Riboseq
my $db_results  = $sqlite_db;
my $dsn_results = "DBI:SQLite:dbname=$db_results";
my $us_results  = "";
my $pw_results  = "";

#Get arguments from arguments table
my ($ensemblversion,$species,$ens_db,$IGENOMES_ROOT) = get_arguments($dsn_results,$us_results,$pw_results);

print "The igenomes_root folder used is                         : $IGENOMES_ROOT\n";
print "Number of cores to use for Mapping                       : $cores\n";
print "The following Ensembl db folder is used                  : $ens_db\n";

#Conversion for species terminology
my $spec = ($species eq "mouse") ? "Mus_musculus" : ($species eq "human") ? "Homo_sapiens" : ($species eq "arabidopsis") ? "Arabidopsis_thaliana" : ($species eq "fruitfly") ? "Drosophila_melanogaster" : "";
my $spec_short = ($species eq "mouse") ? "mmu" : ($species eq "human") ? "hsa" : ($species eq "arabidopsis") ? "ath" : ($species eq "fruitfly") ? "dme" : "";

# Sqlite Ensembl
my $db_ENS  = $ens_db;
my $dsn_ENS = "DBI:SQLite:dbname=$db_ENS";
my $us_ENS  = "";
my $pw_ENS  = "";

#Old mouse assembly = NCBIM37, new one is GRCm38
my $assembly = ($species eq "mouse" && $ensemblversion >= 70 ) ? "GRCm38"
: ($species eq "mouse" && $ensemblversion < 70 ) ? "NCBIM37"
: ($species eq "human") ? "GRCh37"
: ($species eq "arabidopsis") ? "TAIR10"
: ($species eq "fruitfly") ? "BDGP5" : "";

my $chrom_file = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Annotation/Genes/ChromInfo.txt";

############################################
##                                        ##
## Creating transcript based fasta DB     ##
##                                        ##
############################################

# Start time
my $start = time;

# Get the analysis_id that corresponds to the TIS-calling input parameters
my $idsref = get_analysis_ids($dsn_results,$us_results,$pw_results,$tis_ids);

print "\nGet chromosomes... \n\n";
## Get chromosomes based on seq_region_id ##
my $chrs = get_chrs($dsn_ENS,$us_ENS,$pw_ENS,$chrom_file,$assembly,$species);

my $analysis_id;
#Loop over all selected analysis_ids
print "Creating fasta output for sORF transcripts...\n\n";

foreach $analysis_id (@$idsref) {
    
    print "Processing analysis_id $analysis_id ...\n\n";
    assemble_sORF_transcripts($chrs,$TMP,$analysis_id);
    
}

# End time
print "   DONE! \n";
my $end = time - $start;
printf("runtime assembly: %02d:%02d:%02d\n\n",int($end/3600), int(($end % 3600)/60), int($end % 60));


############
# THE SUBS #
############

### Assemble sORF Transcripts ###

sub assemble_sORF_transcripts{

    #Catch
    my $chrs        =   $_[0];
    my $TMP         =   $_[1];
    my $analysis_id =   $_[2];
    
    # Init multi core
    my $pm = new Parallel::ForkManager($cores);
    
    ## Loop over all chromosomes
    foreach my $chr (sort keys %{$chrs}){
        
        ### Start parallel process
        $pm->start and next;
        
        ### DBH per process
        my $dbh = dbh($dsn_results,$us_results,$pw_results);
        
        ### Output fasta and db_csv file per process
        open TMP_db, "+>>".$TMP."/".$analysis_id."_".$chr.".fasta" or die $!;
        
        ### Get sORF transcript table from DB
        my $sorfs = get_sorfs($dbh,$analysis_id,$chr);
    
        ### Output assembled sORF transcripts
        foreach my $sorf (keys %{$sorfs}){
        
            #Init
            my $sorf_id     =   $sorf;
            #chromosoom
            my $sorf_begin  =   $sorfs->{$sorf}{'sorf_begin'};
            my $sorf_end    =   $sorfs->{$sorf}{'sorf_end'};
            my $strand      =   $sorfs->{$sorf}{'strand'};
            my $annotation  =   $sorfs->{$sorf}{'annotation'};
            my $start_codon =   $sorfs->{$sorf}{'start_codon'};
            chop($sorfs->{$sorf}{'aa_seq'});
            my $aa_seq      =   $sorfs->{$sorf}{'aa_seq'};
            my $mass        =   $sorfs->{$sorf}{'mass'};
            my $cons        =   $sorfs->{$sorf}{'PhyloCSF'};
            my $biotype     =   $sorfs->{$sorf}{'biotype'};
            
            my $identifier  =   ">generic|".$sorf_id."_".$chr."_".$strand."_".$sorf_begin."_".$sorf_end."_".$annotation."|".$start_codon."_".$biotype."_".$mass."_".$cons." ";
            #my $identifier  =   ">generic|".$sorf_id."_".$chr."_".$strand."_".$sorf_begin."_".$sorf_end."_".$annotation."|".$start_codon."_".$cons." ";
            #my $identifier  =   ">generic|".$sorf_id."_".$chr."_".$strand."_".$sorf_begin."_".$sorf_end."_".$annotation."|".$start_codon."_".$mass." ";
            
            print TMP_db $identifier."\n";
            print TMP_db $aa_seq."\n";
        
        }
        
        ### Close tmpt output
        close(TMP_db);
        
        ### Finish childs
        print "     * Finished translating chromosome ".$chr."\n";
        $dbh->disconnect();
        $pm->finish;
    }
    #Waiting for all childs to finish
    $pm->wait_all_children();
    
    #Check if outputfolder exists, if not create it...
    my $outputfolder = "$CWD/output";
    if (!-d "$outputfolder") {
        system ("mkdir ". $outputfolder);
    }
    
    #Concatenate TMP sORF fasta DBs
    #open DB, "+>> TIS_sORFs_".$analysis_id."_transcripts.fasta" or die $!;

    ## Loop over all chromosomes
    foreach my $chr (sort keys %{$chrs}){
        
        system("cat ".$TMP."/".$analysis_id."_".$chr.".fasta >> ".$outputfolder."/TIS_sORFs_".$analysis_id."_transcripts.fasta");
    	system("rm -rf ".$TMP."/".$analysis_id."_".$chr.".fasta");
    }
}

### get sORFs ###

sub get_sorfs {
    
    # Catch
    my $dbh         =   $_[0];
    my $id          =   $_[1];
    my $chr         =   $_[2];
    
    # Init
    my $sorfs = {};
    
    # Get genes
    my $query = "SELECT * from TIS_sORFs_".$id."_transcripts where chr = '".$chr."' and exon_overlap < '0.75'";
	my $sth = $dbh->prepare($query);
	$sth->execute();
	$sorfs = $sth->fetchall_hashref('sorf_id');
    
	# Return
	return($sorfs);
}

### Get the analysis ids that need to be processed

sub get_analysis_ids {
    
    # Catch
    my $dsn     =   $_[0];
    my $us      =   $_[1];
    my $pw      =   $_[2];
    my $ids_in  =   $_[3]; #Either comma separated list of identifiers or "all"
    
    #Init
    my $dbh = dbh($dsn,$us,$pw);
    my $idsref;
    if ($ids_in eq "all") {
        my $query_analysis_id = "select ID from TIS_sORFs_overview";
        my @ids = @{$dbh->selectcol_arrayref($query_analysis_id)};
        $idsref = \@ids;
    }
    else {
        my @ids = split(/,/,$ids_in);
        $idsref = \@ids;
        
    }
    
    return $idsref;
}

### GET CHRs ###

sub get_chrs {
    
    # Catch
    my $db          =   $_[0];
    my $us          =   $_[1];
    my $pw          =   $_[2];
    my $chr_file    =   $_[3];
    my $assembly    =   $_[4];
    
    # Init
    my $chrs    =   {};
    my $dbh     =   dbh($db,$us,$pw);
    my ($line,@chr,$coord_system_id,$seq_region_id,@ids,@coord_system);
    
    # Get chrs from Chr_File
    open (Q,"<".$chr_file) || die "Cannot open chr sizes input\n";
    while ($line = <Q>){
        $line =~ /^(\S*)/;
        push (@chr,$1);
    }
    
    # Get correct coord_system_id
    my $query = "SELECT coord_system_id FROM coord_system where name = 'chromosome' and version = '".$assembly."'";
	my $sth = $dbh->prepare($query);
	$sth->execute();
    @coord_system = $sth->fetchrow_array;
    $coord_system_id = $coord_system[0];
    $sth->finish();
   	
    # Get chrs with seq_region_id
    my $chr;
    foreach (@chr){
        if ($species eq "fruitfly"){
            if($_ eq "M"){
                $chr = "dmel_mitochondrion_genome";
            }else{
                $chr = $_;
            }
        }else {
            $chr = $_;
        }
        
        my $query = "SELECT seq_region_id FROM seq_region where coord_system_id = ".$coord_system_id."  and name = '".$chr."' ";
        my $sth = $dbh->prepare($query);
        $sth->execute();
        @ids = $sth->fetchrow_array;
        $seq_region_id = $ids[0];
        $chrs->{$_}{'seq_region_id'} = $seq_region_id;
        $sth->finish();
    }
    
    #Disconnect DBH
    $dbh->disconnect();
	
	# Return
	return($chrs);
    
}

### Get arguments

sub get_arguments{
    
    # Catch
    my $dsn                 =   $_[0];
    my $us                  =   $_[1];
    my $pw                  =   $_[2];
    
    # Init
    my $dbh = dbh($dsn,$us,$pw);
    
    # Get input variables
    my $query = "select value from `arguments` where variable = \'ensembl_version\'";
    my $sth = $dbh->prepare($query);
	$sth->execute();
	my $ensemblversion = $sth->fetch()->[0];
    
    $query = "select value from `arguments` where variable = \'species\'";
    $sth = $dbh->prepare($query);
	$sth->execute();
	my $species = $sth->fetch()->[0];
    
    $query = "select value from `arguments` where variable = \'ens_db\'";
    $sth = $dbh->prepare($query);
	$sth->execute();
	my $ens_db = $sth->fetch()->[0];
    
    $query = "select value from `arguments` where variable = \'igenomes_root\'";
    $sth = $dbh->prepare($query);
	$sth->execute();
	my $igenomes_root = $sth->fetch()->[0];
    
    # Return input variables
    return($ensemblversion,$species,$ens_db,$igenomes_root);
    
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
