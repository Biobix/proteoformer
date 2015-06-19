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
# ./remove_lincRNA.pl --sqliteRES SQLite/results_sORFs_LTM.db --cores 1 --tis_ids 2

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
print "Remove intronic LincRNAs \n\n";

foreach $analysis_id (@$idsref) {
    
    print "Processing analysis_id $analysis_id ...\n\n";
    remove_intronic_lncrna($chrs,$TMP,$analysis_id);
    
}

# End time
print "   DONE! \n";
my $end = time - $start;
printf("runtime assembly: %02d:%02d:%02d\n\n",int($end/3600), int(($end % 3600)/60), int($end % 60));


############
# THE SUBS #
############

### Assemble sORF Transcripts ###

sub remove_intronic_lncrna{

    #Catch
    my $chrs        =   $_[0];
    my $TMP         =   $_[1];
    my $analysis_id =   $_[2];
    
    #init
    my $gene_id_transcript_id;
        
    ### DBH per process
    my $dbh = dbh($dsn_results,$us_results,$pw_results);
        
    ### Get sORF transcript table from DB
    my $sorfs = get_sorfs($dbh,$analysis_id);
    
    ### remove intronic lncRNA sORFS
    foreach my $sorf (keys %{$sorfs}){
        
        #Init
        my $test        =   'No';
        my $sorf_id     =   $sorf;
        my $chr         =   $sorfs->{$sorf}{'chr'};
        my $sorf_begin  =   $sorfs->{$sorf}{'sorf_begin'};
        my $sorf_end    =   $sorfs->{$sorf}{'sorf_end'};
        my $strand      =   $sorfs->{$sorf}{'strand'};
        my $gene_id     =   $sorfs->{$sorf}{'id'};
        
   
        #Get transcript per gene_id
        my $trs = get_gene_id_transcripts($dsn_ENS,$us_ENS,$pw_ENS,$gene_id);
        foreach $gene_id_transcript_id (keys %{$trs}){
            
            #Get exons per transcript id
            my $exons = get_exon_data($gene_id_transcript_id,$dsn_ENS,$us_ENS,$pw_ENS);
            
            foreach my $exon_rank (keys %{$exons}){
                
                if ($sorf_begin >= $exons->{$exon_rank}{'seq_region_start'} and $sorf_end <= $exons->{$exon_rank}{'seq_region_end'}){
                    $test = 'Yes';
                    next;
                }
            }
        
        }
        
        if ($test eq 'No'){
        
            #remove from db
            
            my $query = "delete from TIS_sORFs_".$analysis_id."_transcripts where sorf_id = '".$sorf."'";
            my $sth = $dbh->prepare($query);
            $sth->execute();
            $query = "delete from TIS_sORFs_".$analysis_id."_transcripts_PhyloCSF where sorf_id = '".$sorf."'";
            $sth = $dbh->prepare($query);
            $sth->execute();
            $query = "delete from TIS_sORFs_".$analysis_id."_transcripts_FLOSS where sorfID = '".$sorf."'";
            $sth = $dbh->prepare($query);
            $sth->execute();
            
        
        }

    }

}

sub get_exon_data {
    
    # Catch
    my $tr_id           =   $_[0];
    my $us_ENS          =   $_[1];
    my $pw_ENS          =   $_[2];
    my $gene_id         =   $_[3];
    
    # Init
    my $exons = {};
    
    #Init
    my $dbh_exon = dbh($dsn_ENS,$us_ENS,$pw_ENS);
    
    # Get exons
    my $query = "select tr.transcript_id tr_id, tr.stable_id tr_stable_id,etr.rank,e.seq_region_start,e.seq_region_end,e.seq_region_strand,e.phase,e.end_phase,e.stable_id e_stable_id ".
    "from transcript tr ".
    "   inner join exon_transcript etr on tr.transcript_id = etr.transcript_id ".
    "   inner join exon e on e.exon_id = etr.exon_id ".
    " where tr.transcript_id = '". $tr_id ."'";
    my $sth = $dbh_exon->prepare($query);
    $sth->execute();
    $exons = $sth->fetchall_hashref('rank');
    
    #Return
    return($exons);
    
}

### get gene_id transcripts

sub get_gene_id_transcripts{
    
    #Catch
    my $dsn_ENS     =   $_[0];
    my $us_ENS      =   $_[1];
    my $pw_ENS      =   $_[2];
    my $gene_id     =   $_[3];
    
    #Init
    my $dbh_ENS = dbh($dsn_ENS,$us_ENS,$pw_ENS);
    my $trs = {};
    
    #Get transcripts
    my $query = "SELECT transcript_id,gene_id,seq_region_id,seq_region_strand,seq_region_start,seq_region_end,biotype,stable_id,canonical_translation_id from transcript where gene_id = '".$gene_id."' ";
    my $sth = $dbh_ENS->prepare($query);
    $sth->execute();
    $trs = $sth->fetchall_hashref('transcript_id');
    
    #Return
    return($trs);
    
}

### get sORFs ###

sub get_sorfs {
    
    # Catch
    my $dbh         =   $_[0];
    my $id          =   $_[1];
    
    # Init
    my $sorfs = {};
    
    # Get genes
    my $query = "select * from TIS_sORFs_".$id."_transcripts where biotype = 'lincRNA'";
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
