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
# ./Phylocsf_sORFs.pl --sqlite_db SQLite/results.db --cores 8 --tis_ids 1

# get the command line arguments
my ($sqlite_db,$cores,$maf_root,$tis_ids,$out_sqlite,$tmpfolder,$work_dir);

GetOptions(
"sqlite_db=s"=>\$sqlite_db,                 # The sqlite DB holding all RIBO-pipeline results,                     mandatory argument
"cores=i"=>\$cores,                         # Number of cores to use for Bowtie Mapping,                            mandatory argument
"tis_ids=s"  =>\$tis_ids,                   # list of analysis ids                                                  mandatory argument
"out_sqlite=s"=>\$out_sqlite                # The sqlite DB holding all the output                                  mandatory argument
);

#######################
## PhylyCSF Location ##
#######################

my $PhyloCSF_loc = "/data/steven/PhyloCSF/PhyloCSF";

##Galaxy
#my $PhyloCSF_loc = "/data/steven/PhyloCSF/PhyloCSF";

my $CWD             = getcwd;
my $HOME            = $ENV{'HOME'};
my $TMP             = ($ENV{'TMP'}) ? $ENV{'TMP'} : ($tmpfolder) ? $tmpfolder : "$CWD/tmp"; # (1) get the TMP environment variable, (2) get the $tmpfolder variable, (3) get

print "The following tmpfolder is used                          : $TMP\n";
print "The following results db folder is used                  : $sqlite_db\n";

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
if (!defined($out_sqlite))       {$out_sqlite          = $sqlite_db;}

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

#Phylogeny
my $phylogeny = ($species eq "fruitfly") ? "12flies" : ($species eq "mouse") ? "29mammals" : ($species eq "human") ? "29mammals" : "";
#my $chrom_file = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Annotation/Genes/ChromInfo.txt";
my $BIN_chrom_dir = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Sequence/Chromosomes_BIN";

############################################
##                                        ##
## Creating MAF based PhyloCSF input data ##
##                                        ##
############################################

# Start time
my $start = time;

# Get the analysis_id that corresponds to the TIS-calling input parameters
my $idsref = get_analysis_ids($dsn_results,$us_results,$pw_results,$tis_ids);

## Get chromosome sizes

print "Getting chromosome sizes  ...\n";

my $chr_sizes = get_chr_sizes();

print "\nGet chromosomes... \n\n";
## Get chromosomes based on seq_region_id ##
my $chrs = get_chrs($dsn_ENS,$us_ENS,$pw_ENS,get_chr_sizes(),$assembly,$species);

my $analysis_id;
#Loop over all selected analysis_ids
print "Analysing sORF conservation with PhyloCSF...\n\n";

foreach $analysis_id (@$idsref) {
    
    print "Processing analysis_id $analysis_id ...\n\n";
    PhyloCSF_analyse_sorfs($chrs,$TMP,$analysis_id,$phylogeny);
    
    print "Storing all PhyloCSF scores in DB ... \n\n";
    parse_PhyloCSF_out($chrs,$TMP,$analysis_id,$work_dir)
}

#Move to galaxy history
system ("mv ".$sqlite_db." ".$out_sqlite);

# End time
print "   DONE! \n";
my $end = time - $start;
printf("runtime assembly: %02d:%02d:%02d\n\n",int($end/3600), int(($end % 3600)/60), int($end % 60));


############
# THE SUBS #
############

### Parse PhyloCSF output ###

sub parse_PhyloCSF_out{

    #Catch
    my $chrs        =   $_[0];
    my $TMP         =   $_[1];
    my $analysis_id =   $_[2];
    my $work_dir    =   $_[3];
    
    #Init
    my $dbh = dbh($dsn_results,$us_results,$pw_results);
    my $sorf_id;
    my $score;
    
    ### Get sORF transcript table from DB
    my $sorfs = get_sorfs($dbh,$analysis_id);
    
    ### Parse PhyloCSF scores from files
    foreach my $chr (keys %{$chrs}){
    open(FH, "<", $TMP."/".$analysis_id."_".$chr."_PhyloCSF_out.txt") or die "cannot open < input.txt: $!";
        foreach my $line (<FH>){
            if ($line =~ /^\S/){
                $line =~ /\/(\w*)\.fa\s*\S*\s*((-|)\d*\.\d*)/;
                $sorf_id   =   $1;
                $score     =   $2;
                
                $sorfs->{$sorf_id}{'score'} = $score;
            }
        }
        close FH;
        system ("rm ".$TMP."/".$analysis_id."_".$chr."_PhyloCSF_out.txt");
    }
    
    ### Save in CSV
    store_in_csv($sorfs,$TMP);
    
    ### Save in DB
    store_in_db($dsn_results,$us_results,$pw_results,$analysis_id,$work_dir,$TMP);
}

### Store in DB ###

sub store_in_db{
    
    # Catch
    my $dsn         =   $_[0];
    my $us          =   $_[1];
    my $pw          =   $_[2];
    my $id          =   $_[3];
    my $work_dir    =   $_[4];
    my $TMP         =   $_[5];
    
    # Init
    my $dbh     =   dbh($dsn,$us,$pw);
    my $table   =   "TIS_sORFs_".$id."_transcripts_PhyloCSF";
    
    #drop old table
    my $query = "drop table if exists `".$table."` ";
    $dbh->do($query);
    
    # Create table
    $query = "CREATE TABLE IF NOT EXISTS `".$table."` (
    `sorf_id` int(10),
    `chr` char(50) NOT NULL default '',
    `strand` int(2) NOT NULL default '',
    `sorf_begin` int(10) NOT NULL default '',
    `sorf_end` int(10) NOT NULL default '',
    `PhyloCSF` decimal(11,4) NOT NULL default '')"  ;
    $dbh->do($query);
    
    system("sqlite3 -separator , ".$sqlite_db." \".import ".$TMP."/sorf_transcripts.csv ".$table."\"")== 0 or die "system failed: $?";
    
    #Disconnect dbh
    $dbh->disconnect();
    system("rm ".$TMP."/sorf_transcripts.csv ");

}

### Store sORFs in csv ###

sub store_in_csv{
    
    #Catch
    my $sorfs   =   $_[0];
    my $TMP     =   $_[1];
    
    #Init
    open TMP_sorfs, "+>>".$TMP."/sorf_transcripts.csv" or die $!;
    
    #Save each sorf in csv
    foreach my $sorf (keys %{$sorfs}){
        
        print TMP_sorfs $sorfs->{$sorf}{'sorf_id'}.",".$sorfs->{$sorf}{'chr'}.",".$sorfs->{$sorf}{'strand'}.",".$sorfs->{$sorf}{'sorf_begin'}.",".$sorfs->{$sorf}{'sorf_end'}.",".$sorfs->{$sorf}{'score'}."\n";
    }
    close TMP_sorfs
}

### PhyloCSF analyse sORFs ###

sub PhyloCSF_analyse_sorfs{

    #Catch
    my $chrs        =   $_[0];
    my $tmp         =   $_[1];
    my $analysis_id =   $_[2];
    my $phylogeny   =   $_[3];
    
    # Init multi core
    my $pm = new Parallel::ForkManager($cores);
    
    ## Loop over all chromosomes
    foreach my $chr (sort keys %{$chrs}){
        
        ### Start parallel process
        $pm->start and next;
        
        ### DBH per process
        my $dbh = dbh($dsn_results,$us_results,$pw_results);
        
        ### Create files
        open TMP_phylo, "+>>".$TMP."/".$analysis_id."_".$chr."_PhyloCSF_in.txt" or die $!;
        open TMP_phylo_out, "+>>".$TMP."/".$analysis_id."_".$chr."_PhyloCSF_out.txt" or die $!;
        system("mkdir ".$TMP."/".$chr." ");
        
        ### Get sORF alignments per chromosome from transcript_maf table
        my $sorf_align = get_sorf_align_per_chr($dbh,$analysis_id,$chr);
        
        ### Run over sORFs, create alignments
        foreach my $sorf_align_id (keys %{$sorf_align}){
            
            ## Create and store in files
            store_sorf_align($sorf_align,$sorf_align_id,$chr,$TMP);
        }
        
        #Close db and writable TMP
        $dbh->disconnect();
        close TMP_phylo;
        
        ### Run PhyloCSF on TMP_phylo files
        my $cmd = $PhyloCSF_loc." ".$phylogeny." --files ".$TMP."/".$analysis_id."_".$chr."_PhyloCSF_in.txt --removeRefGaps";
        #my $cmd = "PhyloCSF ".$phylogeny." --files ".$TMP."/".$analysis_id."_".$chr."_PhyloCSF_in.txt --removeRefGaps";
        my $output = qx/$cmd/;
        print TMP_phylo_out $output;
        
        #Finish chromosome analysis
        close TMP_phylo_out;
        system ("rm ".$TMP."/".$analysis_id."_".$chr."_PhyloCSF_in.txt");
        system ("rm -rf ".$TMP."/".$chr);
        print "     * Finished conservation analysis on chromosome ".$chr."\n";
        
        ### Finish child
        $pm->finish;
    }
    
    #Wait all Children
    $pm->wait_all_children();
}

### Store sORF aligns in files ###

sub store_sorf_align{

    #Catch
    my $sorf_align      =   $_[0];
    my $sorf_align_id   =   $_[1];
    my $chr             =   $_[2];
    my $TMP             =   $_[3];
    
    #Init
    open TMP_phylo_sorf, "+>>".$TMP."/".$chr."/".$sorf_align_id.".fa" or die $!;
    my ($PhyloCSF_species,$PhyloCSF_spec);
    if ($species eq "fruitfly"){
        $PhyloCSF_species =  ['dmel','dsim','dsec','dyak','dere','dana','dpse','dper','dwil','dvir','dmoj','dgri'];
        $PhyloCSF_spec = "dmel";
    }elsif ($species eq "human"){
        $PhyloCSF_species = ['Human','Chimp','Rhesus','Bushbaby','TreeShrew','Mouse','Rat','Guinea_Pig','Squirrel','Rabbit','Pika','Alpaca','Dolphin','Cow','Horse','Cat','Dog','Microbat','Megabat','Hedgehog','Shrew','Elephant','Tenrec','Armadillo'];
        $PhyloCSF_spec = "Human";
    }elsif ($species eq "mouse"){
        $PhyloCSF_species =  ['Mouse','Guinea_Pig','Kangaroo_rat','Pika','Rabbit','Rat','Squirrel','TreeShrew','Human','Mouse_lemur','Bushbaby','Chimp','Rhesus','Tarsier','Cow','Dog','Sloth','Armadillo','Tenrec','Horse','Hedgehog','Cat','Elephant','Microbat','Rock_hyrax','Megabat','Shrew','Dolphin','Alpaca'];
        $PhyloCSF_spec = "Mouse";
    }
    
    #Save TMP_phylo_sorf location in TMP_phylo
    print TMP_phylo $TMP."/".$chr."/".$sorf_align_id.".fa\n";
    
    #Save alignment in fasta file
    
    #tmp
    my $length = length($sorf_align->{$sorf_align_id}{$PhyloCSF_spec});
    
    foreach (@{$PhyloCSF_species}){
        
        #tmp
        if (length($sorf_align->{$sorf_align_id}{$_}) == $length){
        
        if ($sorf_align->{$sorf_align_id}{$_} ne ""){
            if ($_ eq $PhyloCSF_spec){
                print TMP_phylo_sorf ">".$_."|".$sorf_align_id."\n";
                print TMP_phylo_sorf $sorf_align->{$sorf_align_id}{$_}."\n";
            }else{
                print TMP_phylo_sorf ">".$_."\n";
                print TMP_phylo_sorf $sorf_align->{$sorf_align_id}{$_}."\n";
            }
        }
         
        #tmp
        }
        
    }
    close TMP_phylo_sorf;
}

### get sORFs ###

sub get_sorfs {
    
    # Catch
    my $dbh         =   $_[0];
    my $id          =   $_[1];
    
    # Init
    my $sorfs = {};
    
    # Get genes
    my $query = "SELECT * from TIS_sORFs_".$id."_transcripts";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    $sorfs = $sth->fetchall_hashref('sorf_id');
    
    
    ## Check if able to delete
    
    foreach my $sorf (keys %{$sorfs}){
        $sorfs->{$sorf}{'score'} = "";
    }
    
    # Return
    return($sorfs);
}

### get sORFs per chromosome ###

sub get_sorf_align_per_chr {
    
    # Catch
    my $dbh         =   $_[0];
    my $id          =   $_[1];
    my $chr         =   $_[2];
    
    # Init
    my $sorf_align = {};
    
    # Get genes
    my $query = "SELECT * from TIS_sORFs_".$id."_transcripts_maf WHERE chr = '".$chr."'";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    $sorf_align = $sth->fetchall_hashref('sorf_id');
    
    # Return
    return($sorf_align);
    
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
    
    # Get correct coord_system_id
    my $query = "SELECT coord_system_id FROM coord_system where name = 'chromosome' and version = '".$assembly."'";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    @coord_system = $sth->fetchrow_array;
    $coord_system_id = $coord_system[0];
    $sth->finish();
    
    # Get chrs with seq_region_id
    my $chr;
    #foreach (@chr){
    foreach my $key (keys(%{$chr_sizes})) {
        if ($species eq "fruitfly"){
            if($key eq "M"){
                $chr = "dmel_mitochondrion_genome";
            }else{
                $chr = $key;
            }
        }else {
            $chr = $key;
        }
        
        my $query = "SELECT seq_region_id FROM seq_region where coord_system_id = ".$coord_system_id."  and name = '".$chr."' ";
        my $sth = $dbh->prepare($query);
        $sth->execute();
        @ids = $sth->fetchrow_array;
        $seq_region_id = $ids[0];
        $chrs->{$chr}{'seq_region_id'} = $seq_region_id;
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
    my $us  = $_[1];
    my $pw  = $_[2];
    
    # Init DB
    my $dbh = DBI->connect($db,$us,$pw,{ RaiseError => 1 },) || die "Cannot connect: " . $DBI::errstr;
    
    return($dbh);
}

### GET CHR SIZES ###

sub get_chr_sizes {
    
    # Catch
    my %chr_sizes;
    my $filename = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Annotation/Genes/ChromInfo.txt";
    if (-e $filename) {
        # Work
        open (Q,"<".$filename) || die "Cannot open chr sizes input\n";
        while (<Q>){
            my @a = split(/\s+/,$_);
            $chr_sizes{$a[0]} = $a[1];
        }
    }
    else
    {
        $filename = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Sequence/WholeGenomeFasta/GenomeSize.xml";
        open (Q,"<".$filename) || die "Cannot open chr sizes input\n";
        while (<Q>){
            if ($_ =~ /contigName=\"(.*)\".*totalBases=\"(\d+)\"/) {
                $chr_sizes{$1} = $2;
            }
        }
    }
    return(\%chr_sizes);
    
}