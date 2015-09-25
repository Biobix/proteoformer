#!/usr/bin/perl -w
use strict;
use warnings;
use DBI;
use DBD::SQLite;
use Data::Dumper;
use Cwd;
use Getopt::Long;
use Storable 'dclone';
use Parallel::ForkManager;
#use lib "/home/steven/perl5/lib/perl5";
use Statistics::R;

#################
# Perl commands #
#################

## Command line ##

# ./FlossProteoformer.pl --sqlite SQLite/results.db  --tis_ids 1
## Galaxy ##


# get the command line arguments
my ($work_dir,$out_sqlite,$tmpfolder,$sqlite_db, $tis_ids);

GetOptions(
"dir=s"=>\$work_dir,                            # Path to the working directory                                                                 optional argument
"tmp:s" =>\$tmpfolder,                          # Folder where temporary files are stored,                                                      optional  argument (default = $TMP env setting)
"sqlite_db:s" =>\$sqlite_db,                    # SQLite results db with mapping and tr_translation                                             mandatory argument
"tis_ids=s"  =>\$tis_ids,                       # list of analysis ids                                                                          mandatory argument
"out_sqlite:s" =>\$out_sqlite                   # Galaxy specific history file location                                                         Galaxy specific
);

my $CWD             = getcwd ||die "CWD failed \n\n";
my $HOME            = $ENV{'HOME'};
my $TMP             = ($ENV{'TMP'}) ? $ENV{'TMP'} : ($tmpfolder) ? $tmpfolder : "$CWD/tmp" ; # First select the TMP environment variable, second select the $tmpfolder variable, lastly select current_working_dir/tmp
print "The following tmpfolder is used                          : $TMP\n";
print "The following results db folder is used                  : $CWD/$sqlite_db\n";

#Check if tmpfolder exists, if not create it...
if (!-d "$TMP") {
    system ("mkdir ".$TMP);
}

# comment on these
if ($work_dir){
    print "The following working directory is used                  : $work_dir\n";
} else {
    $work_dir = $CWD;
    print "The following working directory is used                  : $CWD\n";
}
if ($tis_ids){
    print "Number of TIS ids                                        : $tis_ids\n";
} else {
    die "\nDon't forget to pass number of TIS ids to use for FLOSS calculation using the --tis_ids argument!\n\n";
}

# Create output files for command line script
if (!defined($out_sqlite))       {$out_sqlite          = $work_dir."/".$sqlite_db;}

# R library location (zoo package)
my $RlibLocation = "/data/steven/sorf_pipeline_mmu";

# DB settings
# Sqlite Riboseq
my $db_results  = $sqlite_db;
my $dsn_results = "DBI:SQLite:dbname=$db_results";
my $us_results  = "";
my $pw_results  = "";


#Get arguments from arguments table
my ($ensemblversion,$ens_db,$species,$IGENOMES_ROOT,$cores) = get_arguments($dsn_results,$us_results,$pw_results);

print "The igenomes_root folder used is                         : $IGENOMES_ROOT\n";
print "Number of cores to use for calculating FLOSS             : $cores\n";
print "The following Ensembl db folder is used                  : $ens_db\n";
print "The followin Ensembl version is used                     : $ensemblversion\n";
print "Floss is calculated for the following species            : $species\n";

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

print "--------------------------------------------------------------------------------------\n\n\n\n";

######################
## Floss Calculation #
##                   #
######################



# Start time
my $start = time;

#Make connection with results DB
my $dbh = dbh($dsn_results, $us_results, $pw_results);

#Calculating FLOSS only for CHX sample
my $numberCHX = 1;
my $fastqname = "fastq".$numberCHX;

## Get chromosomes based on seq_region_id ##
my $chrs = get_chrs($dsn_ENS,$us_ENS,$pw_ENS,$chrom_file,$assembly,$species);

#Create reference histogram -> DB (and Length distribution of coding transcripts for cutoff)
my $refFracs = refHist($dsn_results,$us_results,$pw_results,$dsn_ENS,$us_ENS,$pw_ENS, $cores, $chrs, $fastqname);

#Create cutoff -> DB
my $smoother = cutoff($dsn_results,$us_results,$pw_results, $chrs, $fastqname,$db_results);

# Get the analysis_id that corresponds to the TIS-calling input parameters
my $idsref = get_analysis_ids($dsn_results,$us_results,$pw_results,$tis_ids);

foreach my $analysis_id (@$idsref) {
    
    # Init multi core
    my $pm = new Parallel::ForkManager($cores);
    print "Initiating multi core for calculating length distribution and FLOSS scores. \n";
    print "   Using ".$cores." core(s)\n   ---------------\n";
    
    foreach my $chr (sort keys %{$chrs}){

        ### Start parallel process
        $pm->start and next;
        
        # seq_region_id per chromosome
        my $seq_region_id = $chrs->{$chr}{'seq_region_id'};
    
        #Get all data from database and store in hash
        print "Retrieving Length Distribution for chromosome $chr\n";
        my $LD = makeLD($dsn_results, $us_results, $pw_results,$chr, $seq_region_id, $fastqname, $analysis_id, $dsn_ENS, $us_ENS, $pw_ENS);
        
        #Calculate FLOSS for chromosome and classify against cutoffs
        print "Do calculation of FLOSS scores and classification for chromosome $chr\n";
        classify($chr, $fastqname, $smoother, $LD, $refFracs);
    
        ### Finish
        print "* Finished chromosome ".$chr."\n";
        $pm->finish;
    }
    
    # Finish all subprocesses
    $pm->wait_all_children;
    
    #init
    my $classTable = {};
    
    #Convert tmpfiles of all choromosomes to one hash
    foreach my $chr(sort keys %{$chrs}){
        open(classRead,"<".$TMP."/classification/classification_".$fastqname."_".$chr.".txt") || die "ERROR reading classification tmp file \n";
        while(my $read = <classRead>){
            chomp($read);
            my @classChr = split("\t",$read);
            $classTable->{$classChr[0]}->{"nreads"} = $classChr[1];
            $classTable->{$classChr[0]}->{"FLOSS"} = $classChr[2];
            if($classTable->{$classChr[0]}->{"FLOSS"} eq "undef"){
                $classTable->{$classChr[0]}->{"FLOSS"}=undef;
            }
            $classTable->{$classChr[0]}->{"Class"} = $classChr[3];
        }
        close(classRead);
    }

    #Create new table to store FLOSS scores and classification
    my $create_table = "CREATE TABLE `TIS_".$analysis_id."_transcripts_FLOSS` (`Table ID` INTEGER PRIMARY KEY, `TIS ID` VARCHAR(25), `nreads` INTEGER, `FLOSS` DECIMAL(9,7),`classification` VARCHAR(15));";
    $dbh->do($create_table);
    
    #Store data in csv file
    print "export FLOSS data to csv file \n";
    open FLOSSprint, ">>".$TMP."/FLOSSclass.csv" or die "ERROR exporting FLOSS data";
    my $count=1;
    foreach my $ID (keys %{$classTable}){
        if (defined $classTable->{$ID}->{"FLOSS"}){
            print FLOSSprint $count.",".$ID.",".$classTable->{$ID}->{"nreads"}.",".$classTable->{$ID}->{"FLOSS"}.",".$classTable->{$ID}->{"Class"}."\n";
        } else {
            print FLOSSprint $count.",".$ID.",".$classTable->{$ID}->{"nreads"}.",,".$classTable->{$ID}->{"Class"}."\n";
        }
        $count++;
    }
    close(FLOSSprint);
    
    #Fill final table out of csv file
    print "Insert FLOSS data from csv file into table TIS_".$analysis_id."_transcripts_FLOSS\n";
    system("sqlite3 -separator , ".$db_results." \".import ".$TMP."/FLOSSclass.csv TIS_".$analysis_id."_transcripts_FLOSS\"")== 0 or die "system failed: $?";
    
    #Remove csv file
    system("rm ".$TMP."/FLOSSclass.csv");
    
    #Remove tmp files
    system("rm -rf ".$TMP."/classification");
    
}

$dbh->disconnect();

# End time
print "   DONE! \n";
my $end = time - $start;
printf("Runtime FLOSS calculation: %02d:%02d:%02d\n\n",int($end/3600), int(($end % 3600)/60), int($end % 60));



############
# THE SUBS #
############

### DBH ###

sub dbh {
    
    # Catch
    my $dsn  =   $_[0];
    my $us	=   $_[1];
    my $pw	=   $_[2];
    
    # Init DB
    
    my $dbh = DBI->connect($dsn,$us,$pw,{ RaiseError => 1 },) || die "Cannot connect: " . $DBI::errstr;
    
    #Return
    return($dbh);
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
    
    $query = "select value from `arguments` where variable = \'igenomes_root\'";
    $sth = $dbh->prepare($query);
	$sth->execute();
	my $igenomes_root = $sth->fetch()->[0];
    
    $query = "select value from `arguments` where variable = \'nr_of_cores\'";
    $sth = $dbh->prepare($query);
	$sth->execute();
	my $nr_of_cores = $sth->fetch()->[0];
    
    $query = "select value from `arguments` where variable = \'ens_db\'";
    $sth = $dbh->prepare($query);
	$sth->execute();
	my $ens_db = $sth->fetch()->[0];

    
    # Return input variables
    return($ensemblversion,$ens_db,$species,$igenomes_root,$nr_of_cores);
    
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
        my $query_analysis_id = "select ID from TIS_overview";
        my @ids = @{$dbh->selectcol_arrayref($query_analysis_id)};
        $idsref = \@ids;
    }
    else {
        my @ids = split(/,/,$ids_in);
        $idsref = \@ids;
        
    }
    
    return($idsref);
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
        #Only nuclear chromosomes
        if($1 ne "MT"){
            push (@chr,$1);
        }
    }
    
    # Get correct coord_system_id
    my $query = "SELECT coord_system_id FROM coord_system where name = 'chromosome' and version = '".$assembly."'";
	my $sth = $dbh->prepare($query);
	$sth->execute();
    @coord_system = $sth->fetchrow_array();
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
        @ids = $sth->fetchrow_array();
        $seq_region_id = $ids[0];
        $chrs->{$_}{'seq_region_id'} = $seq_region_id;
        $sth->finish();
    }
    
    #Disconnect DBH
    $dbh->disconnect();
	
	# Return
	return($chrs);
    
}


### Make Length Distribution ##

sub makeLD{
    
    # Catch
    my $dsnLD  =   $_[0];
    my $usLD  =   $_[1];
    my $pwLD  =   $_[2];
    my $chr  =  $_[3];
    my $seq_region_id  =  $_[4];
    my $fastqname  =  $_[5];
    my $analysis_id  =  $_[6];
    my $dsn_ENS  =   $_[7];
    my $us_ENS  =   $_[8];
    my $pw_ENS  =   $_[9];
    
    #Make connection with results DB
    my $dbh = dbh($dsnLD,$usLD,$pwLD);
    my $dbhENS = dbh($dsn_ENS, $us_ENS, $pw_ENS);
    
    #Define hash to store all information
    my $LD={};
    
    #First, put all RPFdata for chromosome back in a hash
    my $getRPFdata = "select RPF||start, RPF, count, strand, start from count_".$fastqname."_splitRPF where chr='".$chr."';";
    my $sthRPF = $dbh->prepare($getRPFdata);
    $sthRPF->execute();
    my $RPFs = $sthRPF->fetchall_hashref('RPF||start');
    
    
    #Get all TIStranscripts of that chromosome
    my $query = "SELECT tr_stable_id, strand, start, annotation, tr_seq, aa_seq from TIS_".$analysis_id."_transcripts where chr='".$chr."';";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    while( my ($ENStrID, $strand, $TISstart, $annotation, $tr_seq, $aa_seq) = $sth->fetchrow_array()){
        #Search for the normal transcript ID
        my $query2 = "SELECT transcript_id FROM tr_translation WHERE stable_id = '".$ENStrID."';";
        my $sth2 = $dbh->prepare($query2);
        $sth2->execute();
        my $trID = $sth2->fetchrow_array();
        
        #make id for TIS
        my $TisID = $trID."_".$TISstart;
        
        #Get all exons of TIStranscript
        my $query3 = "SELECT exon_id, rank FROM exon_transcript where transcript_id = '".$trID."';";
        my $sth3 = $dbhENS->prepare($query3);
        $sth3->execute();
        my $exonstruct = $sth3->fetchall_hashref('exon_id');
        $sth3->execute();
        my $exonsByRank = $sth3->fetchall_hashref('rank');
        
        #Get coordinates for each exon
        foreach my $exonID (keys %{$exonstruct}){
            my $query4 = "SELECT seq_region_start, seq_region_end from exon where exon_id = '".$exonID."';";
            my $sth4 = $dbhENS->prepare($query4);
            $sth4->execute();
            my @startstop = $sth4->fetchrow_array();
            $exonstruct->{$exonID}->{'start'}=$startstop[0];
            $exonstruct->{$exonID}->{'stop'}=$startstop[1];
        }
        
        #Init
        my $ORFexons={};
        
        #Calculate number of nucleotides in putative CDS
        my $aminoacids = length($aa_seq);
        my $nucleotides = 3*$aminoacids;
        
        #Init
        my $nucleotidesInRemainingORF = $nucleotides;
        my $startExonId;
        my $rankStart;
        
        #Determine start exon
        foreach my $exonID (keys %{$exonstruct}){
            if($TISstart<=$exonstruct->{$exonID}->{'stop'} && $TISstart>=$exonstruct->{$exonID}->{'start'}){
                $startExonId = $exonID;
                $rankStart = $exonstruct->{$exonID}->{'rank'};
            }
        }
        
        my $currentExon = $startExonId;
        my $currentRank = $rankStart;
        
        #Determination of further exon structure is different along the strand
        if($strand == 1){
            #Set TIS as start coordinate of the first exon in the ORF
            $ORFexons->{$startExonId}->{'start'} = $TISstart;
            while($nucleotidesInRemainingORF > ($exonstruct->{$currentExon}->{'stop'} - $ORFexons->{$currentExon}->{'start'} + 1)){ #+1 because it includes the start and stop place
                #As long as the putatively coding sequence still contains more nucleotides than the current exon, it is needed to introduce an extra intron.
                #First, save the end coordinate of the current exon. And correct the remaining nucleotides.
                $ORFexons->{$currentExon}->{'stop'} = $exonstruct->{$currentExon}->{'stop'};
                $nucleotidesInRemainingORF = $nucleotidesInRemainingORF - ($ORFexons->{$currentExon}->{'stop'} - $ORFexons->{$currentExon}->{'start'} + 1);
                #Go to the next exon
                $currentRank = $currentRank + 1;
                $currentExon = $exonsByRank->{$currentRank}->{'exon_id'};
                #Save the start coordinate of this next exon.
                $ORFexons->{$currentExon}->{'start'} = $exonstruct->{$currentExon}->{'start'};
            }
            #After this loop, it is the last exon of the ORF. Determine stop coordinate
            $ORFexons->{$currentExon}->{'stop'} = $exonstruct->{$currentExon}->{'start'} + $nucleotidesInRemainingORF - 1;
        } elsif($strand == -1){
            #start and stop in exon table are not dependant of the strand
            #Set TIS as stop coordinate of the first exon in the ORF
            $ORFexons->{$startExonId}->{'stop'} = $TISstart;
            while($nucleotidesInRemainingORF > ($ORFexons->{$currentExon}->{'stop'} - $exonstruct->{$currentExon}->{'start'} + 1)){
                $ORFexons->{$currentExon}->{'start'} = $exonstruct->{$currentExon}->{'start'};
                $nucleotidesInRemainingORF = $nucleotidesInRemainingORF - ($ORFexons->{$currentExon}->{'stop'} - $ORFexons->{$currentExon}->{'start'} + 1);
                $currentRank = $currentRank + 1; #the most right exon has the lowest rank (in contrast with the start and stop of the exons itselves)
                $currentExon = $exonsByRank->{$currentRank}->{'exon_id'};
                $ORFexons->{$currentExon}->{'stop'} = $exonstruct->{$currentExon}->{'stop'};
            }
            $ORFexons->{$currentExon}->{'start'} = $exonstruct->{$currentExon}->{'stop'} - $nucleotidesInRemainingORF + 1;
        }
        #In ORFexons are now the start and stop coordinate of every part of the ORF, splitted in exons
        
       
        #Search now for all RPFs in hash that fall in exon structure
        foreach my $exonID (keys %{$ORFexons}){
            for (my $startRPF = $ORFexons->{$exonID}->{'start'}; $startRPF<=$ORFexons->{$exonID}->{'stop'}; $startRPF++){
                for (my $RPFlength = 26; $RPFlength<=34; $RPFlength++){
                    my $RPFid =$RPFlength.$startRPF;
                    if($RPFs->{$RPFid}->{'count'}){
                        if ($RPFs->{$RPFid}->{'strand'} == $strand){
                            if ($LD->{$TisID}->{"RPF"}->{$RPFlength}){
                                $LD->{$TisID}->{"RPF"}->{$RPFlength}=$LD->{$TisID}->{"RPF"}->{$RPFlength}+$RPFs->{$RPFid}->{'count'};
                            } else {
                                $LD->{$TisID}->{"RPF"}->{$RPFlength}=$RPFs->{$RPFid}->{'count'};
                            }
                            if ($LD->{$TisID}->{"nreads"}){
                                $LD->{$TisID}->{"nreads"} = $LD->{$TisID}->{"nreads"} + $RPFs->{$RPFid}->{'count'};
                            } else {
                                $LD->{$TisID}->{"nreads"} = $RPFs->{$RPFid}->{'count'};
                            }
                        }
                    }
                }
            }
        }

    }
    
    $dbh->disconnect();
    $dbhENS->disconnect();
    
    return($LD);
    
}


### Make reference histogram for that chromosome ###

sub refHist{
    
    #Catch
    my $dsn_results  =  $_[0];
    my $us_results  =  $_[1];
    my $pw_results  =  $_[2];
    my $dsn_ENS  =   $_[3];
    my $us_ENS  =   $_[4];
    my $pw_ENS  =   $_[5];
    my $cores  =  $_[6];
    my $chrs  =  $_[7];
    my $fastqname  =  $_[8];
    
    #Make connection with results DB
    my $dbh = dbh($dsn_results,$us_results,$pw_results);
    
    # Check if table already exists
    my $table_name = "FLOSS_ref_fractions";
    my $query = "select name from sqlite_master where type='table' and name like '$table_name'";
    my $sth = $dbh->prepare($query);
    my $table_exists = 0;
    
    if($sth->execute()){
        while( my $t = $sth->fetchrow_array() ) {
            if( $t =~ /^$table_name$/ ){
                $table_exists = 1;
                print "     ...FLOSS reference table already exists. \n";
                
                print "Opening reference fractions\n\n\n";
                #get reference fractions out of DB and store in hash
                $query = "SELECT * FROM FLOSS_ref_fractions";
                $sth = $dbh->prepare($query);
                $sth->execute();
                my $refFracs = $sth->fetchall_hashref('RPF');
                
                return ($refFracs);
            }
        }
    }
    
    # If table does not exist, create it
    if($table_exists == 0){
        
        print "Creating new FLOSS reference table. \n";
        
        #init
        my $refRPFcountTot={};
        my $refRPFfraction={};
        for(my $i=26;$i<=34;$i++){
            $refRPFcountTot->{"$i"}=0;
        }
        
        #To speed up search in RPF database, add indexes
        my $query_idx = "CREATE INDEX IF NOT EXISTS count_".$fastqname."_splitRPF_id_idx ON count_".$fastqname."_splitRPF (chr)";
        $dbh->do($query_idx);
        $query_idx = "CREATE INDEX IF NOT EXISTS count_".$fastqname."_splitRPF_id_idx ON count_".$fastqname."_splitRPF (start)";
        $dbh->do($query_idx);
        $query_idx = "CREATE INDEX IF NOT EXISTS count_".$fastqname."_splitRPF_id_idx ON count_".$fastqname."_splitRPF (strand)";
        $dbh->do($query_idx);
        
        # Attempt to create table
        my $create_table = "CREATE TABLE `FLOSS_ref_fractions` (`RPF` INTEGER PRIMARY KEY,`fraction` DECIMAL(6,5) );";
        $dbh->do($create_table);
        
        print "Initiating multi core for construction of reference histogram \n";
        # Init multi core
        my $pm = new Parallel::ForkManager($cores);
        print "   Using ".$cores." core(s)\n   ---------------\n";
        
        foreach my $chr(sort keys %{$chrs}){

            ### Start parallel process
            $pm->start and next;
        
            #Init hash for length distribution for cutoff calculation for each chromosome
            my $refRPFcountChr = {};
            my $codTranscriptLDChr={};
            
            ### DBH per process
            my $dbhENS = dbh($dsn_ENS,$us_ENS,$pw_ENS);
            my $dbh = dbh($dsn_results,$us_results,$pw_results);
            
            #First, put all RPFdata for chromosome back in a hash
            my $getRPFdata = "select RPF||start, RPF, count, strand, start from count_".$fastqname."_splitRPF where chr='".$chr."';";
            my $sthRPF = $dbh->prepare($getRPFdata);
            $sthRPF->execute();
            my $RPFs = $sthRPF->fetchall_hashref('RPF||start');
        
            #Determine seq region id for chormosome
            my $seq_region_id = $chrs->{$chr}{'seq_region_id'};
            
            ## Get all protein coding genes in chromosome from Ensembl table gene
            my $query = "SELECT gene_id, seq_region_start, seq_region_end FROM gene where seq_region_id = '".$seq_region_id."' AND biotype = 'protein_coding';";
            my $sth = $dbhENS->prepare($query);
            $sth->execute();
            my $genes = $sth->fetchall_hashref('gene_id');
                my $aantalGenen=scalar(keys %{$genes});
            
            # Get all non coding transcripts in chromosome from Ensembl table transcript
            $query = "SELECT transcript_id, seq_region_start, seq_region_end FROM transcript where seq_region_id = '".$seq_region_id."' AND (biotype = 'IG_V_pseudogene' OR biotype = 'Mt_rRNA' OR biotype = 'Mt_tRNA' OR biotype = 'TR_V_pseudogene' OR biotype = 'lincRNA' OR biotype = 'miRNA' OR biotype = 'misc_RNA' OR biotype = 'polymorphic_pseudogene' OR biotype = 'pseudogene' OR biotype = 'rRNA' OR biotype = 'snRNA' OR biotype = 'snoRNA' OR biotype = 'transcribed_unprocessed_pseudogene' OR biotype = 'transcribed_processed_pseudogene' OR biotype = 'translated_processed_pseudogene' OR biotype = 'unitary_pseudogene' OR biotype = 'unprocessed_pseudogene');";
            $sth = $dbhENS->prepare($query);
            $sth->execute();
            my $ncTranscripts = $sth->fetchall_hashref('transcript_id');
            
            #Run over genes
            foreach my $gene_id (keys %{$genes}){
                my $overlap = 'False';
                
                #Check overlap of gene with other non-coding transcripts
                foreach my $nc (keys %{$ncTranscripts}){
                    $overlap = 'False';
                    # transcript               |--------|    or             |----|
                    # gene                |---------|                 |-----------------|
                    if ($genes->{$gene_id}->{'seq_region_start'}<=$ncTranscripts->{$nc}->{'seq_region_start'} && $genes->{$gene_id}->{'seq_region_end'}>=$ncTranscripts->{$nc}->{'seq_region_start'}){
                        $overlap = 'True';
                    }
                    # transcript      |----------|
                    # gene                |---|
                    elsif ($genes->{$gene_id}->{'seq_region_start'}>=$ncTranscripts->{$nc}->{'seq_region_start'} && $genes->{$gene_id}->{'seq_region_end'}<=$ncTranscripts->{$nc}->{'seq_region_end'}){
                        $overlap = 'True';
                    }
                    # transcript         |-----------|
                    # gene                       |----------|
                    elsif ($genes->{$gene_id}->{'seq_region_start'}<=$ncTranscripts->{$nc}->{'seq_region_end'} && $genes->{$gene_id}->{'seq_region_end'}>=$ncTranscripts->{$nc}->{'seq_region_end'}){
                        $overlap = 'True';
                    }
                    
                    #leave loop if there is overlap
                    last if ($overlap eq 'True');
                    
                }
                
                if($overlap eq 'False'){
                    
                    ## Get canonical transcript from ENS db
                    my $canonical_transcript = get_canonical_transcript($dsn_ENS,$us_ENS,$pw_ENS,$gene_id);
                    
                    #Run over transcript of that gene_id (in fact only one id)
                    foreach my $canonical_transcript_id (keys %{$canonical_transcript}){
                        
                        #Get translation data
                        my ($trans,$exons)  =   get_translation_data($dsn_ENS,$us_ENS,$pw_ENS,$canonical_transcript_id,$seq_region_id);
                        my $strand          =   $canonical_transcript->{$canonical_transcript_id}->{'seq_region_strand'};
                        
                        #Check if transcript has translation info and chip off 5' en 3'UTR
                        if($trans->[0]){
                            my $startFirstExon = ($strand eq '1') ? $exons->{$trans->[0]}{'seq_region_start'} + ${$trans}[2] - 1 : $exons->{$trans->[0]}{'seq_region_end'} - ${$trans}[2] + 1;
                            my $stopLastExon =  ($strand eq '1') ? $exons->{$trans->[1]}{'seq_region_start'} + ${$trans}[3] - 1 : $exons->{$trans->[1]}{'seq_region_end'} - ${$trans}[3] + 1;
                            $exons->{$trans->[0]}{'seq_region_start'}=$startFirstExon;
                            $exons->{$trans->[1]}{'seq_region_end'}=$stopLastExon;
                            
                            
                            
                            #Search now for all RPFs in hash that fall in exon structure
                            foreach my $exonID (keys %{$exons}){
                                for (my $startRPF = $exons->{$exonID}->{'seq_region_start'}; $startRPF<=$exons->{$exonID}->{'seq_region_end'}; $startRPF++){
                                    for (my $RPFlength = 26; $RPFlength<=34; $RPFlength++){
                                        my $RPFid =$RPFlength.$startRPF;
                                        if($RPFs->{$RPFid}->{'count'}){
                                            if($RPFs->{$RPFid}->{'strand'} == $strand){
                                                if ($refRPFcountChr->{$RPFlength}){
                                                    $refRPFcountChr->{$RPFlength} = $refRPFcountChr->{$RPFlength} + $RPFs->{$RPFid}->{'count'};
                                                } else {
                                                    $refRPFcountChr->{$RPFlength} = $RPFs->{$RPFid}->{'count'};
                                                }
                                                if ($codTranscriptLDChr->{$canonical_transcript_id}->{$RPFlength}){
                                                    $codTranscriptLDChr->{$canonical_transcript_id}->{$RPFlength} = $codTranscriptLDChr->{$canonical_transcript_id}->{$RPFlength} + $RPFs->{$RPFid}->{'count'};
                                                } else {
                                                    $codTranscriptLDChr->{$canonical_transcript_id}->{$RPFlength} = $RPFs->{$RPFid}->{'count'};
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            
                        }
                        
                    }
                    
                }
                
            }
            
            #write count to file
            system("mkdir -p ".$TMP."/RPFcount/");
            open (TMPcount,">".$TMP."/RPFcount/RPFcount_".$fastqname."_".$chr.".txt") || die "ERROR opening RPF tmp file to write \n";
            foreach my $RPFlength (keys %{$refRPFcountChr}){
                print TMPcount $RPFlength."\t".$refRPFcountChr->{$RPFlength}."\n";
            }
            
            close(TMPcount);
            
            #write coding transcript LD to file for cutoff calculation
            system("mkdir -p ".$TMP."/cutoff/");
            open (cutoffLD,">".$TMP."/cutoff/cutoffLD_".$fastqname."_".$chr.".txt") || die "ERROR opening cutoff LD tmp file to write \n";
            foreach my $transcriptID (keys %{$codTranscriptLDChr}){
                foreach my $RPFlength (keys %{$codTranscriptLDChr->{$transcriptID}}){
                    print cutoffLD $transcriptID."\t".$RPFlength."\t".$codTranscriptLDChr->{$transcriptID}->{$RPFlength}."\n";
                }
            }
            
            close(cutoffLD);
            print "* Finished chromosome ".$chr."\n";
            
            ### Finish child
            $pm->finish;
        }
        
        #Wait all Children
        $pm->wait_all_children();
        
        print "\nCalculate total counts and fractions for reference table. \n";
        #Calculate total counts over all chromosomes
        foreach my $chr(sort keys %{$chrs}){
            open(RPFread,"<".$TMP."/RPFcount/RPFcount_".$fastqname."_".$chr.".txt") || die "ERROR reading RPF tmp file \n";
            while(my $read = <RPFread>){
                chomp($read);
                my @RPFchr = split("\t",$read);
                $refRPFcountTot->{$RPFchr[0]}=$refRPFcountTot->{$RPFchr[0]}+$RPFchr[1];
            }
        }
    
        #Calculate fractions and store in reference table
        my $sumCounts=0;
        foreach my $iRPF (keys %{$refRPFcountTot}){
            $sumCounts=$refRPFcountTot->{$iRPF}+$sumCounts;
            }
        foreach my $iRPF (keys %{$refRPFcountTot}){
            $refRPFfraction->{$iRPF}->{"fraction"} = $refRPFcountTot->{$iRPF}/$sumCounts;
            $query = "INSERT INTO FLOSS_ref_fractions (RPF, fraction) VALUES (".$iRPF.", ".$refRPFfraction->{$iRPF}->{"fraction"}.")";
            $sth = $dbh->prepare($query);
            $sth->execute();;
        }
        
        #Remove tmp files
        system("rm -rf ".$TMP."/RPFcount");
        print "Creating FLOSS reference table completed \n\n\n";
        
        $dbh->disconnect();
        return ($refRPFfraction);
    }
    
}

### get gene_id transcripts

sub get_canonical_transcript{
    
    #Catch
    my $dsn_ENS     =   $_[0];
    my $us_ENS      =   $_[1];
    my $pw_ENS      =   $_[2];
    my $gene_id     =   $_[3];
    
    #Init
    my $dbh_ENS = dbh($dsn_ENS,$us_ENS,$pw_ENS);
    my $trs = {};
    
    #Get canonical transcript id
    my $query = "SELECT canonical_transcript_id FROM gene where gene_id = '".$gene_id."';";
    my $sth =$dbh_ENS->prepare($query);
    $sth->execute();
    my @canTranscriptID = $sth->fetchrow_array();
    
    #Get transcript info
    $query = "SELECT transcript_id, gene_id, seq_region_id, seq_region_strand, seq_region_start, seq_region_end, biotype, stable_id, canonical_translation_id from transcript where transcript_id = '".$canTranscriptID[0]."';";
    $sth = $dbh_ENS->prepare($query);
    $sth->execute();
    $trs = $sth->fetchall_hashref('transcript_id');
    
    #Return
    return($trs);
    
}

### Get transcript translation data ###

sub get_translation_data {
    
    # Catch
	my $dsn_ENS         =   $_[0];
    my $us_ENS          =   $_[1];
    my $pw_ENS          =   $_[2];
    my $tr_id           =   $_[3];
    my $seq_region_id   =   $_[4];
    
    #Init
    my $dbh_ENS = dbh($dsn_ENS,$us_ENS,$pw_ENS);
    my $exons = {};
    my $trans = [];
    my $exonstruct = {};
    my $codingExons = [];
    
    # Get translation data for transcrip_id if available
    my $query = "SELECT start_exon_id,end_exon_id,seq_start,seq_end FROM translation where transcript_id = '".$tr_id."'";
    my $sth = $dbh_ENS->prepare($query);
	$sth->execute();
    $trans = $sth->fetchrow_arrayref();
    
    #Get exon organisation for transcript_id
    $query = "SELECT exon_id, rank FROM exon_transcript where transcript_id = '".$tr_id."';";
    $sth = $dbh_ENS->prepare($query);
    $sth->execute();
    $exonstruct = $sth->fetchall_hashref('exon_id');
    
    
    #Get the exons with translation info
    my $rankstart = $exonstruct -> {$trans->[0]} -> {'rank'};
    my $rankend = $exonstruct -> {$trans->[1]} -> {'rank'};
    foreach my $exonID (keys %{$exonstruct}){
        if($exonstruct->{$exonID}->{'rank'}>=$rankstart && $exonstruct->{$exonID}->{'rank'}<=$rankend){
            push( @$codingExons, $exonID);
        }
    }

    my $statementExons = join( ', ',@$codingExons);
    # Get exon info from ENS DB for all coding exons
    $query = "SELECT exon_id, seq_region_start, seq_region_end FROM exon where exon_id IN (".$statementExons.") AND seq_region_id = '".$seq_region_id."'";
    $sth = $dbh_ENS->prepare($query);
    $sth->execute();
    $exons = $sth->fetchall_hashref('exon_id');
    
        
    #Return
    return($trans,$exons);
    
}


### Calculate cutoff values

sub cutoff {
    
    #Catch
    my $dsn_results  =  $_[0];
    my $us_results  =  $_[1];
    my $pw_results  =  $_[2];
    my $chrs  =  $_[3];
    my $fastqname  =  $_[4];
    my $db_results  =  $_[5];
    
    #Make connection with results DB
    my $dbh = dbh($dsn_results,$us_results,$pw_results);
    
    # Check if table already exists
    my $table_name = "FLOSS_cutoff";
    my $query = "select name from sqlite_master where type='table' and name like '$table_name'";
    my $sth = $dbh->prepare($query);
    my $table_exists = 0;
    
    if($sth->execute()){
        while( my $t = $sth->fetchrow_array() ) {
            if( $t =~ /^$table_name$/ ){
                $table_exists = 1;
                print "     ...FLOSS cutoff table already exists. \n";
                
                print "Getting data from cutoff table \n\n\n";
                #get cutoff data out of DB and store in hash
                $query = "SELECT * FROM ".$table_name.";";
                $sth = $dbh->prepare($query);
                $sth->execute();
                my $smoother={};
                $smoother->{"nreads"} = $sth->fetchall_hashref('nreads');
                $query = "SELECT MIN(nreads) FROM ".$table_name.";";
                $sth = $dbh->prepare($query);
                $sth->execute();
                my ($smootherMin) = $sth->fetchrow_array();
                $query = "SELECT MAX(nreads) FROM ".$table_name.";";
                $sth = $dbh->prepare($query);
                $sth->execute();
                my ($smootherMax) = $sth->fetchrow_array();
                $smoother->{"min"} = $smootherMin;
                $smoother->{"max"} = $smootherMax;
                return($smoother);
            }
        }
    }
    
    # If table does not exist, create it
    if($table_exists == 0){
        print "Creating new FLOSS cutoff table. \n";
        # Attempt to create table
        my $create_table = "CREATE TABLE ".$table_name." (`nreads` INTEGER PRIMARY KEY,`score` VARCHAR(20) );";
        $dbh->do($create_table);
        
        print "Opening reference fractions for cutoff calculation \n";
        #get reference fractions out of DB and store in hash
        my $query = "SELECT * FROM FLOSS_ref_fractions";
        $sth = $dbh->prepare($query);
        $sth->execute();
        my $refFracs = $sth->fetchall_hashref('RPF');
        
        print "Opening length distribution of coding transcripts \n";
        #Get length distribution of coding transcripts out of sub refHist
        my $cutoffLDTot = {};
        foreach my $chr(sort keys %{$chrs}){
            open(cutoffLDread,"<".$TMP."/cutoff/cutoffLD_".$fastqname."_".$chr.".txt") || die "ERROR reading cutoff LD tmp file \n";
            while(my $read = <cutoffLDread>){
                chomp($read);
                my @cutoffLDchr = split("\t",$read);
                if ($cutoffLDTot->{$cutoffLDchr[0]}->{"RPF"}->{$cutoffLDchr[1]}){
                    $cutoffLDTot->{$cutoffLDchr[0]}->{"RPF"}->{$cutoffLDchr[1]} = $cutoffLDTot->{$cutoffLDchr[0]}->{"RPF"}->{$cutoffLDchr[1]} + $cutoffLDchr[2];
                } else {
                    $cutoffLDTot->{$cutoffLDchr[0]}->{"RPF"}->{$cutoffLDchr[1]} = $cutoffLDchr[2];
                }
                if ($cutoffLDTot->{$cutoffLDchr[0]}->{"nReads"}){
                    $cutoffLDTot->{$cutoffLDchr[0]}->{"nReads"} = $cutoffLDTot->{$cutoffLDchr[0]}->{"nReads"} + $cutoffLDchr[2];
                } else {
                    $cutoffLDTot->{$cutoffLDchr[0]}->{"nReads"} = $cutoffLDchr[2];
                }
            }
        }
        
        print "Calculating FLOSS of coding transcripts \n";
        
        #Calculate FLOSS for each coding transcript ID
        foreach my $transcriptID (keys %{$cutoffLDTot}){
            $cutoffLDTot->{$transcriptID}->{"FLOSS"} = flossCalc($transcriptID, $cutoffLDTot, $refFracs);
        }
        
        print "Calculating cutoff statistics \n";
        #Define rolling window width (cfr. Ingolia et al. 2014)
        my $rollWidth = 200;
        
        #Store data for R in csv file
        print "export R input data in csv file \n";
        open toR, ">>".$TMP."/toR.csv" or die "ERROR exporting data to R";
        print toR "transcriptid,nreads,score\n";
        foreach my $transcriptID (keys %{$cutoffLDTot}){
            print toR $transcriptID.",".$cutoffLDTot->{$transcriptID}->{"nReads"}.",".$cutoffLDTot->{$transcriptID}->{"FLOSS"}."\n";
        }
        close(toR);
        
        print "initiate bridge to R and define input variables \n";
        # Create a communication bridge with R and start R
        my $R = Statistics::R->new();
        $R->set('rollwidth', $rollWidth) || die $!;
        my $path = $TMP."/toR.csv";
        my $pathOutput = $TMP."/toPerl.csv";
        $R->set('file', $path) || die $!;
        $R->set('fileOutput', $pathOutput) || die $!;
        $R->set('RlibLocation', $RlibLocation) || die $!;
        
        print "Run R script \n";
        
            $R->run(q`library("zoo", lib.loc=RlibLocation)`);
        
            #Read perl data in
            $R->run(q`unsorted <- read.csv(file, header = TRUE, sep = ",")`);
        
            #Sort for amount of reads
            $R->run(q`fr <- unsorted[order(unsorted$nreads),]`);
            #Calculate statistics
            $R->run(q`nr <- rollapply(fr$nreads, width=rollwidth, function(xs) { mean(xs, na.rm=T) })`);
            $R->run(q`medsc <- rollapply(fr$score, width=rollwidth, function(xs) { median(xs, na.rm=T) })`);
            $R->run(q`q1 <- rollapply(fr$score, width=rollwidth, function(xs) { quantile(xs, probs=0.25, na.rm=T) })`);
            $R->run(q`q3 <- rollapply(fr$score, width=rollwidth, function(xs) { quantile(xs, probs=0.75, na.rm=T) })`);
            #Use Tukey's method
            $R->run(q`rawExtreme <- q3 + 3.0*(q3-q1)`);
            #Smoothing
            $R->run(q`extLoess <- loess(extreme ~ nread, data.frame(nread=log10(nr), extreme=rawExtreme))`);
        
            $R->run(q`testreads <-log10(seq(ceiling(min(nr)),max(nr),1))`);
            $R->run(q`predictions <- predict(extLoess, testreads)`);
            $R->run(q`smootherMin <- min(round(10^testreads))`);
            $R->run(q`smootherMax <- max(round(10^testreads))`);
            #Make data frame to export to perl
            $R->run(q`toPerl <- data.frame(round(10^testreads, digits=0),sprintf("%.14f", round(predictions, digits=14)))`);
            $R->run(q` colnames(toPerl) <- c("nreads", "predictions")`);
        
            #write to csv file
            $R->run(q`write.table(toPerl, file = fileOutput, row.names=FALSE, col.names=FALSE, sep=',')`);
        
            my $smootherMin = $R->get('smootherMin');
            my $smootherMax = $R->get('smootherMax');
        $R->stop();
        
        #get data back from R
        print "Get data back from R \n";
        open fromR, "<".$TMP.'/toPerl.csv' or die "ERROR opening file with data coming back from R!";
        my $smoother={};
        $smoother->{"min"} = $smootherMin;
        $smoother->{"max"} = $smootherMax;
        while(my $read = <fromR>){
            chomp($read);
            my @read=split(",",$read);
            $smoother->{"nreads"}->{$read[0]}->{'nreads'} = $read[0];
            $smoother->{"nreads"}->{$read[0]}->{'score'} = $read[1];
        }
        
        #Insert into cutoff table (for later analyses with the same cutoff table)
        print "Insert data into cutoff table\n";
        system("sqlite3 -separator , ".$db_results." \".import ".$TMP."/toPerl.csv ".$table_name."\"")== 0 or die "system failed: $?";
        
        #Remove data input file for R
        system("rm ".$TMP."/toR.csv");
        system("rm ".$TMP."/toPerl.csv");

        print "\nCutoff calculated \n\n\n";
        
        #Remove data coming from refHist calculation to cutoff calculation
        system("rm -rf ".$TMP."/cutoff");
        
        return ($smoother);
        
    }
    
}

### Calculate FLOSS score for a given id

sub flossCalc {
    
    #Catch
    my $id  =  $_[0];
    my $LD  =  $_[1];
    my $refFractions  =  $_[2];
    
    #Calculate fractions of given id
    my $sumCounts=0;
    my $RPFfractions={};
    foreach my $RPFlength (keys %{$LD->{$id}->{"RPF"}}){
        $sumCounts=$sumCounts+$LD->{$id}->{"RPF"}->{$RPFlength};
    }
    foreach my $RPFlength (keys %{$LD->{$id}->{"RPF"}}){
        $RPFfractions->{$RPFlength}=$LD->{$id}->{"RPF"}->{$RPFlength}/$sumCounts;
    }
    
    #Calculate sum of the absolute difference between fractions for id and ref fractions for each RPF length
    my $sumAbsDiff=0;
    foreach my $RPFlength (keys %{$LD->{$id}->{"RPF"}}){
        $sumAbsDiff =$sumAbsDiff + abs($RPFfractions->{$RPFlength} - $refFractions->{$RPFlength}->{"fraction"});
        
    }
    
    #Calculate floss
    my $floss = $sumAbsDiff / 2;
    
    return($floss);
    
}

### Calculate floss score for all sORFs and classify them for coding potential

sub classify {
    
    #Catch
    my $chr  =  $_[0];
    my $fastqname  =  $_[1];
    my $smoother  =  $_[2];
    my $LD  =  $_[3];
    my $refFracs  =  $_[4];
    
    
    #init hash
    my $score = {};
    
    #Calculate Floss score and classify
    foreach my $ID (keys %{$LD}){
        if(defined $LD->{$ID}->{"nreads"} && $LD->{$ID}->{"nreads"}>0){
            my $nreads = $LD->{$ID}->{"nreads"};
            $score->{$ID}->{"FLOSS"} = flossCalc($ID, $LD, $refFracs);
            #Predict extreme out of Loess smoother out of the amount of reads for that sORF
            if ($nreads<$smoother->{"min"}){
                $score->{$ID}->{"Class"} = "Not in cutoff range";
            }  elsif ($nreads>$smoother->{"max"}){
                $score->{$ID}->{"Class"} = "Not in cutoff range";
            } else {
                #Classify
                my $cutoffBorder = $smoother->{"nreads"}->{$nreads}->{'score'};
                if ($cutoffBorder =~ /^\"(.+)\"$/){
                    $cutoffBorder = $1;
                }
                if ($cutoffBorder eq "NA"){
                    next;
                } elsif ($score->{$ID}->{"FLOSS"}<$cutoffBorder){
                    $score->{$ID}->{"Class"} = "Good"; #Means probably coding
                } else {
                    $score->{$ID}->{"Class"} = "Extreme"; #Means probably non-coding
                }
            }
        } else {
            $LD->{$ID}->{"nreads"}=0;
            my $nreads = $LD->{$ID}->{"nreads"};
            $score->{$ID}->{"FLOSS"}=undef;
            $score->{$ID}->{"Class"} = "No reads";
        }
    }
    
    #write to temporary file (because it is for each chromosome)
    system("mkdir -p ".$TMP."/classification/");
    open classificationFH,">".$TMP."/classification/classification_".$fastqname."_".$chr.".txt" || die "ERROR opening classification tmp file to write \n";
    foreach my $ID (keys %{$score}){
        if (defined $score->{$ID}->{"FLOSS"}){
            print classificationFH $ID."\t".$LD->{$ID}->{"nreads"}."\t".$score->{$ID}->{"FLOSS"}."\t".$score->{$ID}->{"Class"}."\n";
        } else {
            print classificationFH $ID."\t".$LD->{$ID}->{"nreads"}."\tundef\t".$score->{$ID}->{"Class"}."\n";
        }
    }
    close(classificationFH);
    
    return;
}





