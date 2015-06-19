#!/usr/bin/perl -w
use strict;
use warnings;
use DBI;
use DBD::SQLite;
use Data::Dumper;
use Parallel::ForkManager;
use Getopt::Long;
use Storable 'dclone';
use Cwd;

#################
# Perl commands #
#################

## Command line ##

# ./TIScalling_intergenic.pl --sqlite SQLite/results_STAR12.db --cores 22 --local_max 1 --min_count 10 --igenomes_root /data/igenomes/Mus_musculus/Ensembl/GRCm38/ --ens_db /data/jeroenc/ribo_prof_micropeptide_pipeline/SQLite/ENS_mmu_72.db

## Galaxy ##

#####Check###### ./TIScalling_categorised.pl --sqlite_db "${sqlite_db}" --ens_db "${ensembl_db.fields.path}" --cores "${cores}" --igenomes_root "${igenomes_root.fields.path}" --out_sqlite "${out_sqlite}"

# get the command line arguments
my ($work_dir,$local_max,$IGENOMES_ROOT,$out_sqlite,$cores,$tmpfolder,$min_count,$ens_db,$sqlite_db);

GetOptions(
"dir=s"=>\$work_dir,                            # Path to the working directory                                                                 optional argument
"cores=i"=>\$cores,                             # Number of cores to be used                                                                    mandatory argument  # -> (default = number of chromosomes)
"tmp:s" =>\$tmpfolder,                          # Folder where temporary files are stored,                                                      optional  argument (default = $TMP env setting)
"igenomes_root:s" =>\$IGENOMES_ROOT,            # Igenomes root folder                                                                          mandatory argument
"sqlite_db:s" =>\$sqlite_db,                    # SQLite results db with mapping and tr_translation tables                                      mandatory argument
"ens_db:s" =>\$ens_db,                          # SQLite Ensembl DB                                                                             mandatory argument
"local_max:i" =>\$local_max,                    # The range wherein the localmax is (e.g.: 1 means +/- one triplet)                             optional argument (default 1)
"min_count:i" =>\$min_count,                    # The minimum count of riboseq profiles mapping to the aTIS site                                optional argument (default 5)
"out_sqlite:s" =>\$out_sqlite                   # Galaxy specific history file location                                                         Galaxy specific
);

my $CWD             = getcwd;
my $HOME            = $ENV{'HOME'};
$IGENOMES_ROOT      = ($ENV{'IGENOMES_ROOT'}) ? $ENV{'IGENOMES_ROOT'} : $IGENOMES_ROOT;
my $TMP             = ($ENV{'TMP'}) ? $ENV{'TMP'} : ($tmpfolder) ? $tmpfolder : "$CWD/tmp" ; # First select the TMP environment variable, second select the $tmpfolder variable, lastly select current_working_dir/tmp
print "The following tmpfolder is used                          : $TMP\n";
print "The following Ensembl db folder is used                  : $ens_db\n";
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
    #die "\nDon't forget to pass the working directory using the --dir or -d argument!\n\n";
}
if ($IGENOMES_ROOT){
    print "The igenomes_root folder used is                         : $IGENOMES_ROOT\n";
} else {
    die "\nDon't forget to pass the igenomes_root folder --igenomes_root or -ig argument!\n\n";
}
if ($cores){
    print "Number of cores to use for Mapping                       : $cores\n";
} else {
    die "\nDon't forget to pass number of cores to use for mapping using the --cores or -c argument!\n\n";
}
if ($local_max){
    print "The regio for local max is set to                        : $local_max\n";
} else {
    #Choose default value for local_max
    $local_max = 1;
}
if ($min_count){
    print "The minimum coverage of an intergenic TIS is             : $min_count\n";
} else {
    #Choose default value for min_count
    $min_count = 10;
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
my ($ensemblversion,$species) = get_arguments($dsn_results,$us_results,$pw_results);

#Conversion for species terminology
my $spec = ($species eq "mouse") ? "Mus_musculus" : ($species eq "human") ? "Homo_sapiens" : "";
my $spec_short = ($species eq "mouse") ? "mmu" : ($species eq "human") ? "hsa" : "";

# Sqlite Ensembl
my $db_ENS  = $ens_db;
my $dsn_ENS = "DBI:SQLite:dbname=$db_ENS";
my $us_ENS  = "";
my $pw_ENS  = "";

#Old mouse assembly = NCBIM37, new one is GRCm38
my $assembly = ($species eq "mouse" && $ensemblversion >= 70 ) ? "GRCm38"
: ($species eq "mouse" && $ensemblversion < 70 ) ? "NCBIM37"
: ($species eq "human") ? "GRCh37"
: "";

my $chrom_file = $IGENOMES_ROOT."/Annotation/Genes/ChromInfo.txt";
my $BIN_chrom_dir = $IGENOMES_ROOT."/Sequence/Chromosomes_BIN";

#############################################
## Peak Calling on exon_covered transcripts #
##                                          #
#############################################

# Start time
my $start = time;

print "\nGet analysis id...";
## Get Analysis ID
my $id = get_id($dsn_results,$us_results,$pw_results,$local_max,$min_count);
print "                                       : ".$id." \n";

print "Get chromosomes... \n";
## Get chromosomes based on seq_region_id ##
my $chrs = get_chrs($dsn_ENS,$us_ENS,$pw_ENS,$chrom_file,$assembly);

# Create binary chromosomes if they don't exist
print "Checking/Creating binary chrom files ...\n";
if (!-d "$BIN_chrom_dir") {
    create_BIN_chromosomes($BIN_chrom_dir,$cores,$chrs,$work_dir);
}

## Create intergenes
print "Creating intergene table... \n";
create_intergenes($chrom_file,$dsn_ENS,$us_ENS,$pw_ENS,$chrs,$cores,$TMP,$ens_db);

## Tiscalling for all ribo_reads overlapping ENS transcripts
print " \n  Starting transcript peak analysis per chromosome using ".$cores." cores...\n";    
TIS_calling($chrs,$dsn_ENS,$us_ENS,$pw_ENS,$cores,$id,$work_dir,$TMP);

# End time
print "   DONE! \n";
my $end = time - $start;
printf("Runtime TIS calling: %02d:%02d:%02d\n\n",int($end/3600), int(($end % 3600)/60), int($end % 60));


############
# THE SUBS #
############

### TIS CALLING ###

sub TIS_calling{
    
    #Catch
    my $chrs            =   $_[0];
    my $db_ENS          =   $_[1];
    my $us_ENS          =   $_[2];
    my $pw_ENS          =   $_[3];
    my $cores           =   $_[4];
    my $id              =   $_[5];
    my $work_dir        =   $_[6];
    my $TMP             =   $_[7];
    
    # Init multi core
    my $pm = new Parallel::ForkManager($cores);
    
    ## Loop over all chromosomes
    foreach my $chr (sort keys %{$chrs}){
        
        ### Start parallel process
        $pm->start and next;
        
        ### DBH per process
        my $dbh = dbh($dsn_ENS,$us_ENS,$pw_ENS);
        
        ### Open File per process
        my $table = "TIS_intergenic_".$id;
        open TMP, "+>> ".$TMP."/".$table."_".$chr."_tmp.csv" or die $!;
    
        # seq_region_id per chromosome
        my $seq_region_id = $chrs->{$chr}{'seq_region_id'};
        
        # Get intergenes and all reads
        my($intergenes,$CHX_for,$LTM_for,$CHX_rev,$LTM_rev) = get_intergenes_and_reads($seq_region_id,$chr,$dbh);
        
        # Match reads to intergenes
        $intergenes = match_reads_to_intergenes($intergenes,$CHX_for,$CHX_rev,$LTM_for,$LTM_rev);
        
        # Loop over intergenes
        foreach my $ig_id (sort {$a <=> $b} keys %{$intergenes}){
            
            #Get intergene specific reads
            my($LTM_reads,$CHX_reads) = get_intergene_overlapping_reads($dbh,$ig_id,$intergenes,$seq_region_id,$chr,$CHX_for,$CHX_rev,$LTM_for,$LTM_rev);
            
            #Get position checked peaks from LTM reads
            my $TIS = analyze_intergenic_peaks($LTM_reads,$CHX_reads,$ig_id,$intergenes,$min_count,$local_max,$chr);
            
            #Store TIS in TMP csv
            if( keys %{$TIS}){
                store_in_csv($TIS);
            }
        }
        
        #Finish chromosome analysis
        print "     * Finished analyzing chromosome ".$chr."\n";
        
        ### Finish child
        $dbh->disconnect();
        close TMP;
        $pm->finish;  
    }
    
    #Wait all Children
    $pm->wait_all_children();
    
    # Store in db
    print "   Storing all peaks in DB \n";
    store_in_db($dsn_results,$us_results,$pw_results,$id,$work_dir,$chrs,$TMP,$sqlite_db);
}

### Analyze intergenic peaks ###

sub analyze_intergenic_peaks {
    
    #Catch
    my $LTM_reads   =   $_[0];
    my $CHX_reads   =   $_[1];
    my $ig_id       =   $_[2];
    my $intergenes  =   $_[3];
    my $min_count   =   $_[4];
    my $local_max   =   $_[5];
    my $chr         =   $_[6];
    
    #Init
    my $peaks   =   {};
    my $TIS     =   {};
    my $peak    =   "Y";
    my $q       =   "N";
    my $ig      =   $intergenes->{$ig_id};
    my $strand  =   $ig->{'seq_region_strand'};
    my ($ext_cdn,$peak_pos,$start_cdn,$peak_shift,$count,$Rchx,$Rltm,$max,$min,$start,$stop,$i,$rank);
    $local_max = $local_max * 3; # 1 codon = 3 bps.
    
    
    #Only analyze if $LTM_reads not empty
    if(keys %{$LTM_reads}){
        
        #Get peak ext_cdn
        foreach my $key (sort {$a <=> $b} keys %{$LTM_reads}){
            
            $peak           =   "Y";
            $count          =   $LTM_reads->{$key}{'count'};
            $peak_pos       =   $LTM_reads->{$key}{'start'};
            
            # Strand specific
            $ext_cdn = ($strand eq '1') ? get_sequence($chr,$peak_pos-1,$peak_pos+3) : revdnacomp(get_sequence($chr,$peak_pos-3,$peak_pos+1));
            
            # If peak true ATG or near cognate (peak_shift) => TIS
            if(substr($ext_cdn,1,3) =~ m/[ACTG]TG|A[ACTG]G|AT[ACTG]/){
                $start_cdn = substr($ext_cdn,1,3);
                $peak_shift = '0';
            }elsif(substr($ext_cdn,0,3) =~ m/[ACTG]TG|A[ACTG]G|AT[ACTG]/){
                $start_cdn = substr($ext_cdn,0,3);
                $peak_shift = '-1';
                if($strand eq '1'){
                    $peak_pos   =   $peak_pos - 1;
                }elsif($strand eq '-1'){
                    $peak_pos   =   $peak_pos + 1;
                }
                
            }elsif(substr($ext_cdn,2,3) =~ m/[ACTG]TG|A[ACTG]G|AT[ACTG]/){
                $start_cdn = substr($ext_cdn,2,3);
                $peak_shift = '+1';
                if($strand eq '1'){
                    $peak_pos   =   $peak_pos + 1;
                }elsif($strand eq '-1'){
                    $peak_pos   =   $peak_pos - 1;
                }
            }else{
                $start_cdn = $ext_cdn;
                $peak_shift ='0';
                $peak = "N";
            }
            
            #Only save true TIS
            if($peak eq "Y"){
                
                #If already peak on that position, combine
                if(exists $peaks->{$peak_pos}){
                    
                    $peaks->{$peak_pos}->{'peak_shift'}   =         $peaks->{$peak_pos}->{'peak_shift'}.$peak_shift;
                    $peaks->{$peak_pos}->{'count'}        =         $peaks->{$peak_pos}->{'count'}+$count;
                    
                }else{
                    #Add to peaks if not yet a TIS on that position
                    $peaks->{$peak_pos}                   =     $LTM_reads->{$key};
                    $peaks->{$peak_pos}->{'start'}        =     $peak_pos;
                    $peaks->{$peak_pos}->{'peak_shift'}   =     $peak_shift;
                    $peaks->{$peak_pos}->{'start_cdn'}    =     $start_cdn;
                }
            }
            
        }
        
        # Run over all peaks
        foreach my $key (sort {$a <=> $b} keys %{$peaks}){
            
            $count      =       $peaks->{$key}{'count'};
            $peak_pos   =       $peaks->{$key}{'start'};
            $q = "N";
            
            # A TIS has a minimal count of $min_count
            if($count >= $min_count){
                
                # A TIS has to be the local max within $local_max codon(s) up AND downstream of the putative TIS cdna position
                for($i=1;$i<=$local_max;$i++){
                    $max=$key+$i;
                    if(exists $peaks->{$max}){
                        if ($peaks->{$max}{'count'} > $count){
                            $q = "Y";
                            next;
                        }
                    }
                    $min=$key-$i;
                    if(exists $peaks->{$min}){
                        if ($peaks->{$min}{'count'} > $count){
                            $q = "Y";
                            next;
                        }
                    }
                    
                }
                
                # If $q = "Y",  peak not local max within $local_max basepairs
                if($q eq "Y"){
                    next;
                }
                
                # Store Peak in Hash TIS
                $TIS->{$key}                            =   $peaks->{$key};
                $TIS->{$key}->{'intergenic_id'}         =   $ig_id;
                $TIS->{$key}->{'annotation'}            =   "intergenic";
                
                #Get distances to closest up- and downstream genes
                $TIS->{$key}->{'upstream_gene_distance'}        =   ($strand eq '1') ? $intergenes->{$ig_id}{'seq_region_end'} - $key :  $key - $intergenes->{$ig_id}{'seq_region_start'};
                $TIS->{$key}->{'downstream_gene_distance'}      =   ($strand eq '1') ? $key - $intergenes->{$ig_id}{'seq_region_start'} :  $intergenes->{$ig_id}{'seq_region_end'} - $key;
                
            }else{next;}
        }
    }
    
    #Return
    return($TIS);
}

### Store in tmp csv ##

sub store_in_csv{
    
    # Catch
    my $TIS         =   $_[0];
    
    # Append TISs
    foreach my $key (sort keys %{$TIS}){
        
        #if (!defined $TIS->{$key}{'start'}){
        #    print Dumper($TIS->{$key});
        #}
        
        print TMP $TIS->{$key}{'intergenic_id'}.",".$TIS->{$key}{'chr'}.",".$TIS->{$key}{'strand'}.",".$TIS->{$key}{'start'}.",".$TIS->{$key}{'annotation'}.",".$TIS->{$key}->{'upstream_gene_distance'}.",".$TIS->{$key}->{'downstream_gene_distance'}.",".$TIS->{$key}{'start_cdn'}.",".$TIS->{$key}{'peak_shift'}.",".$TIS->{$key}{'count'}."\n";
    }
}

### Store in DB ##

sub store_in_db{
    
    # Catch
    my $dsn             =   $_[0];
    my $us              =   $_[1];
    my $pw              =   $_[2];
    my $id              =   $_[3];
    my $work_dir        =   $_[4];
    my $chrs            =   $_[5];
    my $TMP             =   $_[6];
    my $sqlite_db       =   $_[7];
    
    # Init
    my $dbh     =   dbh($dsn,$us,$pw);
    my $table   =   "TIS_intergenic_".$id;
    
    # Create table 
    my $query = "CREATE TABLE IF NOT EXISTS `".$table."` (
    `intergene_id` varchar(128) NOT NULL default '',
    `chr` char(50) NOT NULL default '',
    `strand` int(2) NOT NULL default '',
    `start` int(10) NOT NULL default '',
    `annotation` varchar(128) NOT NULL default 'NA',
    `upstream_gene_distance` int(10) NOT NULL default '',
    `downstream_gene_distance` int(10) NOT NULL default '',
    `start_codon` varchar(128) NOT NULL default '',
    `peak_shift` int(2) NOT NULL default '',
    `count` float default NULL)"  ;
    $dbh->do($query);
    
    # Add indexes
    my $query_idx = "CREATE INDEX IF NOT EXISTS ".$table."_intergene_idx ON ".$table." (intergene_id)";
    $dbh->do($query_idx);
    $query_idx = "CREATE INDEX IF NOT EXISTS ".$table."_chr_idx ON ".$table." (chr)";
    $dbh->do($query_idx);    
    
    # Store
    foreach my $chr (sort keys %{$chrs}){
        system("sqlite3 -separator , ".$sqlite_db." \".import ".$TMP."/".$table."_".$chr."_tmp.csv ".$table."\"")== 0 or die "system failed: $?";
    }
    
    #Move to galaxy history
    system ("mv ".$sqlite_db." ".$out_sqlite);
    
    #Disconnect dbh
    $dbh->disconnect();
    
    # Unlink tmp csv files
    foreach my $chr (sort keys %{$chrs}){
        unlink $TMP."/".$table."_".$chr."_tmp.csv";
    }
}

### Get intergene overlapping reads ###

sub get_intergene_overlapping_reads{
    
    # Catch
	my $dbh             =   $_[0];
    my $ig_id           =   $_[1];
	my $intergenes      =   $_[2];
    my $seq_region_id   =   $_[3];
    my $chr             =   $_[4];
    my $CHX_for         =   $_[5];
    my $CHX_rev         =   $_[6];
    my $LTM_for         =   $_[7];
    my $LTM_rev         =   $_[8];
    
    #Init
    my %LTM_reads = ();
    my %CHX_reads = ();
    my %LTM_tr_reads = ();
    my %CHX_tr_reads = ();
    
    my $intergene = $intergenes->{$ig_id};
    my $strand = $intergene->{'seq_region_strand'};
    my ($CHX_all,$LTM_all);
    
    # Sense specific
    if ($strand eq '1') {
        $CHX_all = $CHX_for;
        $LTM_all = $LTM_for;
    }
    if ($strand eq '-1') {
        $CHX_all = $CHX_rev;
        $LTM_all = $LTM_rev;
    }
    
    #Get intergene specific reads if exist
    if($intergenes->{$ig_id}{'CHX'}){
        my @keys = @{$intergenes->{$ig_id}{'CHX'}};
        @CHX_reads{@keys} = @{$CHX_all}{@keys};
    }
    
    if($intergenes->{$ig_id}{'LTM'}){
        my @keys = @{$intergenes->{$ig_id}{'LTM'}};
        @LTM_reads{@keys} = @{$LTM_all}{@keys};
    }
    
    #Return
    return(dclone \%LTM_reads,\%CHX_reads);
}

### Match reads to intergenes ##

sub match_reads_to_intergenes{
    
    #Catch
    my $intergenes  =   $_[0];
    my $CHX_for     =   $_[1];
    my $CHX_rev     =   $_[2];
    my $LTM_for     =   $_[3];
    my $LTM_rev     =   $_[4];
    
    #Init
    my @window = ();
    my ($ig_for,$ig_rev,$ig_for_LTM,$ig_rev_LTM);
    
    #Split intergenes in forward and reverse arrays
    foreach my $ig_id (sort { $intergenes->{$a}{'seq_region_start'} <=> $intergenes->{$b}{'seq_region_start'} } keys %{$intergenes}){
        if ($intergenes->{$ig_id}{'seq_region_strand'} eq '1'){
            push (@$ig_for,$ig_id);
            push (@$ig_for_LTM,$ig_id);
        }else{
            push (@$ig_rev,$ig_id);
            push (@$ig_rev_LTM,$ig_id);
        }
    }
    
    # Loop over CHX_forward
    foreach my $key (sort {$a <=> $b} keys %{$CHX_for}){
        
        # Push all ig_ids to @window where ig_start < window_pos
        foreach my $ig_for_id (@$ig_for){
            if($intergenes->{$ig_for_id}{'seq_region_start'} <= $key){
                push(@window,$ig_for_id);
            }else{last;}
        }
        # Get rid of ig_for elements already in @window
        @$ig_for = grep { $intergenes->{$_}{'seq_region_start'} > $key} @$ig_for;
        
        # Get rid of ig_ids in @$window where ig_end < window_pos
        @window = grep { $intergenes->{$_}{'seq_region_end'} >= $key} @window;
        
        # Loop over window and add read position to window_intergenes
        foreach my $window_id (@window){
            push(@{$intergenes->{$window_id}{'CHX'}},$key);
        }
    }
    
    #Empty @window
    @window = ();
    
    # Loop over CHX_reverse
    foreach my $key (sort {$a <=> $b} keys %{$CHX_rev}){
        
        # Push all ig_ids to @window where ig_start < window_pos
        foreach my $ig_rev_id (@$ig_rev){
            if($intergenes->{$ig_rev_id}{'seq_region_start'} <= $key){
                push(@window,$ig_rev_id);
            }else{last;}
        }
        # Get rid of ig_for elements already in @window
        @$ig_rev = grep { $intergenes->{$_}{'seq_region_start'} > $key} @$ig_rev;
        
        # Get rid of ig_ids in @$window where ig_end < window_pos
        @window = grep { $intergenes->{$_}{'seq_region_end'} >= $key} @window;
        
        # Loop over window and add read position to window_intergenes
        foreach my $window_id (@window){
            push(@{$intergenes->{$window_id}{'CHX'}},$key);
        }
    }
    
    #Empty @window
    @window = ();
    
    # Loop over LTM_forward
    foreach my $key (sort {$a <=> $b} keys %{$LTM_for}){
        
        # Push all ig_ids to @window where ig_start < window_pos
        foreach my $ig_for_LTM_id (@$ig_for_LTM){
            if($intergenes->{$ig_for_LTM_id}{'seq_region_start'} <= $key){
                push(@window,$ig_for_LTM_id);
            }else{last;}
        }
        # Get rid of ig_for elements already in @window
        @$ig_for_LTM = grep { $intergenes->{$_}{'seq_region_start'} > $key} @$ig_for_LTM;
        
        # Get rid of ig_ids in @$window where ig_end < window_pos
        @window = grep { $intergenes->{$_}{'seq_region_end'} >= $key} @window;
        
        # Loop over window and add read position to window_intergenes
        foreach my $window_id (@window){
            push(@{$intergenes->{$window_id}{'LTM'}},$key);
        }
    }
    
    #Empty @window
    @window = ();
    
    # Loop over LTM_reverse
    foreach my $key (sort {$a <=> $b} keys %{$LTM_rev}){
        
        # Push all ig_ids to @window where ig_start < window_pos
        foreach my $ig_rev_LTM_id (@$ig_rev_LTM){
            if($intergenes->{$ig_rev_LTM_id}{'seq_region_start'} <= $key){
                push(@window,$ig_rev_LTM_id);
            }else{last;}
        }
        # Get rid of ig_for elements already in @window
        @$ig_rev_LTM = grep { $intergenes->{$_}{'seq_region_start'} > $key} @$ig_rev_LTM;
        
        # Get rid of ig_ids in @$window where ig_end < window_pos
        @window = grep { $intergenes->{$_}{'seq_region_end'} >= $key} @window;
        
        # Loop over window and add read position to window_intergenes
        foreach my $window_id (@window){
            push(@{$intergenes->{$window_id}{'LTM'}},$key);
        }
    }
    
    #Return
    return($intergenes);
}

### Get Intergenes and reads ###

sub get_intergenes_and_reads {
    
    # Catch
    my $seq_region_id   =   $_[0];
    my $chr             =   $_[1];
    my $dbh_ENS         =   $_[2];
    
    # Init
    my $intergenes  = {};
    my $CHX_for     = {};
    my $LTM_for     = {};
    my $CHX_rev     = {};
    my $LTM_rev     = {};
    my $dbh = dbh($dsn_results,$us_results,$pw_results);
    
    # Get intergenes
    my $query = "SELECT intergene_id,seq_region_id,seq_region_strand,seq_region_start,seq_region_end FROM intergene WHERE seq_region_id = $seq_region_id";
	my $sth = $dbh_ENS->prepare($query);
	$sth->execute();
	$intergenes = $sth->fetchall_hashref('intergene_id');
    
    # Get Reads
    $query = "SELECT * FROM count_fastq1 WHERE chr = '$chr' and strand = '1'";
    $sth = $dbh->prepare($query);
	$sth->execute();
	$CHX_for = $sth->fetchall_hashref('start');
    
    $query = "SELECT * FROM count_fastq2 WHERE chr = '$chr' and strand = '1'";
	$sth = $dbh->prepare($query);
	$sth->execute();
	$LTM_for = $sth->fetchall_hashref('start');
    
    $query = "SELECT * FROM count_fastq1 WHERE chr = '$chr' and strand = '-1'";
	$sth = $dbh->prepare($query);
	$sth->execute();
	$CHX_rev = $sth->fetchall_hashref('start');
    
    $query = "SELECT * FROM count_fastq2 WHERE chr = '$chr' and strand = '-1'";
	$sth = $dbh->prepare($query);
	$sth->execute();
	$LTM_rev = $sth->fetchall_hashref('start');
    
    #Disconnect DBH
    $dbh->disconnect();
    
	# Return
	return($intergenes,$CHX_for,$LTM_for,$CHX_rev,$LTM_rev);
}

### Create Intergenes ###

sub create_intergenes{
    
    # Catch
    my $chrom_file      =   $_[0];
    my $dsn_ENS         =   $_[1];
    my $us_ENS          =   $_[2];
    my $pw_ENS          =   $_[3];
    my $chrs            =   $_[4];
    my $cores           =   $_[5];
    my $TMP             =   $_[6];
    my $ens_db          =   $_[7];
    
    # Init
    my $chr_sizes = get_chr_sizes($chrom_file);
    my $table_name = "intergene";
    
    # Check if table already exists
    my $dbh = dbh($dsn_ENS,$us_ENS,$pw_ENS);
    my $query = "select name from sqlite_master where type='table' and name like '$table_name'";
    my $sth = $dbh->prepare($query);
    my $table_exists = 0;
    
    if($sth->execute()){
        while( my $t = $sth->fetchrow_array() ) {
            if( $t =~ /^$table_name$/ ){
                $table_exists = 1;
                print "     ...intergene table already exists. \n";
                return;
            }
        }
    }
    
    # If table does not exist, create it
    if($table_exists == 0){
        
        # Attempt to create table
        my $create_table_tmp = "CREATE TABLE `intergene_tmp` (
        `seq_region_id` int(11) NOT NULL,
        `seq_region_start` int(11) NOT NULL,
        `seq_region_end` int(11) NOT NULL,
        `seq_region_strand` int(11) NOT NULL )";
        
        my $create_table = "CREATE TABLE `intergene` (
        `intergene_id` INTEGER PRIMARY KEY AUTOINCREMENT,
        `seq_region_id` int(11) NOT NULL,
        `seq_region_start` int(11) NOT NULL,
        `seq_region_end` int(11) NOT NULL,
        `seq_region_strand` int(11) NOT NULL )";
        
        $dbh->do($create_table);
        $dbh->do($create_table_tmp);
    }
    
    # Init multi core
    my $pm = new Parallel::ForkManager($cores);
    
    ## Loop over all chromosomes
    foreach my $chr (sort keys %{$chrs}){
        
        ### Start parallel process
        $pm->start and next;
        
        ### DBH per process
        my $dbh = dbh($dsn_ENS,$us_ENS,$pw_ENS);
        
        ### Open File per process
        open TMP, "+>> ".$TMP."/intergenes_".$chr."_tmp.csv" or die $!;
        
        ### seq_region_id for chromosome
        my $seq_region_id = $chrs->{$chr}{'seq_region_id'};
        
        ## Get genes from table
        my $query = "SELECT gene_id,seq_region_id,seq_region_start,seq_region_end,seq_region_strand FROM gene where seq_region_id = '".$seq_region_id."'";
        my $sth = $dbh->prepare($query);
        $sth->execute();
        my $gene = $sth->fetchall_hashref('gene_id');
        
        #Split genes in forward and reverse arrays
        my ($gene_for,$gene_rev);
        
        foreach my $gene_id (sort { $gene->{$a}{'seq_region_start'} <=> $gene->{$b}{'seq_region_start'} } keys %{$gene}){
            if ($gene->{$gene_id}{'seq_region_strand'} eq '1'){
                push (@$gene_for,$gene_id);
            }else{
                push (@$gene_rev,$gene_id);
            }
        }
        
        ## Make intergenes forward strand
        my $strand = 1;
        my $intergene_start = 1;
        my $intergene_end = 1;
        foreach (@$gene_for){
            
            if ($intergene_start >= $gene->{$_}{'seq_region_start'} && $intergene_start > $gene->{$_}{'seq_region_end'}){
                next;
            }
            
            if ($intergene_start >= $gene->{$_}{'seq_region_start'} && $intergene_start <= $gene->{$_}{'seq_region_end'}){
                $intergene_start = $gene->{$_}{'seq_region_end'} + 1;
                next;
            }
            
            $intergene_end = $gene->{$_}{'seq_region_start'} - 1;
            
            #Save in csv
            print TMP $seq_region_id.",".$intergene_start.",".$intergene_end.",".$strand."\n";
            
            $intergene_start = $gene->{$_}{'seq_region_end'} + 1;
            
        }
        
        #Save last part in csv
        $intergene_end = $chr_sizes->{$chr};
        print TMP $seq_region_id.",".$intergene_start.",".$intergene_end.",".$strand."\n";
        
        ## Make intergenes reverse strand
        $strand = -1;
        $intergene_start = 1;
        $intergene_end = 1;
        
        foreach (@$gene_rev){
            
            if ($intergene_start >= $gene->{$_}{'seq_region_start'} && $intergene_start > $gene->{$_}{'seq_region_end'}){
                next;
            }
            
            if ($intergene_start >= $gene->{$_}{'seq_region_start'} && $intergene_start <= $gene->{$_}{'seq_region_end'}){
                $intergene_start = $gene->{$_}{'seq_region_end'} +1;
                next;
            }
            
            $intergene_end = $gene->{$_}{'seq_region_start'} - 1;
            
            #Save in csv
            print TMP $seq_region_id.",".$intergene_start.",".$intergene_end.",".$strand."\n";
            
            $intergene_start = $gene->{$_}{'seq_region_end'} + 1;
            
        }
        #Save last part in csv
        $intergene_end = $chr_sizes->{$chr};
        print TMP $seq_region_id.",".$intergene_start.",".$intergene_end.",".$strand."\n";
        
        ### Finish child
        $dbh->disconnect();
        $pm->finish;
    }
    #Wait all Children
    $pm->wait_all_children();
    
    #Store in tmp intergene db
    foreach my $chr (sort keys %{$chrs}){
        system("sqlite3 -separator , ".$ens_db." \".import ".$TMP."/intergenes_".$chr."_tmp.csv ".$table_name."_tmp\"")== 0 or die "system failed: $?";
    }
    
    #Save in intergene with auto_increment and drop tmp
    $query = "INSERT INTO intergene (seq_region_id,seq_region_start,seq_region_end,seq_region_strand) SELECT * FROM intergene_tmp";
    $dbh->do($query);
    $query = "drop table intergene_tmp";
    $dbh->do($query);
    
    # Add index
    my $query_idx = "CREATE INDEX IF NOT EXISTS ".$table_name."_chr_idx ON ".$table_name." (seq_region_id)";
    $dbh->do($query_idx);
    $query_idx = "CREATE INDEX IF NOT EXISTS ".$table_name."_intergene_idx ON ".$table_name." (intergene_id)";
    $dbh->do($query_idx);
    
    # Unlink tmp csv files
    foreach my $chr (sort keys %{$chrs}){
        unlink $TMP."/intergenes_".$chr."_tmp.csv";
    }
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
    foreach (@chr){
        
        my $query = "SELECT seq_region_id FROM seq_region where coord_system_id = ".$coord_system_id."  and name = '".$_."'";
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

### Create Bin Chromosomes ##

sub create_BIN_chromosomes {
    
    # Catch
    my $BIN_chrom_dir = $_[0];
    my $cores = $_[1];
    my $chrs = $_[2];
    my $work_dir = $_[3];
    
    #Init
    my $TMP = $work_dir."tmp/";
    
    # Create BIN_CHR directory
    system ("mkdir -p ".$BIN_chrom_dir);
    
    # Create binary chrom files
    ## Init multi core
    my $pm = new Parallel::ForkManager($cores);
    #print "   Using ".$cores." core(s)\n   ---------------\n";
    
    foreach my $chr (keys %$chrs){
        
        ## Start parallel process
        $pm->start and next;
        
        open (CHR,"<".$IGENOMES_ROOT."/Sequence/Chromosomes/".$chr.".fa") || die "Cannot open chr fasta input\n";
        open (CHR_BIN,">".$TMP."/".$chr.".fa");
        
        while (<CHR>){
            #Skip first header line
            chomp($_);
            if ($_ =~ m/^>/) { next; }
            print CHR_BIN $_;
        }
        
        close(CHR);
        close(CHR_BIN);
        system ("mv ".$TMP."/".$chr.".fa ".$BIN_chrom_dir."/".$chr.".fa");
        $pm->finish;
    }
    
    # Finish all subprocesses
    $pm->wait_all_children;
    
    
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
    
    # Return input variables
    #return($run_name,$ensemblversion,$species,$mapper,$unique,$adaptorSeq,$readlength,$readtype);
    return($ensemblversion,$species);
    
}

### Get Analysis ID

sub get_id{
    
    # Catch
    my $dsn                 =   $_[0];
    my $us                  =   $_[1];
    my $pw                  =   $_[2];
    my $local_max           =   $_[3];
    my $min_count           =   $_[4];
    
    # Init
    my $dbh = dbh($dsn,$us,$pw);
    
    # Create table 
    my $query = "CREATE TABLE IF NOT EXISTS `TIS_intergenic_overview` (
    `ID` INTEGER primary key,
    `local_max` int(10) NOT NULL default '',
    `min_count` int(10) NOT NULL default '')"  ;
    $dbh->do($query);
    
    # Add parameters to overview table and get ID
    $dbh->do("INSERT INTO `TIS_intergenic_overview`(local_max,min_count)VALUES($local_max,$min_count)");
    my $id= $dbh->func('last_insert_rowid');
    
    # Return
    return($id);
}

### get reverse complement sequence ###

sub revdnacomp {
    my $dna     =   shift;
    my $revcomp =   reverse($dna);
    $revcomp    =~  tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}

### GET CHR SIZES ###

sub get_chr_sizes {
    
    # Catch
    my $chromosome_sizes = $_[0];
    
    # Work
    my %chr_sizes;
    open (Q,"<".$chromosome_sizes) || die "Cannot open chr sizes input\n";
    while (<Q>){
        my @a = split(/\s+/,$_);
        $chr_sizes{$a[0]} = $a[1];
    }
    
    return(\%chr_sizes);
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
