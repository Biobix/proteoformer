#!/usr/bin/perl -w
use strict;
use warnings;
use DBI;
use DBD::SQLite;
use Data::Dumper;
use Parallel::ForkManager;
use Getopt::Long;
use Storable 'dclone';

################
# Perl command #
################

## TO DO:       
##          Add frame_shift calculator if @translation

# ./TIScalling.pl --dir /data/RIBO_runs/RIBO_Ingolia_GerbenM/ --file lane --name mESC_GA --species mouse --ensembl 72 --cores 22 --mapper STAR --tmp /data/RIBO_runs/RIBO_Ingolia_GerbenM/tmp/ --unique N (--intergenic Y --mincount 10 --localmax 1 --R 0.05 --CHX 1 --LTM 2)

# get the command line arguments
my ($run_name,$species,$CHX_lane,$LTM_lane,$work_dir,$min_count,$local_max,$R,$redo,$cores,$tmpfolder,$ensemblversion,$intergenic_TIS,$mapper,$unique,$seqFileName);
GetOptions(
"dir=s"=>\$work_dir,            	# Path to the working directory                                     mandatory argument
"file=s"=>\$seqFileName,            # Base name of the fastq file, e.g. "lane",                         mandatory argument
"name=s"=>\$run_name, 				# Name of the run                                                   mandatory argument
"species=s"=>\$species, 			# Species name, eg mouse/human/fruitfly                             mandatory argument
"ensembl=i"=>\$ensemblversion,      # Ensembl annotation version, eg 70 (2013),                          mandatory argument
"cores=i"=>\$cores, 				# Number of cores to be used                                        optional argument  # Check if possible -> (default = number of chromosomes)
"mapper:s"=>\$mapper,               # The mapper used for alignment (Bowtie,Bowtie2,STAR,TopHat2)       optional argument   (default = Bowtie) => layered approach (first cDNA, second genomic)
"unique=s" =>\$unique,              # Retain the uniquely (and multiple) mapping reads (Y or N),        mandatory  argument
"tmp:s" =>\$tmpfolder,         		# Folder where temporary files are stored,                          optional  argument (default = $TMP env setting)
"CHX:i" =>\$CHX_lane,               # The lane with the CHX riboseq data (numeric: 1 or 2, ...)         optional argument (default 1)
"LTM:i" =>\$LTM_lane,               # The lane with the LTM/HARR riboseq data (numeric: 1 or 2, ...)    optional argument (default 2)
"intergenic:s"=>\$intergenic_TIS,   # Calling intergenic TISses (Y or N)                                optional argument (default = Y)                                                           
"mincount:i" =>\$min_count,         # The minimum count of riboseq profiles mapping to the TIS site     optional argument (default 10)
"localmax:i" =>\$local_max,         # The range wherein the localmax is (e.g.: 1 means +/- one triplet) optional argument (default 1)
"R:f" => \$R                        # The Rltm - Rchx value calculated based on both CHX and LTM data   optional argument (default .05)
);

open(STDERR, ">&STDOUT");

# comment on these
if ($work_dir){
    print "Working directory                                        : $work_dir\n";
} else {
    die "\nDon't forget to pass the working directory using the --dir or -d argument!\n\n";
}
if ($seqFileName){
    print "Base fastq filename                                      : $seqFileName\n";
} else {
    die "\nDon't forget to pass the base name for the FastQ files using the --file or -f argument!\n\n";
}
if ($run_name){
    print "Run name                                                 : $run_name\n";
} else {
    die "\nDon't forget to pass the Run Name using the --name or -n argument!\n\n";
}
if ($species){
    print "Species                                                  : $species\n";
} else {
    die "\nDon't forget to pass the Species name using the --species or -sp argument!\n\n";
}
if ($ensemblversion){
    print "Ensembl Version                                          : $ensemblversion\n";
} else {
    die "\nDon't forget to pass the Ensembl Version using the --ensembl or -ens argument!\n\n";
}
if ($cores){
    print "Number of cores to use for Mapping                       : $cores\n";
} else {
    die "\nDon't forget to pass number of cores to use for mapping using the --cores or -c argument!\n\n";
}
if ($mapper){
    print "The mapper used is                                       : $mapper\n";
} else {
    #Choose default value for mapper
    $mapper = "Bowtie";
}
if ($CHX_lane){
    print "The lane with CHX info is                                : $CHX_lane\n";
} else {
    #Choose default value for mapper
    $CHX_lane = 1;
}
if ($LTM_lane){
    print "The lane with LTM info is                                : $LTM_lane\n";
} else {
    #Choose default value for mapper
    $LTM_lane = 2;
}
if ($min_count){
    print "The minimum coverage of a TIS is                         : $min_count\n";
} else {
    #Choose default value for mapper
    $min_count = 10;
}
if ($local_max){
    print "The regio for local max is set to                        : $local_max\n";
} else {
    #Choose default value for mapper
    $local_max = 1;
}
if ($R){
    print "The R value used is                                      : $R\n";
} else {
    #Choose default value for mapper
    $R = 0.05;
}
if ($unique){
    print "Unique mapped reads                                      : $unique\n";
} else {
    die "\nDon't forget to pass the unique or multiple read retention parameter --unique or -u argument!\n\n";
    
}
if ($intergenic_TIS){
    print "Calling intergenic Tisses                                : $intergenic_TIS\n";
} else {
    #Choose default value for intergenic TIS calling
    $intergenic_TIS = "Y";
}

my $HOME            = $ENV{'HOME'};
my $IGENOMES_ROOT   = ($ENV{'IGENOMES_ROOT'}) ? $ENV{'IGENOMES_ROOT'} : '/data/igenomes';
my $TMP             = ($ENV{'TMP'}) ? $ENV{'TMP'} : $tmpfolder;
print "The following tmpfolder is used                          : $TMP\n";
print "The following igenomes folder is used                    : $IGENOMES_ROOT\n";

#Check if tmpfolder exists, if not create it...
if (!-d "$TMP") {
    system ("mkdir ". $TMP);
}

#ADDED FOR TEST ISSUES
my $clusterPhoenix = "N";

#Set program run_name
$run_name = $run_name."_".$mapper."_".$unique."_".$ensemblversion;

#Conversion for species terminology
my $spec = ($species eq "mouse") ? "Mus_musculus" : ($species eq "human") ? "Homo_sapiens" : "";
my $spec_short = ($species eq "mouse") ? "mmu" : ($species eq "human") ? "hsa" : "";

#Old mouse assembly = NCBIM37, new one is GRCm38
my $assembly = ($species eq "mouse" && $ensemblversion >= 70 ) ? "GRCm38"
: ($species eq "mouse" && $ensemblversion < 70 ) ? "NCBIM37"
: ($species eq "human") ? "GRCh37"
: "";

my $chrom_file = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Annotation/Genes/ChromInfo.txt";
my $BIN_chrom_dir = ($clusterPhoenix eq "Y") ? $HOME."/igenomes_extra/".$spec."/Ensembl/".$assembly."/Sequence/Chromosomes_BIN" : $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Sequence/Chromosomes_BIN";

# DB settings
# Sqlite Riboseq  
my $db_results  = $work_dir."SQLite/".$run_name."_results.db";
my $dsn_results = "DBI:SQLite:dbname=$db_results";
my $us_results  = "";
my $pw_results  = "";

# Sqlite Ensembl
my $db_ENS  = $work_dir."SQLite/ENS_".$spec_short."_".$ensemblversion.".db";
my $dsn_ENS = "DBI:SQLite:dbname=$db_ENS";
my $us_ENS  = "";
my $pw_ENS  = "";


#############################################
## Peak Calling on exon_covered transcripts #
##                                          #
#############################################

# Start time
my $start = time;

print "\nGet analysis id...";
## Get Analysis ID
my $id = get_id($dsn_results,$us_results,$pw_results,$min_count,$local_max,$R,$run_name,$CHX_lane,$LTM_lane,$ensemblversion,$spec,$seqFileName,$intergenic_TIS);
print "                                       : ".$id." \n";

print "Get chromosomes... \n";
## Get chromosomes based on seq_region_id ##
my $chrs = get_chrs($dsn_ENS,$us_ENS,$pw_ENS,$chrom_file,$assembly);

# Create binary chromosomes if they don't exist
print "Checking/Creating binary chrom files ...\n";
if (!-d "$BIN_chrom_dir") {
    create_BIN_chromosomes($BIN_chrom_dir,$cores,$chrs,$work_dir);
}

## Create intergenes  | if intergenic TIS = Yes 
if ($intergenic_TIS eq "Y"){
    print "Creating intergene table... \n";
    create_intergenes($chrom_file,$dsn_ENS,$us_ENS,$pw_ENS,$chrs,$cores,$run_name,$TMP,$ensemblversion,$work_dir,$spec_short);
}

## Tiscalling for all ribo_reads overlapping ENS transcripts
print " \n  Starting transcript peak analysis per chromosome using ".$cores." cores...\n";    
TIS_calling($chrs,$dsn_ENS,$us_ENS,$pw_ENS,$run_name,$CHX_lane,$LTM_lane,$cores,$id,$work_dir,$TMP);

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
    my $chrs        =   $_[0];
    my $db_ENS      =   $_[1];
    my $us_ENS      =   $_[2];
    my $pw_ENS      =   $_[3];
    my $run_name    =   $_[4];
    my $CHX_lane    =   $_[5];
    my $LTM_lane    =   $_[6];
    my $cores       =   $_[7];
    my $id          =   $_[8];
    my $work_dir    =   $_[9];
    my $TMP         =   $_[10];
    
    # Init multi core
    my $pm = new Parallel::ForkManager($cores);
    
    ## Loop over all chromosomes
    foreach my $chr (sort keys %{$chrs}){
        
        ### Start parallel process
        $pm->start and next;
        
        ### DBH per process
        my $dbh = dbh($dsn_ENS,$us_ENS,$pw_ENS);
        
        ### Open File per process
        my $table = $run_name."_TIS_".$id;
        open TMP, "+>> ".$TMP."/".$table."_".$chr."_tmp.csv" or die $!;
    
        # seq_region_id per chromosome
        my $seq_region_id = $chrs->{$chr}{'seq_region_id'};
    
        # Get transcriptsand all reads
        my($trs,$CHX_for,$LTM_for,$CHX_rev,$LTM_rev) = get_transcripts_and_reads($seq_region_id,$chr,$run_name,$CHX_lane,$LTM_lane,$seqFileName);
    
        # Match reads to trancripts
        $trs = match_reads_to_transcripts($trs,$CHX_for,$CHX_rev,$LTM_for,$LTM_rev); 
    
        # Loop over transcripts
        foreach my $tr_id (sort {$a <=> $b} keys %{$trs}){
        
            #Get transcript specific reads and translation data
            my($trans,$exons,$tr_seq,$LTM_reads,$CHX_reads,$trs) = get_transcript_translation_data_and_overlapping_reads($dbh,$tr_id,$trs,$seq_region_id,$chr,$CHX_for,$CHX_rev,$LTM_for,$LTM_rev);  
        
            #Get position checked peaks from LTM reads
            my $TIS = analyze_transcript_peaks($trans,$exons,$tr_seq,$LTM_reads,$CHX_reads,$tr_id,$trs,$min_count,$local_max,$R,$chr);
        
            #Store TIS in TMP csv
            if( keys %{$TIS}){
            store_in_csv($TIS);
            }
        }
        
        ### TIScalling for all ribo profiles overlapping intergenes | if intergenic TIS = Yes
        if ($intergenic_TIS eq "Y"){
            
            # Get intergenes
            my $intergenes = get_intergenes($seq_region_id,$dbh);
            
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
    store_in_db($dsn_results,$us_results,$pw_results,$run_name,$id,$work_dir,$chrs,$TMP);
}

### Analyze peaks ###

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
    
    # Total number of LTM and CHX reads in transcript
    my $total_LTM = 0;
    my $total_CHX = 0;
    
    #Only analyze if $LTM_reads not empty
    if(keys %{$LTM_reads}){
        
        foreach my $key (keys %{$CHX_reads}){
            $total_CHX += $CHX_reads->{$key}{'count'};
        }
        
        #Get peak ext_cdn
        foreach my $key (sort {$a <=> $b} keys %{$LTM_reads}){
            
            $peak           =   "Y";
            $count          =   $LTM_reads->{$key}{'count'};
            $total_LTM      +=  $count;
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
                    $peaks->{$peak_pos}                   =       $LTM_reads->{$key};
                    $peaks->{$peak_pos}->{'start'}        =       $peak_pos;
                    $peaks->{$peak_pos}->{'peak_shift'}   =       $peak_shift;
                    $peaks->{$peak_pos}->{'start_cdn'}    =       $start_cdn;
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
                $TIS->{$key}->{'transcript_id'}         =   $ig_id;
                $TIS->{$key}->{'stable_id'}             =   "NA";
                $TIS->{$key}->{'biotype'}               =   "NA";
                $TIS->{$key}->{'dist_to_trans_start'}   =   "NA";
                $TIS->{$key}->{'dist_to_atis'}          =   "NA";
                $TIS->{$key}->{'annotation'}            =   "intergenic";
                $TIS->{$key}->{'Rltm_min_Rchx'}         =   "NA";
                
            }else{next;}
        }
    }
    
    #Return
    return($TIS);
}

### Analyze transcript peaks ###

sub analyze_transcript_peaks{
    
    #Catch
    my $trans       =   $_[0];
    my $exons       =   $_[1];
    my $tr_seq      =   $_[2];
    my $LTM_reads   =   $_[3];
    my $CHX_reads   =   $_[4];
    my $tr_id       =   $_[5];
    my $trs         =   $_[6];
    my $min_count   =   $_[7];
    my $local_max   =   $_[8];
    my $R           =   $_[9];
    my $chr         =   $_[10];
    
    #Init
    my $peaks   =   {};
    my $TIS     =   {};
    my $peak    =   "Y";
    my $q       =   "N";
    my $tr      =   $trs->{$tr_id};
    my $strand  =   $tr->{'seq_region_strand'};
    my $end_ex_rank     =   $tr->{'end_exon_rank'};
    my $start_ex_rank   =   $tr->{'start_exon_rank'};
    my ($ext_cdn,$cdna_pos,$peak_pos,$start_cdn,$peak_shift,$count,$Rchx,$Rltm,$max,$min,$start,$stop,$i,$rank);
    $local_max = $local_max * 3; # 1 codon = 3 bps.
    
    # Total number of LTM and CHX reads in transcript
    my $total_LTM = 0;
    my $total_CHX = 0;
    
    #Only analyze if $LTM_reads not empty
    if(keys %{$LTM_reads}){
        
        foreach my $key (keys %{$CHX_reads}){
            $total_CHX += $CHX_reads->{$key}{'count'};
        }
        
        #Get peak cDNA positions and ext_cdn
        foreach my $key (sort {$a <=> $b} keys %{$LTM_reads}){
            
            $peak           =   "Y";
            $count          =   $LTM_reads->{$key}{'count'};
            $total_LTM      +=  $count;
            $peak_pos       =   $LTM_reads->{$key}{'start'};
            $rank           =   $LTM_reads->{$key}{'exon_rank'};
            $cdna_pos       =   cdna_position($peak_pos,$exons,$rank);
            
            # Strand specific
            # in a string the first position is called 0, however, it's cDNA position = 1. ==> $cdna_pos - 1 = correct string position
            # substr does not work for negative values, so first get rid of small cdna_pos
            
            if($cdna_pos < 4){
                $ext_cdn = ($strand eq '1') ? get_sequence($chr,$peak_pos-1,$peak_pos+3) : revdnacomp(get_sequence($chr,$peak_pos-3,$peak_pos+1));
            }else{
                $ext_cdn = ($strand eq '1') ? substr($tr_seq,$cdna_pos -2,5) : revdnacomp(substr($tr_seq,$cdna_pos-4,5));
            }
            
            # If peak at end of transcript sequence, $ext_cdn will be < 5, get genomic sequence
            if(length($ext_cdn)<5){
                $ext_cdn = ($strand eq '1') ? get_sequence($chr,$peak_pos-1,$peak_pos+3) : revdnacomp(get_sequence($chr,$peak_pos-3,$peak_pos+1));
            } 
            
            # If peak true ATG or near cognate (peak_shift) => TIS
            if(substr($ext_cdn,1,3) =~ m/[ACTG]TG|A[ACTG]G|AT[ACTG]/){
                $start_cdn = substr($ext_cdn,1,3);
                $peak_shift = '0';
            }elsif(substr($ext_cdn,0,3) =~ m/[ACTG]TG|A[ACTG]G|AT[ACTG]/){
                $start_cdn = substr($ext_cdn,0,3);
                $peak_shift = '-1';
                if($strand eq '1'){
                    $peak_pos   =   $peak_pos - 1;
                    $cdna_pos   =   $cdna_pos - 1;
                }elsif($strand eq '-1'){
                    $peak_pos   =   $peak_pos + 1;
                    $cdna_pos   =   $cdna_pos + 1;
                }
                
            }elsif(substr($ext_cdn,2,3) =~ m/[ACTG]TG|A[ACTG]G|AT[ACTG]/){
                $start_cdn = substr($ext_cdn,2,3);
                $peak_shift = '+1';
                if($strand eq '1'){
                    $peak_pos   =   $peak_pos + 1;
                    $cdna_pos   =   $cdna_pos + 1;
                }elsif($strand eq '-1'){
                    $peak_pos   =   $peak_pos - 1;
                    $cdna_pos   =   $cdna_pos - 1;
                }
            }else{
                $start_cdn = $ext_cdn;
                $peak_shift ='0';
                $peak = "N";
            }
            
            #Only save true TIS
            if($peak eq "Y"){
                
                #If already peak on that position, combine
                if(exists $peaks->{$cdna_pos}){
                    
                    $peaks->{$cdna_pos}->{'peak_shift'}   =         $peaks->{$cdna_pos}->{'peak_shift'}.$peak_shift;
                    $peaks->{$cdna_pos}->{'count'}        =         $peaks->{$cdna_pos}->{'count'}+$count;
                    
                }else{
                    #Add to peaks if not yet a TIS on that position!! Whatch out for non-ATG or non-near-cognate peaks!!!
                    $peaks->{$cdna_pos}                   =       $LTM_reads->{$key};
                    $peaks->{$cdna_pos}->{'start'}        =       $peak_pos;
                    $peaks->{$cdna_pos}->{'peak_shift'}   =       $peak_shift;
                    $peaks->{$cdna_pos}->{'start_cdn'}    =       $start_cdn;
                    $peaks->{$cdna_pos}->{'cdna_pos'}     =       $cdna_pos;
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
                
                #Calculate RLTM - RCHX
                #Rk = (Xk/Nk) x 10 (k = LTM, CHX), Xk number of reads on that position in data k, Nk total number of reads for transcript.
                
                # Calculate Rltm
                $Rltm = ($count / $total_LTM) * 10;
                if(exists $CHX_reads->{$peak_pos}{'count'}){
                    $Rchx = ($CHX_reads->{$peak_pos}{'count'}/$total_CHX);
                    $Rchx = $Rchx*10;
                }else{
                    $Rchx=0;
                }
                my $R_test = $Rltm - $Rchx;
                
                # A TIS must have an RLTM - RCHX value > $R
                if(($Rltm - $Rchx) < $R){next;}
                
                # Store Peak in Hash TIS
                $TIS->{$key}                        =   $peaks->{$key};
                $TIS->{$key}->{'Rltm_min_Rchx'}     =   $Rltm-$Rchx;
                $TIS->{$key}->{'transcript_id'}     =   $tr_id;
                
            }else{next;}
        }
        
        #Run over all TIS
        
        # If empty translation ==> biotype without translation info
        if(@{$trans}){
            
            #Get aTIS from first exon
            $start = ($strand eq '1') ? $exons->{$start_ex_rank}{'e_start'} + ${$trans}[2] -1 : $exons->{$start_ex_rank}{'e_end'} - ${$trans}[2] + 1;
            
            #Get stop codon position
            $stop =  ($strand eq '1') ?  $exons->{$end_ex_rank}{'e_start'} + ${$trans}[3] -1 : $exons->{$end_ex_rank}{'e_end'} - ${$trans}[3] +1;
            
            # Compare start peak position with aTis position
            # Compare start peak position with transcript start position
            # Annotate peak position (5'UTR,aTIS,CDS,3'UTR)
            foreach my $key (sort {$a <=> $b} keys %{$TIS}){
                
                $TIS->{$key}->{'biotype'}   =   $tr->{'biotype'};
                $TIS->{$key}->{'stable_id'} =   $tr->{'stable_id'};
                $cdna_pos                   =   $TIS->{$key}{'cdna_pos'};
                
                if($strand eq '1'){
                    
                    $TIS->{$key}->{'dist_to_atis'}          =       $cdna_pos - $start;
                    $TIS->{$key}->{'dist_to_trans_start'}   =       $cdna_pos;
                    
                    if($cdna_pos < $start){
                        $TIS->{$key}->{'annotation'}        =       "5'UTR";
                    }elsif($cdna_pos == $start){
                        $TIS->{$key}->{'annotation'}        =       "aTIS";
                    }elsif($cdna_pos > $stop){
                        $TIS->{$key}->{'annotation'}        =       "3'UTR";
                    }else{
                        $TIS->{$key}->{'annotation'}        =       "CDS";
                    }
                }elsif($strand eq '-1'){
                    
                    $TIS->{$key}->{'dist_to_atis'}          =       $start - $cdna_pos;
                    $TIS->{$key}->{'dist_to_trans_start'}   =       length($tr_seq) - $cdna_pos;
                    
                    if($cdna_pos > $start){
                        $TIS->{$key}->{'annotation'}        =       "5'UTR";
                    }elsif($cdna_pos == $start){
                        $TIS->{$key}->{'annotation'}        =       "aTIS";
                    }elsif($cdna_pos < $stop){
                        $TIS->{$key}->{'annotation'}        =       "3'UTR";
                    }else{
                        $TIS->{$key}{'annotation'} = "CDS";
                    }
                }
            }    
        }else{
            
            # Compare start peak position with transcript start position
            foreach my $key (sort {$a <=> $b} keys %{$TIS}){
                
                $TIS->{$key}{'annotation'}          =   "no_translation";
                $TIS->{$key}{'dist_to_atis'}        =   "NA";
                $TIS->{$key}{'biotype'}             =   $tr->{'biotype'};
                $TIS->{$key}->{'stable_id'}         =   $tr->{'stable_id'};
                $cdna_pos                           =   $TIS->{$key}{'cdna_pos'};
                $TIS->{$key}{'dist_to_trans_start'} =   ($strand eq '1') ? $cdna_pos : length($tr_seq) - $cdna_pos;
                
            }
        }
    }
    
    #Return
    return($TIS);
}

### Get transcript specific translation data ###

sub get_transcript_translation_data_and_overlapping_reads{
    
    # Catch
	my $dbh             =   $_[0];
    my $tr_id           =   $_[1];
	my $trs             =   $_[2];
    my $seq_region_id   =   $_[3];
    my $chr             =   $_[4];
    my $CHX_for         =   $_[5];
    my $CHX_rev         =   $_[6];
    my $LTM_for         =   $_[7];
    my $LTM_rev         =   $_[8];
    
    #Init
    my $exons = {};
    my $trans = [];
    my %LTM_reads = ();
    my %CHX_reads = ();
    my %LTM_tr_reads = ();
    my %CHX_tr_reads = ();
    
    my $e_seq = "";
    my $tr_seq = "";
    my $tr = $trs->{$tr_id};
    my $strand = $tr->{'seq_region_strand'};
    my ($e_start,$e_end,$e_id,$tr_start,$tr_end,$CHX_all,$LTM_all);
    
    # Sense specific
    if ($strand eq '1') { 
        $CHX_all = $CHX_for;
        $LTM_all = $LTM_for;
    }
    if ($strand eq '-1') {
        $CHX_all = $CHX_rev;
        $LTM_all = $LTM_rev;
    }
    
    # Get translation data for transcrip_id if available
    my $query = "SELECT start_exon_id,end_exon_id,seq_start,seq_end FROM translation where transcript_id = '".$tr_id."'";
    my $sth = $dbh->prepare($query);
	$sth->execute();
    $trans = $sth->fetchrow_arrayref();
    

    # Get exons from ENS DB
    $query = "SELECT a.rank,a.exon_id,b.seq_region_start,b.seq_region_end,b.phase,b.end_phase FROM exon_transcript a join exon b on a.exon_id = b.exon_id where a.transcript_id = '".$tr_id."' AND b.seq_region_id = '".$seq_region_id."'";
    $sth = $dbh->prepare($query);
    $sth->execute();
    $exons = $sth->fetchall_hashref('rank');
        
    #Get transcript specific reads if exist 
    if($trs->{$tr_id}{'CHX'}){
        my @keys = @{$trs->{$tr_id}{'CHX'}};
        @CHX_tr_reads{@keys} = @{$CHX_all}{@keys};
    }
    
    if($trs->{$tr_id}{'LTM'}){
        my @keys = @{$trs->{$tr_id}{'LTM'}};
        @LTM_tr_reads{@keys} = @{$LTM_all}{@keys};
    }
    
    #Get transcript sequence based on exons
    #Construct sequence by concatenating exon sequences, always forward strand
    if($strand eq '1'){
        
        # Concatenate sequence based on exons
        foreach my $exon (sort {$a <=> $b} keys %{$exons}){
            
            #Get exon from exons hash
            my $ex      =   $exons->{$exon};
            
            $e_start    =   $ex->{'seq_region_start'};
            $e_end      =   $ex->{'seq_region_end'};
            $e_id       =   $ex->{'exon_id'};
            
            $e_seq = get_sequence($chr,$e_start,$e_end);
            $exons->{$exon}->{'e_start'} = length($tr_seq) + 1;
            $tr_seq = $tr_seq.$e_seq;
            $exons->{$exon}->{'e_end'} = length($tr_seq);
                
            if($trans->[0]){
                
                # Add rank of start and end exon to $trans                                          
                if($ex->{'exon_id'} == ${$trans}[0]){                                     
                    $trs->{$tr_id}{'start_exon_rank'} = $ex->{'rank'};                            
                }                                                                                      
                if($ex->{'exon_id'} == ${$trans}[1]){                                         
                    $trs->{$tr_id}{'end_exon_rank'} = $ex->{'rank'};                                     
                }                
            }     

            #Grep exon specific reads
            my @keys = grep { $CHX_tr_reads{$_}{'start'} > $e_start -1} keys %CHX_tr_reads;
            @keys = grep { $_ < $e_end +1 } @keys;
            @CHX_reads{@keys} = @CHX_tr_reads{@keys};
            
                
            @keys = grep { $LTM_tr_reads{$_}{'start'} > $e_start -1} keys %LTM_tr_reads;
            @keys = grep { $_ < $e_end +1 } @keys;
            @LTM_reads{@keys} = @LTM_tr_reads{@keys};
            
            foreach my $key (@keys){
                $LTM_reads{$key}{'exon_rank'} = $ex->{'rank'};
            }
        }
    }elsif($strand eq '-1'){
            
        # Concatenate sequence based on exons
        foreach my $exon (sort {$b <=> $a} keys %{$exons}){
            
            #Get exon from exons hash
            my $ex      =   $exons->{$exon};
            
            $e_start    =   $ex->{'seq_region_start'};
            $e_end      =   $ex->{'seq_region_end'};
            $e_id       =   $ex->{'exon_id'};
            
            $e_seq = get_sequence($chr,$e_start,$e_end);
            $exons->{$exon}->{'e_start'} = length($tr_seq) + 1;
            $tr_seq = $tr_seq.$e_seq;
            $exons->{$exon}->{'e_end'} = length($tr_seq);
            
            if($trans->[0]){   
                 
                # Add rank of start and end exon to $trans                                         
                if($ex->{'exon_id'} == ${$trans}[0]){                                     
                    $trs->{$tr_id}->{'start_exon_rank'} = $ex->{'rank'};                           
                }                                                                                        
                if($ex->{'exon_id'} == ${$trans}[1]){                                        
                    $trs->{$tr_id}->{'end_exon_rank'} = $ex->{'rank'};                                     
                }                
            }   
                
            #Grep exon specific reads
            my @keys = grep { $CHX_tr_reads{$_}{'start'} > $e_start -1} keys %CHX_tr_reads;
            @keys = grep { $_ < $e_end +1 } @keys;
            @CHX_reads{@keys} = @CHX_tr_reads{@keys};
            
            @keys = grep { $LTM_tr_reads{$_}{'start'} > $e_start -1} keys %LTM_tr_reads;
            @keys = grep { $_ < $e_end +1 } @keys;
            @LTM_reads{@keys} = @LTM_tr_reads{@keys};
            
            foreach my $key (@keys){
                $LTM_reads{$key}{'exon_rank'} = $ex->{'rank'};
            }
        }
    }

    #Return
    return($trans,$exons,$tr_seq,dclone(\%LTM_reads),\%CHX_reads,$trs);
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

### Match reads to transcripts ##

sub match_reads_to_transcripts{
    
    #Catch
    my $trs         =   $_[0];
    my $CHX_for     =   $_[1];
    my $CHX_rev     =   $_[2];
    my $LTM_for     =   $_[3];
    my $LTM_rev     =   $_[4];
    
    #Init
    my @window = ();
    my ($tr_for,$tr_rev,$tr_for_LTM,$tr_rev_LTM);         
    
    #Split transcripts in forward and reverse arrays
    foreach my $tr_id (sort { $trs->{$a}{'seq_region_start'} <=> $trs->{$b}{'seq_region_start'} } keys %{$trs}){
        if ($trs->{$tr_id}{'seq_region_strand'} eq '1'){
            push (@$tr_for,$tr_id);
            push (@$tr_for_LTM,$tr_id);
        }else{
            push (@$tr_rev,$tr_id);
            push (@$tr_rev_LTM,$tr_id);
        }
    }
    
    # Loop over CHX_forward
    foreach my $key (sort {$a <=> $b} keys %{$CHX_for}){
        
        # Push all tr_ids to @window where tr_start < window_pos
        foreach my $tr_for_id (@$tr_for){
            if($trs->{$tr_for_id}{'seq_region_start'} <= $key){
                push(@window,$tr_for_id);
            }else{last;}
        }
        # Get rid of tr_for elements already in @window
        @$tr_for = grep { $trs->{$_}{'seq_region_start'} > $key} @$tr_for;
        
        # Get rid of tr_ids in @$window where tr_end < window_pos
        @window = grep { $trs->{$_}{'seq_region_end'} >= $key} @window;
        
        # Loop over window and add read position to window_transcripts
        foreach my $window_id (@window){
            push(@{$trs->{$window_id}{'CHX'}},$key);
        }
    }
    
    #Empty @window
    @window = ();
    
    # Loop over CHX_reverse
    foreach my $key (sort {$a <=> $b} keys %{$CHX_rev}){
        
        # Push all tr_ids to @window where tr_start < window_pos
        foreach my $tr_rev_id (@$tr_rev){
            if($trs->{$tr_rev_id}{'seq_region_start'} <= $key){
                push(@window,$tr_rev_id);
            }else{last;}
        }
        # Get rid of tr_for elements already in @window
        @$tr_rev = grep { $trs->{$_}{'seq_region_start'} > $key} @$tr_rev;
        
        # Get rid of tr_ids in @$window where tr_end < window_pos
        @window = grep { $trs->{$_}{'seq_region_end'} >= $key} @window;
        
        # Loop over window and add read position to window_transcripts
        foreach my $window_id (@window){
            push(@{$trs->{$window_id}{'CHX'}},$key);
        }
    }
    
    #Empty @window
    @window = ();
    
    # Loop over LTM_forward
    foreach my $key (sort {$a <=> $b} keys %{$LTM_for}){
        
        # Push all tr_ids to @window where tr_start < window_pos
        foreach my $tr_for_LTM_id (@$tr_for_LTM){
            if($trs->{$tr_for_LTM_id}{'seq_region_start'} <= $key){
                push(@window,$tr_for_LTM_id);
            }else{last;}
        }
        # Get rid of tr_for elements already in @window
        @$tr_for_LTM = grep { $trs->{$_}{'seq_region_start'} > $key} @$tr_for_LTM;
        
        # Get rid of tr_ids in @$window where tr_end < window_pos
        @window = grep { $trs->{$_}{'seq_region_end'} >= $key} @window;
        
        # Loop over window and add read position to window_transcripts
        foreach my $window_id (@window){
            push(@{$trs->{$window_id}{'LTM'}},$key);
        }
    }
    
    #Empty @window
    @window = ();
    
    # Loop over LTM_reverse
    foreach my $key (sort {$a <=> $b} keys %{$LTM_rev}){
        
        # Push all tr_ids to @window where tr_start < window_pos
        foreach my $tr_rev_LTM_id (@$tr_rev_LTM){
            if($trs->{$tr_rev_LTM_id}{'seq_region_start'} <= $key){
                push(@window,$tr_rev_LTM_id);
            }else{last;}
        }
        # Get rid of tr_for elements already in @window
        @$tr_rev_LTM = grep { $trs->{$_}{'seq_region_start'} > $key} @$tr_rev_LTM;
        
        # Get rid of tr_ids in @$window where tr_end < window_pos
        @window = grep { $trs->{$_}{'seq_region_end'} >= $key} @window;
        
        # Loop over window and add read position to window_transcripts
        foreach my $window_id (@window){
            push(@{$trs->{$window_id}{'LTM'}},$key);
        }
    }
    
    #Return
    return($trs);
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

### Get Transcripts from transcript calling ###

sub get_transcripts_and_reads {
    
    # Catch
    my $seq_region_id   =   $_[0];
    my $chr             =   $_[1];
    my $run_name        =   $_[2];
    my $CHX_lane        =   $_[3];
    my $LTM_lane        =   $_[4];
    my $seqFileName     =   $_[5];
    
    # Init
    my $dbh = dbh($dsn_results,$us_results,$pw_results);
    my $trs         = {};
    my $CHX_for     = {};
    my $LTM_for     = {};
    my $CHX_rev     = {};
    my $LTM_rev     = {};
    
    # Get transcripts
    my $query = "SELECT transcript_id,seq_region_id,seq_region_strand,seq_region_start,seq_region_end,read_counts,biotype,stable_id FROM ".$run_name."_".$seqFileName.$CHX_lane."_transcript_translation WHERE seq_region_id = $seq_region_id AND exon_coverage = 'Yes'";
   	my $sth = $dbh->prepare($query);
	$sth->execute();
	$trs = $sth->fetchall_hashref('transcript_id');
    
    # Get Reads
    $query = "SELECT * FROM bins_".$run_name."_".$seqFileName.$CHX_lane." WHERE chr = '$chr' and strand = '1'";
    $sth = $dbh->prepare($query);
	$sth->execute();
	$CHX_for = $sth->fetchall_hashref('start');
    
    $query = "SELECT * FROM bins_".$run_name."_".$seqFileName.$LTM_lane." WHERE chr = '$chr' and strand = '1'";
	$sth = $dbh->prepare($query);
	$sth->execute();
	$LTM_for = $sth->fetchall_hashref('start');
    
    $query = "SELECT * FROM bins_".$run_name."_".$seqFileName.$CHX_lane." WHERE chr = '$chr' and strand = '-1'";
	$sth = $dbh->prepare($query);
	$sth->execute();
	$CHX_rev = $sth->fetchall_hashref('start');
    
    $query = "SELECT * FROM bins_".$run_name."_".$seqFileName.$LTM_lane." WHERE chr = '$chr' and strand = '-1'";
	$sth = $dbh->prepare($query);
	$sth->execute();
	$LTM_rev = $sth->fetchall_hashref('start');
    
    #Disconnect DBH
    $dbh->disconnect();
    
	# Return
	return($trs,$CHX_for,$LTM_for,$CHX_rev,$LTM_rev);    
}

### Get Intergenes and reads ###

sub get_intergenes {
    
    # Catch
    my $seq_region_id   =   $_[0];
    my $dbh_ENS         =   $_[1];
    
    # Init
    my $intergenes  = {};
    
    # Get intergenes
    my $query = "SELECT intergene_id,seq_region_id,seq_region_strand,seq_region_start,seq_region_end FROM intergene WHERE seq_region_id = $seq_region_id ";
	my $sth = $dbh_ENS->prepare($query);
	$sth->execute();
	$intergenes = $sth->fetchall_hashref('intergene_id');
    
    #Disconnect DBH
    $dbh_ENS->disconnect();
    
	# Return
	return($intergenes);    
}

### Store in DB ##

sub store_in_db{
    
    # Catch
    my $dsn         =   $_[0];
    my $us          =   $_[1];
    my $pw          =   $_[2];
    my $run_name    =   $_[3];
    my $id          =   $_[4];
    my $work_dir    =   $_[5];
    my $chrs        =   $_[6];
    my $TMP         =   $_[7];
    
    # Init
    my $dbh     =   dbh($dsn,$us,$pw);
    my $table   =   $run_name."_TIS_".$id;
    
    # Create table 
    my $query = "CREATE TABLE IF NOT EXISTS `".$table."` (
    `transcript_id` varchar(128) NOT NULL default '',
    `stable_id` varchar(128) NOT NULL default '',
    `biotype` varchar(128) NOT NULL default '',
    `chr` char(50) NOT NULL default '',
    `strand` int(2) NOT NULL default '',
    `start` int(10) NOT NULL default '',
    `dist_to_transcript_start` int(10) NOT NULL default '',
    `dist_to_aTIS` int(10) NOT NULL default 'NA',
    `annotation` varchar(128) NOT NULL default 'NA',
    `start_codon` varchar(128) NOT NULL default '',
    `peak_shift` int(2) NOT NULL default '',
    `count` float default NULL,
    `Rltm_min_Rchx` decimal(11,8) NOT NULL default '0')"  ;
    $dbh->do($query);
    
    # Add indexes
    my $query_idx = "CREATE INDEX IF NOT EXISTS ".$table."_transcript_idx ON ".$table." (transcript_id)";
    $dbh->do($query_idx);
    $query_idx = "CREATE INDEX IF NOT EXISTS ".$table."_chr_idx ON ".$table." (chr)";
    $dbh->do($query_idx);    
    
    # Store
    foreach my $chr (sort keys %{$chrs}){
        system("sqlite3 -separator , ".$work_dir."SQLite/".$run_name."_results.db \".import ".$TMP."/".$table."_".$chr."_tmp.csv ".$table."\"")== 0 or die "system failed: $?";
    }
    #Disconnect dbh
    $dbh->disconnect();
    
    # Unlink tmp csv files
    foreach my $chr (sort keys %{$chrs}){
        unlink $TMP."/".$table."_".$chr."_tmp.csv";
    }
}

### Store in tmp csv ##

sub store_in_csv{
    
    # Catch
    my $TIS         =   $_[0];
    
    # Append TISs
    foreach my $key (sort keys %{$TIS}){
        
        print TMP $TIS->{$key}{'transcript_id'}.",".$TIS->{$key}{'stable_id'}.",".$TIS->{$key}{'biotype'}.",".$TIS->{$key}{'chr'}.",".$TIS->{$key}{'strand'}.",".$TIS->{$key}{'start'}.",".$TIS->{$key}{'dist_to_trans_start'}.",".$TIS->{$key}{'dist_to_atis'}.",".$TIS->{$key}{'annotation'}.",".$TIS->{$key}{'start_cdn'}.",".$TIS->{$key}{'peak_shift'}.",".$TIS->{$key}{'count'}.",".$TIS->{$key}{'Rltm_min_Rchx'}."\n"; 
    }
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
    my $run_name        =   $_[6];
    my $TMP             =   $_[7];
    my $ensemblversion  =   $_[8];
    my $work_dir        =   $_[9];
    my $spec_short      =   $_[10];
    
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
        system("sqlite3 -separator , ".$work_dir."SQLite/ENS_".$spec_short."_".$ensemblversion.".db \".import ".$TMP."/intergenes_".$chr."_tmp.csv ".$table_name."_tmp\"")== 0 or die "system failed: $?";
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
        
        open (CHR,"<".$IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Sequence/Chromosomes/".$chr.".fa") || die "Cannot open chr fasta input\n";
        open (CHR_BIN, ">".$TMP."/".$chr.".fa");
        
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

### Get Analysis ID

sub get_id{
    
    # Catch
    my $dsn                 =   $_[0];
    my $us                  =   $_[1];
    my $pw                  =   $_[2];
    my $min_count           =   $_[3];
    my $local_max           =   $_[4];
    my $R                   =   $_[5];
    my $run_name            =   $_[6];
    my $CHX                 =   $_[7];
    my $LTM                 =   $_[8];
    my $ensemblversion      =   $_[9];
    my $spec                =   $_[10];
    my $seqFileName         =   $_[11];
    my $intergenic_TIS      =   $_[12];
    
    # Init
    my $dbh = dbh($dsn,$us,$pw);
    
    # Create table 
    my $query = "CREATE TABLE IF NOT EXISTS `".$run_name."_TIS_overview` (
    `ID` INTEGER primary key,
    `CHX_lane` varchar(32) NOT NULL default '',
    `LTM_lane` varchar(32) NOT NULL default '',
    `ensembl_version` varchar(32) NOT NULL default '',
    `intergenic` varchar(32) NOT NULL default '',
    `species` varchar(32) NOT NULL default '',
    `min_count` int(10) NOT NULL default '',
    `local_max` int(10) NOT NULL default '',
    `r` decimal(11,8) NOT NULL default '')"  ;
    $dbh->do($query);
    
    # Add parameters to overview table and get ID
    $dbh->do("INSERT INTO `".$run_name."_TIS_overview`(CHX_lane,LTM_lane,ensembl_version,intergenic,species,min_count,local_max,r) VALUES('$seqFileName$CHX','$seqFileName$LTM','$ensemblversion','$intergenic_TIS','$spec',$min_count,$local_max,$R)");
    my $id= $dbh->func('last_insert_rowid');
    
    # Return
    return($id);
}

### get peak cdna position ###

sub cdna_position {
    
    #Catch
    my $pos     =   $_[0];
    my $exons   =   $_[1];
    my $rank    =   $_[2];
    
    #Init
    my $cdna_pos;
    
    #cDNA position always as forward strand position 
    $cdna_pos = $pos - $exons->{$rank}->{'seq_region_start'} + $exons->{$rank}->{'e_start'};
    
    #Return
    return($cdna_pos);
    
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