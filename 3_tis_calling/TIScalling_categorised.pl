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

# ./TIScalling_categorised.pl --sqlite SQLite/results.db --cores 22 --local_max 1 --R_aTIS 0.01 --min_count_aTIS 5 --R_5 0.05 --min_count_5 10 --R_CDS 0.15 --min_count_CDS 15 --R_3 0.05 --min_count_3 10 --R_ntr 0.05 --min_count_ntr 10

## Galaxy ##

# ./TIScalling_categorised.pl --sqlite_db "${sqlite_db}" --ens_db "${ensembl_db.fields.path}" --cores "${cores}" --igenomes_root "${igenomes_root.fields.path}" --out_sqlite "${out_sqlite}"

# get the command line arguments
my ($work_dir,$local_max,$out_sqlite,$tmpfolder,$min_count_aTIS,$R_aTIS,$min_count_5,$R_5,$min_count_CDS,$R_CDS,$min_count_3,$R_3,$min_count_ntr,$R_ntr,$sqlite_db,$transcriptfilter);
GetOptions(
"dir=s"=>\$work_dir,                            # Path to the working directory                                                                 optional argument
"tmp:s" =>\$tmpfolder,                          # Folder where temporary files are stored,                                                      optional  argument (default = $TMP env setting)
"sqlite_db:s" =>\$sqlite_db,                    # SQLite results db with mapping and tr_translation tables                                      mandatory argument
"local_max:i" =>\$local_max,                    # The range wherein the localmax is (e.g.: 1 means +/- one triplet)                             optional argument (default 1)
"min_count_aTIS:i" =>\$min_count_aTIS,          # The minimum count of riboseq profiles mapping to the aTIS site                                optional argument (default 5)
"R_aTIS:f" => \$R_aTIS,                         # The Rltm - Rchx value calculated based on both CHX and LTM data for a aTIS                    optional argument (default .01)
"min_count_5:i" =>\$min_count_5,                # The minimum count of riboseq profiles mapping to the 5'UTR site                               optional argument (default 10)
"R_5:f" => \$R_5,                               # The Rltm - Rchx value calculated based on both CHX and LTM data for a 5'UTR TIS               optional argument (default .05)
"min_count_CDS:i" =>\$min_count_CDS,            # The minimum count of riboseq profiles mapping to the CDS site                                 optional argument (default 15)
"R_CDS:f" => \$R_CDS,                           # The Rltm - Rchx value calculated based on both CHX and LTM data for a CDS TIS                 optional argument (default .15)
"min_count_3:i" =>\$min_count_3,                # The minimum count of riboseq profiles mapping to the 3'UTR site                               optional argument (default 10)
"R_3:f" => \$R_3,                               # The Rltm - Rchx value calculated based on both CHX and LTM data for a 3'UTR TIS               optional argument (default .05)
"min_count_ntr:i" =>\$min_count_ntr,  # The minimum count of riboseq profiles mapping to the no translation site                                optional argument (default 10)
"R_ntr:f" => \$R_ntr,                 # The Rltm - Rchx value calculated based on both CHX and LTM data for a no translation TIS                optional argument (default .05)
"out_sqlite:s" =>\$out_sqlite,                  # Galaxy specific history file location                                                         Galaxy specific
"transcriptfilter:s" =>\$transcriptfilter      # Use certain filters at transcript level to reduce complexity (CCDS-id, canonical)             optional argument (default none), this means that protein_coding and non protein_coding with exon_coverage = 'YES' will be exported
);

my $CWD             = getcwd;
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
if ($local_max){
    print "The regio for local max is set to                        : $local_max\n";
} else {
    #Choose default value for local_max
    $local_max = 1;
}
if ($min_count_aTIS){
    print "The minimum coverage of an aTIS is                       : $min_count_aTIS\n";
} else {
    #Choose default value for min_count
    $min_count_aTIS = 5;
}
if ($min_count_5){
    print "The minimum coverage of a 5'UTR TIS is                   : $min_count_5\n";
} else {
    #Choose default value for min_count
    $min_count_5 = 10;
}
if ($min_count_CDS){
    print "The minimum coverage of a CDS TIS is                     : $min_count_CDS\n";
} else {
    #Choose default value for min_count
    $min_count_CDS = 15;
}
if ($min_count_3){
    print "The minimum coverage of a 3'UTR TIS is                   : $min_count_3\n";
} else {
    #Choose default value for min_count
    $min_count_3 = 10;
}
if ($min_count_ntr){
    print "The minimum coverage of a ntrlation TIS is               : $min_count_ntr\n";
} else {
    #Choose default value for min_count
    $min_count_ntr = 10;
}
if ($R_aTIS){
    print "The R value used for an aTIS is                          : $R_aTIS\n";
} else {
    #Choose default value for R
    $R_aTIS = 0.01;
}
if ($R_5){
    print "The R value used for a 5'UTR TIS is                      : $R_5\n";
} else {
    #Choose default value for R
    $R_5 = 0.05;
}
if ($R_CDS){
    print "The R value used for an CDS TIS is                       : $R_CDS\n";
} else {
    #Choose default value for R
    $R_CDS = 0.15;
}
if ($R_3){
    print "The R value used for an 3'UTR TIS is                     : $R_3\n";
} else {
    #Choose default value for R
    $R_3 = 0.05;
}
if ($R_ntr){
    print "The R value used for an no translation TIS is            : $R_ntr\n";
} else {
    #Choose default value for R
    $R_ntr = 0.05;
}
if ($transcriptfilter){
    print "The transcripts are filtered based on            : $transcriptfilter\n";
} else {
    #Choose a transcript filter
    $transcriptfilter = 'none';
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
my ($ensemblversion,$species,$ens_db,$IGENOMES_ROOT,$cores,$mean_length_fastq1,$mean_length_fastq2) = get_arguments($dsn_results,$us_results,$pw_results);

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

#Old mouse assembly = NCBIM37, new one is GRCm38. Old human assembly is GRCh37, new one is GRCh38.
my $assembly = ($species eq "mouse" && $ensemblversion >= 70 ) ? "GRCm38"
: ($species eq "mouse" && $ensemblversion < 70 ) ? "NCBIM37"
: ($species eq "human" && $ensemblversion >= 76) ? "GRCh38"
: ($species eq "human" && $ensemblversion < 76) ? "GRCh37"
: ($species eq "arabidopsis") ? "TAIR10"
: (uc($species) eq "FRUITFLY" && $ensemblversion < 79) ? "BDGP5"
: (uc($species) eq "FRUITFLY" && $ensemblversion >= 79) ? "BDGP6" : "";

my $chrom_file = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Annotation/Genes/ChromInfo.txt";
my $BIN_chrom_dir = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Sequence/Chromosomes_BIN";

#############################################
## Peak Calling on exon_covered transcripts #
##                                          #
#############################################

# Start time
my $start = time;

print "\nGet analysis id...";
## Get Analysis ID
my $id = get_id($dsn_results,$us_results,$pw_results,$local_max,$min_count_aTIS,$R_aTIS,$min_count_5,$R_5,$min_count_CDS,$R_CDS,$min_count_3,$R_3,$min_count_ntr,$R_ntr,$transcriptfilter);
print "                                       : ".$id." \n";

print "Get chromosomes... \n";
## Get chromosomes based on seq_region_id ##
my $chrs = get_chrs($dsn_ENS,$us_ENS,$pw_ENS,get_chr_sizes(),$assembly);

# Create binary chromosomes if they don't exist
print "Checking/Creating binary chrom files ...\n";
if (!-d "$BIN_chrom_dir") {
    create_BIN_chromosomes($BIN_chrom_dir,$cores,$chrs,$work_dir,$TMP);
}

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
        my $table = "TIS_".$id;
        open TMP, "+>> ".$TMP."/".$table."_".$chr."_tmp.csv" or die $!;

        # seq_region_id per chromosome
        my $seq_region_id = $chrs->{$chr}{'seq_region_id'};

        # Get transcriptsand all reads
        my($trs,$CHX_for,$LTM_for,$CHX_rev,$LTM_rev) = get_transcripts_and_reads($seq_region_id,$chr);

        # Match reads to trancripts
        $trs = match_reads_to_transcripts($trs,$CHX_for,$CHX_rev,$LTM_for,$LTM_rev);

        # Loop over transcripts
        foreach my $tr_id (sort {$a <=> $b} keys %{$trs}){

            #Get transcript specific reads and translation data
            my($trans,$exons,$tr_seq,$LTM_reads,$CHX_reads,$trs) = get_transcript_translation_data_and_overlapping_reads($dbh,$tr_id,$trs,$seq_region_id,$chr,$CHX_for,$CHX_rev,$LTM_for,$LTM_rev);

            #Get position checked peaks from LTM reads
            my $TIS = analyze_transcript_peaks($trans,$exons,$tr_seq,$LTM_reads,$CHX_reads,$tr_id,$trs,$local_max,$chr,$min_count_aTIS,$min_count_5,$min_count_CDS,$min_count_3,$min_count_ntr,$R_aTIS,$R_5,$R_CDS,$R_3,$R_ntr,$mean_length_fastq1,$mean_length_fastq2);

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

### Analyze transcript peaks ###

sub analyze_transcript_peaks{

    #Catch
    my $trans               =   $_[0];
    my $exons               =   $_[1];
    my $tr_seq              =   $_[2];
    my $LTM_reads           =   $_[3];
    my $CHX_reads           =   $_[4];
    my $tr_id               =   $_[5];
    my $trs                 =   $_[6];
    my $local_max           =   $_[7];
    my $chr                 =   $_[8];
    my $min_count_aTIS      =   $_[9];
    my $min_count_5         =   $_[10];
    my $min_count_CDS       =   $_[11];
    my $min_count_3         =   $_[12];
    my $min_count_ntr       =   $_[13];
    my $R_aTIS              =   $_[14];
    my $R_5                 =   $_[15];
    my $R_CDS               =   $_[16];
    my $R_3                 =   $_[17];
    my $R_ntr               =   $_[18];
    my $mean_length_fastq1  =   $_[19];
    my $mean_length_fastq2  =   $_[20];

    #Init
    my $peaks           =   {};
    my $CHX_cdna_reads  =   {};
    my $TIS             =   {};
    my $peak            =   "Y";
    my $q               =   "N";
    my $tr              =   $trs->{$tr_id};
    my $strand          =   $tr->{'seq_region_strand'};
    my $end_ex_rank     =   $tr->{'end_exon_rank'};
    my $start_ex_rank   =   $tr->{'start_exon_rank'};
    my ($ext_cdn,$cdna_pos,$peak_pos,$start_cdn,$peak_shift,$count,$Rchx,$Rltm,$max,$min,$start,$stop,$i,$rank,$annotation,$CHX_peak_pos,$exon_rank);
    $local_max = ($local_max * 3); # 1 codon = 3 bps.

    # Total number of LTM and CHX reads in transcript
    my $total_LTM = 0;
    my $total_CHX = 0;


    #Only analyze if $LTM_reads not empty
    if(keys %{$LTM_reads}){

        foreach my $key (keys %{$CHX_reads}){
            $total_CHX += $CHX_reads->{$key}{'count'};

            #Get peak cDNA positions for CHX_reads
            $rank                                           =   $CHX_reads->{$key}{'exon_rank'};
            $cdna_pos                                       =   cdna_position($key,$exons,$rank);
            $CHX_cdna_reads->{$cdna_pos}                    =   $CHX_reads->{$key};
            $CHX_cdna_reads->{$cdna_pos}->{'cdna_pos'}      =   $cdna_pos;
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
                    $cdna_pos   =   $cdna_pos - 1;
                    $peak_pos = peak_position($cdna_pos,$exons);
                }elsif($strand eq '-1'){
                    $cdna_pos   =   $cdna_pos + 1;
                    $peak_pos = peak_position($cdna_pos,$exons);
                }

            }elsif(substr($ext_cdn,2,3) =~ m/[ACTG]TG|A[ACTG]G|AT[ACTG]/){
                $start_cdn = substr($ext_cdn,2,3);
                $peak_shift = '+1';
                if($strand eq '1'){
                    $cdna_pos   =   $cdna_pos + 1;
                    $peak_pos = peak_position($cdna_pos,$exons);
                }elsif($strand eq '-1'){
                    $cdna_pos   =   $cdna_pos - 1;
                    $peak_pos   =   peak_position($cdna_pos,$exons);
                }
            }else{
                $start_cdn = $ext_cdn;
                $peak_shift ='0';
                $peak = "N";
            }

            # Get rid of peaks at the -1 and ( transcript length +1) transcript positions.
            if($cdna_pos == 0){
                next;
            }elsif($cdna_pos > length($tr_seq)){
                next;
            }

            #Only save true TIS
            if($peak eq "Y"){

                #If already peak on that position, combine
                if(exists $peaks->{$cdna_pos}){

                    $peaks->{$cdna_pos}->{'peak_shift'}   =         $peaks->{$cdna_pos}->{'peak_shift'}." ".$peak_shift;
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

        # Run over all peaks to get annotation

        # If empty translation ==> biotype without translation info
        if(@{$trans}){

            #Get aTIS from first exon
            $start = ($strand eq '1') ? $exons->{$start_ex_rank}{'e_start'} + ${$trans}[2] -1 : $exons->{$start_ex_rank}{'e_end'} - ${$trans}[2] + 1;

            #Keep track of start position
            $tr->{'aTIS'}       =   $start;
            $tr->{'aTIS_call'}  =   "False";

            #Get stop codon position
            $stop =  ($strand eq '1') ?  $exons->{$end_ex_rank}{'e_start'} + ${$trans}[3] -1 : $exons->{$end_ex_rank}{'e_end'} - ${$trans}[3] +1;

            # Compare start peak position with aTis position
            # Compare start peak position with transcript start position
            # Annotate peak position (5'UTR,aTIS,CDS,3'UTR)
            foreach my $key (sort {$a <=> $b} keys %{$peaks}){

                $peaks->{$key}->{'biotype'}     =   $tr->{'biotype'};
                $peaks->{$key}->{'aTIS_call'}   =   "NA";
                $peaks->{$key}->{'stable_id'}   =   $tr->{'stable_id'};
                $cdna_pos                       =   $peaks->{$key}{'cdna_pos'};

                if($strand eq '1'){

                    $peaks->{$key}->{'dist_to_atis'}            =   $cdna_pos - $start;
                    $peaks->{$key}->{'dist_to_trans_start'}     =   $cdna_pos;

                    if($cdna_pos < $start){
                        $peaks->{$key}->{'annotation'}          =   "5UTR";
                        $peaks->{$key}->{'min_count'}           =   $min_count_5;
                        $peaks->{$key}->{'R'}                   =   $R_5;
                        $peaks->{$key}->{'exon'}                =   "NA";
                    }elsif($cdna_pos == $start){
                        $peaks->{$key}->{'annotation'}          =   "aTIS";
                        $peaks->{$key}->{'min_count'}           =   $min_count_aTIS;
                        $peaks->{$key}->{'R'}                   =   $R_aTIS;
                        $peaks->{$key}->{'aTIS_call'}           =   "True";
                        $tr->{'aTIS_call'}                      =   "True";
                        $peaks->{$key}->{'exon'}                =   "1";
                    }elsif($cdna_pos > $stop){
                        $peaks->{$key}->{'annotation'}          =   "3UTR";
                        $peaks->{$key}->{'min_count'}           =   $min_count_3;
                        $peaks->{$key}->{'R'}                   =   $R_3;
                        $peaks->{$key}->{'exon'}                =   "NA";
                    }else{
                        $peaks->{$key}->{'annotation'}          =   "CDS";
                        $peaks->{$key}->{'min_count'}           =   $min_count_CDS;
                        $peaks->{$key}->{'R'}                   =   $R_CDS;

                        #Get exon rank for TIS
                        foreach my $ex (keys %{$exons}){

                            if ($peaks->{$key}{'start'} >= $exons->{$ex}{'seq_region_start'} && $peaks->{$key}{'start'} <= $exons->{$ex}{'seq_region_end'}){
                                $exon_rank = $exons->{$ex}{'rank'};
                            }
                        }
                        $peaks->{$key}->{'exon'}                =   $exon_rank - $start_ex_rank + 1;
                    }
                }elsif($strand eq '-1'){

                    $peaks->{$key}->{'dist_to_atis'}            =   $start - $cdna_pos;
                    $peaks->{$key}->{'dist_to_trans_start'}     =   length($tr_seq) - $cdna_pos;

                    if($cdna_pos > $start){
                        $peaks->{$key}->{'annotation'}          =   "5UTR";
                        $peaks->{$key}->{'min_count'}           =   $min_count_5;
                        $peaks->{$key}->{'R'}                   =   $R_5;
                        $peaks->{$key}->{'exon'}                =   "NA";
                    }elsif($cdna_pos == $start){
                        $peaks->{$key}->{'annotation'}          =   "aTIS";
                        $peaks->{$key}->{'min_count'}           =   $min_count_aTIS;
                        $peaks->{$key}->{'R'}                   =   $R_aTIS;
                        $peaks->{$key}->{'aTIS_call'}           =   "True";
                        $tr->{'aTIS_call'}                      =   "True";
                        $peaks->{$key}->{'exon'}                =   "1";
                    }elsif($cdna_pos < $stop){
                        $peaks->{$key}->{'annotation'}          =   "3UTR";
                        $peaks->{$key}->{'min_count'}           =   $min_count_3;
                        $peaks->{$key}->{'R'}                   =   $R_3;
                        $peaks->{$key}->{'exon'}                =   "NA";
                    }else{
                        $peaks->{$key}{'annotation'}            =   "CDS";
                        $peaks->{$key}->{'min_count'}           =   $min_count_CDS;
                        $peaks->{$key}->{'R'}                   =   $R_CDS;

                        #Get exon rank for TIS
                        foreach my $ex (keys %{$exons}){

                            if ($peaks->{$key}{'start'} >= $exons->{$ex}{'seq_region_start'} && $peaks->{$key}{'start'} <= $exons->{$ex}{'seq_region_end'}){
                                $exon_rank = $exons->{$ex}{'rank'};
                            }
                        }
                        $peaks->{$key}->{'exon'}                =   $exon_rank - $start_ex_rank + 1;
                    }
                }
            }
        }else{

            # Compare start peak position with transcript start position
            foreach my $key (sort {$a <=> $b} keys %{$peaks}){

                $peaks->{$key}{'annotation'}            =   "ntr";
                $peaks->{$key}{'dist_to_atis'}          =   "NA";
                $peaks->{$key}{'biotype'}               =   $tr->{'biotype'};
                $peaks->{$key}->{'stable_id'}           =   $tr->{'stable_id'};
                $cdna_pos                               =   $peaks->{$key}{'cdna_pos'};
                $peaks->{$key}{'dist_to_trans_start'}   =   ($strand eq '1') ? $cdna_pos : length($tr_seq) - $cdna_pos;
                $peaks->{$key}->{'min_count'}           =   $min_count_ntr;
                $peaks->{$key}->{'R'}                   =   $R_ntr;
                $peaks->{$key}->{'aTIS_call'}           =   "NA";
                $peaks->{$key}->{'exon'}                =   "NA";

            }
        }

        # Run over all TIS peaks categorically to get true TISes
        foreach my $key (sort {$a <=> $b} keys %{$peaks}){

            $count              =   $peaks->{$key}{'count'};
            $peak_pos           =   $peaks->{$key}{'start'};
            $annotation         =   $peaks->{$key}{'annotation'};
            $q                  =   "N";
            my @CHX_split       =   ();
            my $CHX_split_total =   0;

            # A TIS has a minimal count of $min_count
            if($count >= $peaks->{$key}{'min_count'} || ($count < $peaks->{$key}{'min_count'} && $peaks->{$key}{'annotation'} eq "aTIS")){

                #If aTIS < min_count -> aTIS_call
                if ($count < $peaks->{$key}{'min_count'} && $peaks->{$key}{'annotation'} eq "aTIS"){
                    $peaks->{$key}->{'aTIS_call'} = "False";
                }

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

                # If $q = "Y",  peak not local max within $local_max basepairs | check for false aTIS
                 if($q eq "Y"){
                     if($annotation eq "aTIS"){
                         $peaks->{$key}->{'aTIS_call'} = "False";
                     }else{
                        next;
                     }
                 }

                #Calculate RLTM - RCHX
                #Rk = (Xk/(Nk/Mk)) x 10 (k = LTM, CHX), Xk number of reads on that position in data k, Nk total number of reads for transcript, Mk mean length of mapped ribo profile.

                #Combine CHX peaks number of ribo profiles
                @CHX_split = split(/ /,$peaks->{$key}->{'peak_shift'});

                #Calculate number of CHX reads
                if($strand eq '1'){
                    foreach(@CHX_split){
                        if (exists $CHX_cdna_reads->{$key - $_}{'count'}){
                            $CHX_split_total += $CHX_cdna_reads->{$key - $_}{'count'};
                        }
                    }
                }elsif($strand eq '-1'){
                    foreach(@CHX_split){
                        if (exists $CHX_cdna_reads->{$key + $_}{'count'}){
                            $CHX_split_total += $CHX_cdna_reads->{$key + $_}{'count'};
                        }
                    }
                }

                # Calculate Rltm
                $Rltm = ($count / ($total_LTM / $mean_length_fastq2)) * 10;
                if($CHX_split_total == 0){
                    $Rchx=0;
                }else{
                    $Rchx = ($CHX_split_total/ ($total_CHX / $mean_length_fastq1) );
                    $Rchx = $Rchx*10;
                }
                my $R_test = $Rltm - $Rchx;

                # A TIS must have an RLTM - RCHX value > $R | Check aTIS
                if(($Rltm - $Rchx) < $peaks->{$key}{'R'}){
                    if($annotation eq "aTIS"){
                        $peaks->{$key}->{'aTIS_call'} = "False";
                    }else{
                        next;
                    }
                }

                # Store Peak in Hash TIS
                $TIS->{$key}                        =   $peaks->{$key};
                $TIS->{$key}->{'Rltm_min_Rchx'}     =   $Rltm-$Rchx;
                $TIS->{$key}->{'transcript_id'}     =   $tr_id;

            }else{next;}

        }

        # Check if aTIS available for protein_coding, otherwise add annotated TIS position to TIS hash
        if(@{$trans}){
            if ($tr->{'aTIS_call'} eq "False"){

                my $genomic_atis_position = peak_position($tr->{'aTIS'},$exons);

                $TIS->{$tr->{'aTIS'}}->{'transcript_id'}        =   $tr_id;
                $TIS->{$tr->{'aTIS'}}->{'stable_id'}            =   $tr->{'stable_id'};
                $TIS->{$tr->{'aTIS'}}->{'biotype'}              =   $tr->{'biotype'};
                $TIS->{$tr->{'aTIS'}}->{'chr'}                  =   $chr;
                $TIS->{$tr->{'aTIS'}}->{'strand'}               =   $strand;
                $TIS->{$tr->{'aTIS'}}->{'start'}                =   $genomic_atis_position;
                $TIS->{$tr->{'aTIS'}}->{'dist_to_trans_start'}  =   ($strand eq '1') ?  $tr->{'aTIS'} : length($tr_seq) - $tr->{'aTIS'};
                $TIS->{$tr->{'aTIS'}}->{'dist_to_atis'}         =   "0";
                $TIS->{$tr->{'aTIS'}}->{'annotation'}           =   "aTIS";
                $TIS->{$tr->{'aTIS'}}->{'start_cdn'}            =   ($strand eq '1') ? substr($tr_seq,$tr->{'aTIS'} -1,3) : revdnacomp(substr($tr_seq,$tr->{'aTIS'}-3,3));
                $TIS->{$tr->{'aTIS'}}->{'peak_shift'}           =   "NA";
                $TIS->{$tr->{'aTIS'}}->{'count'}                =   "NA";
                $TIS->{$tr->{'aTIS'}}->{'Rltm_min_Rchx'}        =   "NA";
                $TIS->{$tr->{'aTIS'}}->{'aTIS_call'}            =   "no_data";
                $TIS->{$tr->{'aTIS'}}->{'exon'}                 =   "1";
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

            foreach my $key (@keys){
                $CHX_reads{$key}{'exon_rank'} = $ex->{'rank'};
            }

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

            foreach my $key (@keys){
                $CHX_reads{$key}{'exon_rank'} = $ex->{'rank'};
            }

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

### Get Transcripts from transcript calling ###

sub get_transcripts_and_reads {

    # Catch
    my $seq_region_id   =   $_[0];
    my $chr             =   $_[1];

    # Init
    my $dbh = dbh($dsn_results,$us_results,$pw_results);
    my $trs         = {};
    my $CHX_for     = {};
    my $LTM_for     = {};
    my $CHX_rev     = {};
    my $LTM_rev     = {};

    # Get transcripts
    my $query;
    if (uc($transcriptfilter) eq "NONE") {
        $query = "SELECT transcript_id,seq_region_id,seq_region_strand,seq_region_start,seq_region_end,read_counts,biotype,stable_id from tr_translation where seq_region_id = $seq_region_id AND (biotype = 'protein_coding' or (biotype != 'protein_coding' and exon_coverage = 'Yes'))";
    } elsif (uc($transcriptfilter) eq "CANONICAL") {
        $query = "SELECT transcript_id,seq_region_id,seq_region_strand,seq_region_start,seq_region_end,read_counts,biotype,stable_id from tr_translation where seq_region_id = $seq_region_id AND (biotype = 'protein_coding' or (biotype != 'protein_coding' and exon_coverage = 'Yes')) and canonical = 'Yes'";
    } elsif (uc($transcriptfilter) eq "CCDS") {
         $query = "SELECT transcript_id,seq_region_id,seq_region_strand,seq_region_start,seq_region_end,read_counts,biotype,stable_id from tr_translation where seq_region_id = $seq_region_id AND (biotype = 'protein_coding' or (biotype != 'protein_coding' and exon_coverage = 'Yes')) and ccds <> 'No'";
    }
    my $sth = $dbh->prepare($query);
	$sth->execute();
	$trs = $sth->fetchall_hashref('transcript_id');

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
	return($trs,$CHX_for,$LTM_for,$CHX_rev,$LTM_rev);
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
    my $table   =   "TIS_".$id;

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
    `exon` varchar(128) NOT NULL default 'NA',
    `aTIS_call` varchar(128) NOT NULL default 'NA',
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
        system("sqlite3 -separator , ".$sqlite_db." \".import ".$TMP."/".$table."_".$chr."_tmp.csv ".$table."\"")== 0 or die "system failed: $?";
    }

    #Move to galaxy history
    system ("mv ".$sqlite_db." ".$out_sqlite);
    #system("rm -rf ".$sqlite_db);

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

        print TMP $TIS->{$key}{'transcript_id'}.",".$TIS->{$key}{'stable_id'}.",".$TIS->{$key}{'biotype'}.",".$TIS->{$key}{'chr'}.",".$TIS->{$key}{'strand'}.",".$TIS->{$key}{'start'}.",".$TIS->{$key}{'dist_to_trans_start'}.",".$TIS->{$key}{'dist_to_atis'}.",".$TIS->{$key}{'annotation'}.",".$TIS->{$key}->{'exon'}.",".$TIS->{$key}->{'aTIS_call'}.",".$TIS->{$key}{'start_cdn'}.",".$TIS->{$key}{'peak_shift'}.",".$TIS->{$key}{'count'}.",".$TIS->{$key}{'Rltm_min_Rchx'}."\n";
    }
}

### GET CHRs ###

sub get_chrs {

    # Catch
    my $db          =   $_[0];
    my $us          =   $_[1];
    my $pw          =   $_[2];
    my $chr_sizes   =   $_[3];
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

### Create Bin Chromosomes ##

sub create_BIN_chromosomes {

    # Catch
    my $BIN_chrom_dir   =   $_[0];
    my $cores           =   $_[1];
    my $chrs            =   $_[2];
    my $work_dir        =   $_[3];
    my $TMP             =   $_[4];

    # Create BIN_CHR directory
    system ("mkdir -p ".$BIN_chrom_dir);

    # Create binary chrom files
    ## Init multi core
    my $pm = new Parallel::ForkManager($cores);

    foreach my $chr (keys %$chrs){

        ## Start parallel process
        $pm->start and next;

        open (CHR,"<".$IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Sequence/Chromosomes/".$chr.".fa") || die "Cannot open chr fasta input\n";
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

    $query = "select value from `arguments` where variable = \'ens_db\'";
    $sth = $dbh->prepare($query);
	$sth->execute();
	my $ens_db = $sth->fetch()->[0];

    $query = "select value from `arguments` where variable = \'igenomes_root\'";
    $sth = $dbh->prepare($query);
	$sth->execute();
	my $igenomes_root = $sth->fetch()->[0];

    $query = "select value from `arguments` where variable = \'nr_of_cores\'";
    $sth = $dbh->prepare($query);
	$sth->execute();
	my $nr_of_cores = $sth->fetch()->[0];

    $query = "select value from `arguments` where variable = \'mean_length_fastq1\'";
    $sth = $dbh->prepare($query);
	$sth->execute();
	my $mean_length_fastq1 = $sth->fetch()->[0];

    $query = "select value from `arguments` where variable = \'mean_length_fastq2\'";
    $sth = $dbh->prepare($query);
	$sth->execute();
	my $mean_length_fastq2 = $sth->fetch()->[0];

    # Return input variables
    return($ensemblversion,$species,$ens_db,$igenomes_root,$nr_of_cores,$mean_length_fastq1,$mean_length_fastq2);

}

### Get Analysis ID

sub get_id{

    # Catch
    my $dsn                 =   $_[0];
    my $us                  =   $_[1];
    my $pw                  =   $_[2];
    my $local_max           =   $_[3];
    my $min_count_aTIS      =   $_[4];
    my $R_aTIS              =   $_[5];
    my $min_count_5         =   $_[6];
    my $R_5                 =   $_[7];
    my $min_count_CDS       =   $_[8];
    my $R_CDS               =   $_[9];
    my $min_count_3         =   $_[10];
    my $R_3                 =   $_[11];
    my $min_count_ntr  =   $_[12];
    my $R_ntr          =   $_[13];

    # Init
    my $dbh = dbh($dsn,$us,$pw);

    # Create table
    my $query = "CREATE TABLE IF NOT EXISTS `TIS_overview` (
    `ID` INTEGER primary key,
    `local_max` int(10) NOT NULL default '',
    `min_count_aTIS` int(10) NOT NULL default '',
    `R_aTis` decimal(11,8) NOT NULL default '',
    `min_count_5UTR` int(10) NOT NULL default '',
    `R_5UTR` decimal(11,8) NOT NULL default '',
    `min_count_CDS` int(10) NOT NULL default '',
    `R_CDS` decimal(11,8) NOT NULL default '',
    `min_count_3UTR` int(10) NOT NULL default '',
    `R_3UTR` decimal(11,8) NOT NULL default '',
    `min_count_ntr` int(10) NOT NULL default '',
    `R_ntr` decimal(11,8) NOT NULL default '',
    `SNP` varchar(20) default '',
    `filter` varchar(20) default '')";
    $dbh->do($query);

    # Add parameters to overview table and get ID
    print "INSERT INTO `TIS_overview`(local_max,min_count_aTIS,R_aTIS,`min_count_5UTR`,`R_5UTR`,min_count_CDS,R_CDS,`min_count_3UTR`,`R_3UTR`,min_count_ntr,R_ntr,filter) VALUES($local_max,$min_count_aTIS,$R_aTIS,$min_count_5,$R_5,$min_count_CDS,$R_CDS,$min_count_3,$R_3,$min_count_ntr,$R_ntr,'$transcriptfilter')\n";
    $dbh->do("INSERT INTO `TIS_overview`(local_max,min_count_aTIS,R_aTIS,`min_count_5UTR`,`R_5UTR`,min_count_CDS,R_CDS,`min_count_3UTR`,`R_3UTR`,min_count_ntr,R_ntr,filter) VALUES($local_max,$min_count_aTIS,$R_aTIS,$min_count_5,$R_5,$min_count_CDS,$R_CDS,$min_count_3,$R_3,$min_count_ntr,$R_ntr,'$transcriptfilter')");
    my $id= $dbh->func('last_insert_rowid');

    # Return
    return($id);
}

### get peak genomic position based on cDNA position ###

sub peak_position {

    #Catch
    my $cdna_pos        =   $_[0];
    my $exons           =   $_[1];

    #Init
    my $peak_pos;

    #Run over all exons
    foreach my $exon (sort {$a <=> $b} keys %{$exons}){

        #Check if in which exon cDNA position is located
        if ($cdna_pos ~~ [$exons->{$exon}{"e_start"}..$exons->{$exon}{"e_end"}]) {

            #Get peak genomic position based on cdna position and exon start
            $peak_pos = $exons->{$exon}{'seq_region_start'} - $exons->{$exon}{"e_start"} + $cdna_pos;
        }
    }

    #Return
    return($peak_pos);

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
