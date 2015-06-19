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

##############
##Command-line
# ./assembly_intergenic_sORFs.pl  --sqliteRES SQLite/results_STAR12.db --sqliteENS /data/jeroenc/ribo_prof_micropeptide_pipeline/SQLite/ENS_mmu_72.db --cores 20 --tis_ids 1 --tmp /data/jeroenc/ribo_prof_micropeptide_pipeline/tmp/ --min_length 10 --max_length 100 --R .05 --coverage 0.75 --igenomes_root /data/igenomes/Mus_musculus/Ensembl/GRCm38/


# get the command line arguments
my ($cores,$tmpfolder,$R,$snp,$min,$max,$coverage,$resultDB,$ensDB,$work_dir,$IGENOMES_ROOT,$tis_ids,$out_sqlite);

GetOptions(
"sqliteRES=s"=>\$resultDB,                  # The sqlite DB holding all RIBO-pipeline results,                      mandatory argument
"sqliteENS=s"=>\$ensDB,                     # The sqlite DB holding the ensembl information,                        mandatory argument
"dir:s"=>\$work_dir,                        # Path to the working directory,                                        optional  argument
"cores=i"=>\$cores,                         # Number of cores to use for Bowtie Mapping,                            mandatory argument
"tmp:s" =>\$tmpfolder,                      # Folder where temporary files are stored,                              optional  argument (default = $TMP env setting
"min_length=i"=>\$min,                      # Minimal sORF length in Amino Acids                                    mandatory argument
"max_length=i"=>\$max,                      # Maximal sORF length in Amino Acids                                    mandatory argument
"R:f" => \$R,                               # The Rltm - Rchx value calculated based on both CHX and LTM data       optional argument (default .05)
"coverage:f" => \$coverage,                 # The minimal coverage for sORF sequence based on total ribo_reads      optional argument (default .75)
"snp:s"=>\$snp,                             # The snp calling algorithm applied                                     optional argument (default "NO", others can be "samtools", "GATK")
"igenomes_root=s" =>\$IGENOMES_ROOT,        # IGENOMES ROOT FOLDER                                                  mandatory argument
"tis_ids=s"  =>\$tis_ids,                   # list of analysis ids                                                  mandatory argument
"out_sqlite=s"=>\$out_sqlite,               # The sqlite DB holding all the output                                  mandatory argument
);

my $CWD             = getcwd;
my $HOME            = $ENV{'HOME'};
$IGENOMES_ROOT      = ($ENV{'IGENOMES_ROOT'}) ? $ENV{'IGENOMES_ROOT'} : $IGENOMES_ROOT;
print "The following igenomes folder is used                    : $IGENOMES_ROOT\n";
my $TMP             = ($ENV{'TMP'}) ? $ENV{'TMP'} : ($tmpfolder) ? $tmpfolder : "$CWD/tmp"; # (1) get the TMP environment variable, (2) get the $tmpfolder variable, (3) get current_working_dir/tmp
print "The following tmpfolder is used                          : $TMP\n";

#Check if tmpfolder exists, if not create it...
if (!-d "$TMP") {
    system ("mkdir ". $TMP);
}

# Check and create protein output directory
#my $prot_dir = "$CWD/protein";
#if (!-d "$prot_dir") {
#   system ("mkdir -p ".$prot_dir);
#}

#comment on these
if ($resultDB){
    print "SqliteDB used is                                         : $CWD/$resultDB\n";
} else {
    die "\nDon't forget to pass the SQLite DB using the --sqliteDB argument!\n\n";
}
if ($ensDB){
    print "EnsDB used is                                            : $ensDB\n";
} else {
    die "\nDon't forget to pass the Ensembl custom DB using the --sqliteENS argument!\n\n";
}
if ($out_sqlite){
    print "Output SqliteDB used is                                  : $CWD/$out_sqlite\n";
} else {
    $out_sqlite = $resultDB;
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
if ($min){
    print "Minimal sORF length (AA)                                 : $min\n";
} else {
    die "\nDon't forget to pass the minimal sORF length using the --min_length or -min argument!\n\n";
}
if ($max){
    print "Maximal sORF length (AA)                                 : $max\n";
} else {
    die "\nDon't forget to pass the maximal sORF length using the --max_length or -max argument!\n\n";
}
if ($R){
    print "The R value used is                                      : $R\n";
} else {
    #Choose default value for R
    $R = 0.05;
}
if ($coverage){
    print "The minimal sORF coverage based on total ribo_reads      : $coverage\n";
} else {
    #Choose default value for coverage
    $coverage = 0.75;
}
if ($snp){
    print "The snp calling algorithm used is                        : $snp\n";
} else {
    #Choose default value for mapper
    $snp = "NO";
    print "The snp calling algorithm used is                        : $snp\n";
}

# DB settings
# Sqlite Riboseq
my $db_results  = $resultDB;
my $dsn_results = "DBI:SQLite:dbname=$db_results";
my $us_results  = "";
my $pw_results  = "";

# Sqlite Ensembl
my $db_ENS  = $ensDB;
my $dsn_ENS = "DBI:SQLite:dbname=$db_ENS";
my $us_ENS  = "";
my $pw_ENS  = "";

#Get the input variables (Check nog welke eruit mogen)
my $dbh_results = dbh($dsn_results,$us_results,$pw_results);
my ($ensemblversion,$species,$Mreads) = get_input_vars($dbh_results);

 print "Total number of CHX reads in millions of reads           : $Mreads\n";

#Conversion for species terminology
my $spec = ($species eq "mouse") ? "Mus_musculus" : ($species eq "human") ? "Homo_sapiens" : "";
my $spec_short = ($species eq "mouse") ? "mmu" : ($species eq "human") ? "hsa" : "";
#Old mouse assembly = NCBIM37, new one is GRCm38
my $assembly = ($species eq "mouse" && $ensemblversion >= 70 ) ? "GRCm38"
: ($species eq "mouse" && $ensemblversion < 70 ) ? "NCBIM37"
: ($species eq "human") ? "GRCh37"
: "";

my $chrom_file = $IGENOMES_ROOT."/Annotation/Genes/ChromInfo.txt";
my $BIN_chrom_dir = $IGENOMES_ROOT."/Sequence/Chromosomes_BIN";

#################################################
##                                             ##
## ASSEMBLY OF INTERGENIC TRANSLATION PRODUCTS ##
##                                             ##
#################################################

# Start time
my $start = time;

# Get the analysis_id that corresponds to the TIS-calling input parameters
my $idsref = get_analysis_ids($dbh_results,$tis_ids); #$tis_ids is input variable

print "\nGet chromosomes... \n";
## Get chromosomes based on seq_region_id ##
my $chrs = get_chrs($dsn_ENS,$us_ENS,$pw_ENS,$chrom_file,$assembly);

# Create binary chromosomes if they don't exist
print "Checking/Creating binary chrom files ...\n";
if (!-d "$BIN_chrom_dir") {
    create_BIN_chromosomes($BIN_chrom_dir,$cores,$chrs,$TMP);
}

my $analysis_id;
#Loop over all selected analysis_ids
print "Assembly of translation products with FPKM and coverage analysis ...\n";
foreach $analysis_id (@$idsref) {
    
    print "Processing analysis_id $analysis_id ...\n";
    construct_intergenic_sORFs($chrs,$TMP,$analysis_id,$min,$max,$Mreads,$R,$coverage);

    ## Store in DB
    print "   Storing all sORFs in DB \n";
    store_in_db($dsn_results,$us_results,$pw_results,$analysis_id,$work_dir,$chrs,$TMP);
}

#Move to galaxy history
system ("mv ".$resultDB." ".$out_sqlite);

# End time
print "   DONE! \n";
my $end = time - $start;
printf("runtime assembly: %02d:%02d:%02d\n\n",int($end/3600), int(($end % 3600)/60), int($end % 60));

############
# THE SUBS #
############

### Construct Intergenic sORFs ###

sub construct_intergenic_sORFs{

    #Catch
    my $chrs        =   $_[0];
    my $tmp         =   $_[1];
    my $analysis_id =   $_[2];
    my $min_length  =   $_[3];
    my $max_length  =   $_[4];
    my $Mreads      =   $_[5];
    my $R           =   $_[6];
    my $coverage    =   $_[7];
    
    # Init multi core
    my $pm = new Parallel::ForkManager($cores);
    
    ## Loop over all chromosomes
    foreach my $chr (sort keys %{$chrs}){
        
        ### Start parallel process
        $pm->start and next;
        
        ### DBH per process
        my $dbh = dbh($dsn_results,$us_results,$pw_results);
        
        ### Output fasta and db_csv file per process
        open TMP_db, "+>>".$TMP."/".$analysis_id."_".$chr."_sORF_tmp.csv" or die $!;

        #my $temp_fasta = $TMP."/".$analysis."_".$chr."_sORF_tmp.fasta";
        #open TMP_fasta, "+>>".$temp_fasta or die $!;
        
        ## Get intergenic_ids per chr from transcript calling
        my $ig_sORF_starts = get_intergenic_sORF_starts_per_chromosome($dbh,$analysis_id,$chr);
        
        # Init
        my ($sORF_seq,$tr_sORF_seq,$tmp_sORF_seq,$ig_sORF_id,$TIS,$strand,$start_codon,$peak_shift,$count,$sORF_length,$FPKM,$coverage_sORF,$sorf_end,$R_sORF,$sorf_chr,$sorf_begin,$sorf_strand,$annotation,$upstream_gene_distance,$downstream_gene_distance);
        my $ig_sorfs = {};
        my $max_DNA_length = ($max_length*3) -1;
            
        foreach my $ig_sORF_start ( keys %{$ig_sORF_starts}){
            #my %tr_SNPs;
            $sORF_seq = '';
            $tmp_sORF_seq = '';
            
            # Split transcript and TIS
            $ig_sORF_id                     =   $ig_sORF_starts->{$ig_sORF_start}{'intergene_id'};
            $TIS                            =   $ig_sORF_starts->{$ig_sORF_start}{'start'};
            $strand                         =   $ig_sORF_starts->{$ig_sORF_start}{'strand'};
            $start_codon                    =   $ig_sORF_starts->{$ig_sORF_start}{'start_codon'};
            $peak_shift                     =   $ig_sORF_starts->{$ig_sORF_start}{'peak_shift'};;
            $count                          =   $ig_sORF_starts->{$ig_sORF_start}{'count'};
            $downstream_gene_distance       =   $ig_sORF_starts->{$ig_sORF_start}{'downstream_gene_distance'};
            $upstream_gene_distance         =   $ig_sORF_starts->{$ig_sORF_start}{'upstream_gene_distance'};
            
            # Create max_length sequence from TIS (strand specific)
            $tmp_sORF_seq = ($strand eq '1') ? get_sequence($chr,$TIS,$TIS + $max_DNA_length) : revdnacomp(get_sequence($chr,$TIS - $max_DNA_length,$TIS));
            
            # Create translated sequence
            ($tr_sORF_seq,$sORF_seq,$sORF_length) = translate($tmp_sORF_seq,$min,$max);
            
            #next if no sorf sequence
            if($sORF_seq eq 'X'){next;}
            
            #Save in ig_sorfs
            $sorf_end   =   ($strand eq '1') ? $TIS + (($sORF_length * 3) -1) : $TIS;
            $sorf_begin =   ($strand eq '1') ? $TIS : $TIS - (($sORF_length * 3) +1);
                
            $ig_sorfs->{$ig_sORF_start}{'intergene_id'}             =   $ig_sORF_id;
            $ig_sorfs->{$ig_sORF_start}{'sorf_chr'}                 =   $chr;
            $ig_sorfs->{$ig_sORF_start}{'sorf_begin'}               =   $sorf_begin;
            $ig_sorfs->{$ig_sORF_start}{'FPKM_begin'}               =   $sorf_begin - 15;
            $ig_sorfs->{$ig_sORF_start}{'sorf_end'}                 =   $sorf_end;
            $ig_sorfs->{$ig_sORF_start}{'FPKM_end'}                 =   $sorf_end +15;
            $ig_sorfs->{$ig_sORF_start}{'sorf_strand'}              =   $strand;
            $ig_sorfs->{$ig_sORF_start}{'start_codon'}              =   $start_codon;
            $ig_sorfs->{$ig_sORF_start}{'peak_shift'}               =   $peak_shift;
            $ig_sorfs->{$ig_sORF_start}{'count'}                    =   $count;
            $ig_sorfs->{$ig_sORF_start}{'sorf_length'}              =   $sORF_length;
            $ig_sorfs->{$ig_sORF_start}{'sorf_seq'}                 =   $sORF_seq;
            $ig_sorfs->{$ig_sORF_start}{'tr_sorf_seq'}              =   $tr_sORF_seq;
            $ig_sorfs->{$ig_sORF_start}{'downstream_gene_distance'} =   $downstream_gene_distance;
            $ig_sorfs->{$ig_sORF_start}{'upstream_gene_distance'}   =   $upstream_gene_distance;
            $ig_sorfs->{$ig_sORF_start}{'annotation'}               =   'intergenic';
            
            #For LTM reads start from sorf_begin -1 to calculate total LTM
            #Otherwize, if sORF has only LTM_reads on sorf_begin -1 position total_LTM = 0 -> illegal division
            $ig_sorfs->{$ig_sORF_start}{'LTM_begin'}   =    ($strand eq '1') ?  $sorf_begin - 1                          :   $ig_sorfs->{$ig_sORF_start}{'sorf_begin'};
            $ig_sorfs->{$ig_sORF_start}{'LTM_end'}     =    ($strand eq '1') ?  $ig_sorfs->{$ig_sORF_start}{'sorf_end'}  :   $sorf_end + 1;
        }
        
        #Get CHX and LTM reads
        my ($CHX_for,$LTM_for,$CHX_rev,$LTM_rev) = get_reads($chr);
        
        # Match reads to intergenic sORFs
        $ig_sorfs = match_reads_to_intergenic_sORFs($ig_sorfs,$CHX_for,$CHX_rev,$LTM_for,$LTM_rev);
        
        # Loop over sORFs in $ig_sORFs 
        foreach my $ig_sorf_id (keys %{$ig_sorfs}){
        
            # Get sORF specific reads
            my ($LTM_reads,$CHX_reads,$FPKM_CHX_reads) = get_overlapping_reads($ig_sorf_id,$ig_sorfs,$chr,$CHX_for,$CHX_rev,$LTM_for,$LTM_rev);
            
            # Calculate FPKM, coverage and R
            ($FPKM,$coverage_sORF,$R_sORF) = calculate_FPKM_and_coverage($ig_sorfs,$ig_sorf_id,$Mreads,$LTM_reads,$CHX_reads,$FPKM_CHX_reads);
        
            #Check if R value and coverage are higher or equal as $R and $coverage
            if($R_sORF < $R){next;}
            if($coverage_sORF < $coverage){next;}
            
            # Save sORF in csv
            $ig_sORF_id                 =   $ig_sorfs->{$ig_sorf_id}{'intergene_id'};
            $sorf_chr                   =   $ig_sorfs->{$ig_sorf_id}{'sorf_chr'};
            $sorf_begin                 =   $ig_sorfs->{$ig_sorf_id}{'sorf_begin'};
            $sorf_end                   =   $ig_sorfs->{$ig_sorf_id}{'sorf_end'};
            $sorf_strand                =   $ig_sorfs->{$ig_sorf_id}{'sorf_strand'};
            $start_codon                =   $ig_sorfs->{$ig_sorf_id}{'start_codon'};
            $peak_shift                 =   $ig_sorfs->{$ig_sorf_id}{'peak_shift'};
            $count                      =   $ig_sorfs->{$ig_sorf_id}{'count'};
            $sORF_length                =   $ig_sorfs->{$ig_sorf_id}{'sorf_length'};
            $sORF_seq                   =   $ig_sorfs->{$ig_sorf_id}{'sorf_seq'};
            $tr_sORF_seq                =   $ig_sorfs->{$ig_sorf_id}{'tr_sorf_seq'};
            $downstream_gene_distance   =   $ig_sorfs->{$ig_sorf_id}{'downstream_gene_distance'};
            $upstream_gene_distance     =   $ig_sorfs->{$ig_sorf_id}{'upstream_gene_distance'};
            $annotation                 =   $ig_sorfs->{$ig_sorf_id}{'annotation'};
            
            # Check how to implement
            my $out_cnt        =   '1';
            # Should becom $_->{'SNP'}
            my $snp_tmp = '0'; 
            
            print TMP_db $ig_sORF_id.",".$sorf_chr.",".$sorf_strand.",".$sorf_begin.",".$sorf_end.",".$sORF_length.",".$start_codon.",".$downstream_gene_distance.",".$upstream_gene_distance.",".$annotation.",".$peak_shift.",".$count.",".$R_sORF.",".$coverage_sORF.",".$FPKM.",".$snp_tmp.",".$sORF_seq.",".$tr_sORF_seq."\n";
            
            #print TMP_fasta ">generic|".$ig_sORF_id."_".$sorf_chr."_".$sorf_begin."_".$sorf_strand."_".$out_cnt."|".$start_codon."_".$downstream_gene_distance."_".$upstream_gene_distance."_".$annotation."_".$peak_shift."_".$count."_".$R_sORF."_".$coverage_sORF."_".$FPKM."_".$snp_tmp."\n".$tr_sORF_seq."\n";
        
        }    
        ### Finish childs
        print "     * Finished translating chromosome ".$chr."\n";
        $dbh->disconnect();
        $pm->finish;
    }
    #Waiting for all childs to finish
    $pm->wait_all_children();

    # Gather all fasta files into one file in the protein subdirectory
    #my $temp_fasta_all = $TMP."/".$run_name."_sORF_tmp.fasta";
    #system("touch ".$temp_fasta_all);
    
    #foreach my $chr (sort keys %{$chrs}){
    #    my $temp_fasta = $TMP."/".$run_name."_".$chr."_sORF_tmp.fasta";
    #    system ("cat ".$temp_fasta." >>". $temp_fasta_all);
    #}
    #Remove _tmp.fasta files
    #system("rm -rf ".$TMP."/".$run_name."_*_sORF_tmp.fasta");
    #system ("mv ".$temp_fasta_all." ".$prot_dir."/".$run_name."_sORF.fasta");
}

### Store in DB ###

sub store_in_db{
    
    # Catch
    my $dsn         =   $_[0];
    my $us          =   $_[1];
    my $pw          =   $_[2];
    my $id          =   $_[3];
    my $work_dir    =   $_[4];
    my $chrs        =   $_[5];
    my $TMP         =   $_[6];
    
    # Init
    my $dbh     =   dbh($dsn,$us,$pw);
    my $table   =   "TIS_intergenic_".$id."_transcripts";
    
    # Create table 
    my $query = "CREATE TABLE IF NOT EXISTS `".$table."` (
    `intergene_id` varchar(128) NOT NULL default '',
    `chr` char(50) NOT NULL default '',
    `strand` int(2) NOT NULL default '',
    `sorf_begin` int(10) NOT NULL default '',
    `sorf_end` int(10) NOT NULL default '',
    `sorf_length` int(10) NOT NULL default '',
    `start_codon` varchar(128) NOT NULL default '',
    `downstream_gene_distance` int(10) NOT NULL default '',
    `upstream_gene_distance` int(10) NOT NULL default 'NA',
    `annotation` varchar(128) NOT NULL default 'NA',
    `peak_shift` int(2) NOT NULL default '',
    `count` float default NULL,
    `Rltm_min_Rchx` decimal(11,8) NOT NULL default '0',
    `coverage` decimal(11,8) NOT NULL default '0',
    `RPKM` decimal(11,8) NOT NULL default '0',
    `SNP_` decimal(11,8) NOT NULL default '0',
    `tr_seq` TEXT NOT NULL default '',
    `aa_seq` TEXT NOT NULL default '' )"  ;
    $dbh->do($query);  
    
    # Store
    foreach my $chr (sort keys %{$chrs}){
        system("sqlite3 -separator , ".$work_dir."/".$resultDB." \".import ".$TMP."/".$id."_".$chr."_sORF_tmp.csv ".$table."\"")== 0 or die "system failed: $?";
    }
    
    #Disconnect dbh
    $dbh->disconnect();
    
    # Unlink tmp csv files
    foreach my $chr (sort keys %{$chrs}){
        unlink $TMP."/".$id."_".$chr."_sORF_tmp.csv";
    }

}

### Calculate FPKM and coverage ###

sub calculate_FPKM_and_coverage{
    
    # Catch
    my $ig_sorfs        =   $_[0];
    my $ig_sorf_id      =   $_[1];
    my $Mreads          =   $_[2];
    my $LTM_reads       =   $_[3];
    my $CHX_reads       =   $_[4];
    my $FPKM_CHX_reads  =   $_[5];
    
    # Init
    my $sorf_end    =   $ig_sorfs->{$ig_sorf_id}{'sorf_end'};
    my $sorf_begin  =   $ig_sorfs->{$ig_sorf_id}{'sorf_begin'};
    my $count       =   $ig_sorfs->{$ig_sorf_id}{'count'};
    my ($FPKM,$coverage_sORF,$read_begin,$read_end,$Rltm,$Rchx,$R_sORF);
    my $sorf_length =   $sorf_end-$sorf_begin +1;
    my $dbh = dbh($dsn_results,$us_results,$pw_results);
    my $hits = 0;
    my $hits_total = 0;
    my $covered_positions = 0;
    my $position = 0;
    my $coverage = {};
    
    # Total number of LTM and CHX reads in intergene
    my $total_LTM = 0;
    my $total_CHX = 0;

    foreach my $key (keys %{$CHX_reads}){
        $total_CHX += $CHX_reads->{$key}{'count'};
    }
    
    foreach my $key (keys %{$LTM_reads}){
        $total_LTM += $LTM_reads->{$key}{'count'};
    }
    # if ($total_LTM == 0){print Dumper($ig_sorfs->{$ig_sorf_id});}
    
    # Run over CHX ribo_profiles
    foreach my $key (keys %{$FPKM_CHX_reads}) {
        
        $position   =   $FPKM_CHX_reads->{$key}{'start'};
        #$hits       =   $FPKM_CHX_reads->{$key}{'count'};
        $read_begin =   $position - 14;
        $read_end   =   $position +14;
        #$hits_total =   $hits + $hits_total;
        
        #Run over read and create ribo_read coverage
        for (my $i=$read_begin;$i<$read_end;$i++){
            $coverage->{$i} = '1';
        }					
    }
    
    #Calculate percent coverage and FPKM for sORF
    for (my $i= $sorf_begin; $i<$sorf_end+1; $i++){
        if(defined $coverage->{$i}){
            if($coverage->{$i} == 1){
                $covered_positions++;	
            }
        }
    }
    $coverage_sORF = $covered_positions/$sorf_length;
    
    #$FPKM = $hits_total / (($sorf_length / 1000) * $Mreads);
    $FPKM = $total_CHX / (($sorf_length / 1000) * $Mreads);
    
    #Calculate RLTM - RCHX
    #Rk = (Xk/Nk) x 10 (k = LTM, CHX), Xk number of reads on that position in data k, Nk total number of reads for transcript.
    
    # Calculate Rltm
    $Rltm = ($count / $total_LTM) * 10;
    if(exists $CHX_reads->{$sorf_begin}{'count'}){
        $Rchx = ($CHX_reads->{$sorf_begin}{'count'}/$total_CHX);
        $Rchx = $Rchx*10;
    }else{
        $Rchx=0;
    }
    $R_sORF = $Rltm - $Rchx;
    
    #Return
    return($FPKM,$coverage_sORF,$R_sORF);
}

### translation subroutine ###
    
sub translate {
    
    # Add SNPs if available 
    #       -> Watch sequence length (if indel, add sequence)
    #       -> Watch positions, sort SNP Hash, if dels, downstream position changes
    
    #Catch
    my $tmp_sORF_seq    =   $_[0];
    my $min             =   $_[1];
    my $max             =   $_[2];
    
    #Init
    my ($i,$triplet,$tr_sORF_seq,$sORF_seq);
    my $AA_seq      =   '';
    my $DNA_seq     =   '';
    my $AA          =   '';
    my $sORF_length =   0;
       

    for($i=0;$i<=length($tmp_sORF_seq)-2;$i+=3){
        $sORF_length++;
        $triplet = substr($tmp_sORF_seq,$i,3);
        
        #Translate DNA to AA
        $AA = DNAtoAA($triplet);
        $AA_seq     =   $AA_seq . $AA;
        $DNA_seq    =   $DNA_seq . $triplet;
        
        #Stop translation if STOP-codon is reached
        if ($AA eq '*'){last;}
    }
    
    # Check length of sORF_AA_seq
    if($sORF_length >= 10){
        if ($AA_seq =~ /\*$/ && $AA_seq ne '') {
           
            #Replace near-cognate start to cognate methionine...
            $AA_seq = (substr($AA_seq,0,1) ne 'M') ? 'M'.substr($AA_seq,1) : $AA_seq;
            
            $tr_sORF_seq    =   $AA_seq;
            $sORF_seq       =   $DNA_seq;
            
            #Return
            return($tr_sORF_seq,$sORF_seq,$sORF_length);
        }
    }
    
    #Bad sORF sequence (too short or no stop codon)
    $tr_sORF_seq    =   'X';
    $sORF_seq       =   'X';
    
    #Return
    return($tr_sORF_seq,$sORF_seq,$sORF_length);
}
    
### Get intergene overlapping reads ###

sub get_overlapping_reads{
    
    # Catch
    my $ig_sorf_id      =   $_[0];
	my $ig_sorfs        =   $_[1];
    my $chr             =   $_[2];
    my $CHX_for         =   $_[3];
    my $CHX_rev         =   $_[4];
    my $LTM_for         =   $_[5];
    my $LTM_rev         =   $_[6];
    
    #Init
    my %LTM_reads = ();
    my %CHX_reads = ();
    my %FPKM_CHX_reads = ();
    
    my $sorf = $ig_sorfs->{$ig_sorf_id};
    my $strand = $sorf->{'sorf_strand'};
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
    if($ig_sorfs->{$ig_sorf_id}{'CHX'}){
        my @keys = @{$ig_sorfs->{$ig_sorf_id}{'CHX'}};
        @CHX_reads{@keys} = @{$CHX_all}{@keys};
    }
    
    if($ig_sorfs->{$ig_sorf_id}{'LTM'}){
        my @keys = @{$ig_sorfs->{$ig_sorf_id}{'LTM'}};
        @LTM_reads{@keys} = @{$LTM_all}{@keys};
    }
    
    if($ig_sorfs->{$ig_sorf_id}{'FPKM_CHX'}){
        my @keys = @{$ig_sorfs->{$ig_sorf_id}{'FPKM_CHX'}};
        @FPKM_CHX_reads{@keys} = @{$CHX_all}{@keys};
    }
    
    #Return
    return(\%LTM_reads,\%CHX_reads,\%FPKM_CHX_reads);
}

### Match reads to intergenes ##

sub match_reads_to_intergenic_sORFs{
    
    #Catch
    my $ig_sorfs    =   $_[0];
    my $CHX_for     =   $_[1];
    my $CHX_rev     =   $_[2];
    my $LTM_for     =   $_[3];
    my $LTM_rev     =   $_[4];
    
    #Init
    my @window = ();
    my ($ig_sorf_for,$ig_sorf_rev,$ig_sorf_for_LTM,$ig_sorf_rev_LTM,$ig_sorf_for_FPKM,$ig_sorf_rev_FPKM);         
    
    #Split ig_sorfs in forward and reverse arrays
    foreach my $ig_sorf_id (sort { $ig_sorfs->{$a}{'sorf_begin'} <=> $ig_sorfs->{$b}{'sorf_begin'} } keys %{$ig_sorfs}){
        if ($ig_sorfs->{$ig_sorf_id}{'sorf_strand'} eq '1'){
            push (@$ig_sorf_for,$ig_sorf_id);
            push (@$ig_sorf_for_LTM,$ig_sorf_id);
            push (@$ig_sorf_for_FPKM,$ig_sorf_id);
        }else{
            push (@$ig_sorf_rev,$ig_sorf_id);
            push (@$ig_sorf_rev_LTM,$ig_sorf_id);
            push (@$ig_sorf_rev_FPKM,$ig_sorf_id);
        }
    }
    
    # Loop over CHX_forward
    foreach my $key (sort {$a <=> $b} keys %{$CHX_for}){
        
        # Push all ig_sorf_ids to @window where ig_sorf_start < window_pos
        foreach my $ig_sorf_for_id (@$ig_sorf_for){
            if($ig_sorfs->{$ig_sorf_for_id}{'sorf_begin'} <= $key){
                push(@window,$ig_sorf_for_id);
            }else{last;}
        }
        # Get rid of ig_sorf_for elements already in @window
        @$ig_sorf_for = grep { $ig_sorfs->{$_}{'sorf_begin'} > $key} @$ig_sorf_for;
        
        # Get rid of ig_sorf_ids in @$window where ig_end < window_pos
        @window = grep { $ig_sorfs->{$_}{'sorf_end'} >= $key} @window;
        
        # Loop over window and add read position to window_ig_sorfs
        foreach my $window_id (@window){
            push(@{$ig_sorfs->{$window_id}{'CHX'}},$key);
        }
    }
    
    #Empty @window
    @window = ();
    
    # Loop over CHX_reverse
    foreach my $key (sort {$a <=> $b} keys %{$CHX_rev}){
        
        # Push all ig_sorf_ids to @window where ig_sorf_start < window_pos
        foreach my $ig_sorf_rev_id (@$ig_sorf_rev){
            if($ig_sorfs->{$ig_sorf_rev_id}{'sorf_begin'} <= $key){
                push(@window,$ig_sorf_rev_id);
            }else{last;}
        }
        # Get rid of ig_sorf_for elements already in @window
        @$ig_sorf_rev = grep { $ig_sorfs->{$_}{'sorf_begin'} > $key} @$ig_sorf_rev;
        
        # Get rid of ig_sorf_ids in @$window where ig_sorf_end < window_pos
        @window = grep { $ig_sorfs->{$_}{'sorf_end'} >= $key} @window;
        
        # Loop over window and add read position to window_ig_sorfs
        foreach my $window_id (@window){
            push(@{$ig_sorfs->{$window_id}{'CHX'}},$key);
        }
    }
    
    #Empty @window
    @window = ();
    
    # Loop over LTM_forward
    foreach my $key (sort {$a <=> $b} keys %{$LTM_for}){
        
        # Push all ig_sorf_ids to @window where ig_sorf_start < window_pos
        foreach my $ig_sorf_for_LTM_id (@$ig_sorf_for_LTM){
            if($ig_sorfs->{$ig_sorf_for_LTM_id}{'LTM_begin'} <= $key){
                push(@window,$ig_sorf_for_LTM_id);
            }else{last;}
        }
        # Get rid of ig_sorf_for elements already in @window
        @$ig_sorf_for_LTM = grep { $ig_sorfs->{$_}{'LTM_begin'} > $key} @$ig_sorf_for_LTM;
        
        # Get rid of ig_sorf_ids in @$window where ig_sorf_end < window_pos
        @window = grep { $ig_sorfs->{$_}{'LTM_end'} >= $key} @window;
        
        # Loop over window and add read position to window_ig_sorfs
        foreach my $window_id (@window){
            push(@{$ig_sorfs->{$window_id}{'LTM'}},$key);
        }
    }
    
    #Empty @window
    @window = ();
    
    # Loop over LTM_reverse
    foreach my $key (sort {$a <=> $b} keys %{$LTM_rev}){
        
        # Push all ig_sorf_ids to @window where ig_sorf_start < window_pos
        foreach my $ig_sorf_rev_LTM_id (@$ig_sorf_rev_LTM){
            if($ig_sorfs->{$ig_sorf_rev_LTM_id}{'LTM_begin'} <= $key){
                push(@window,$ig_sorf_rev_LTM_id);
            }else{last;}
        }
        # Get rid of ig_sorf_for elements already in @window
        @$ig_sorf_rev_LTM = grep { $ig_sorfs->{$_}{'LTM_begin'} > $key} @$ig_sorf_rev_LTM;
        
        # Get rid of ig_sorf_ids in @$window where ig_sorf_end < window_pos
        @window = grep { $ig_sorfs->{$_}{'LTM_end'} >= $key} @window;
        
        # Loop over window and add read position to window_ig_sorfs
        foreach my $window_id (@window){
            push(@{$ig_sorfs->{$window_id}{'LTM'}},$key);
        }
    }
    
    ## Get CHX reads from sorf_begin -15 (=FPKM_begin) to sorf_end +15 (=FPKM_end) for FPKM and coverage analysis
    
    #Empty @window
    @window = ();
    
    # Loop over CHX_forward
    foreach my $key (sort {$a <=> $b} keys %{$CHX_for}){
        
        # Push all ig_sorf_ids to @window where ig_sorf_start < window_pos
        foreach my $ig_sorf_for_id (@$ig_sorf_for_FPKM){
            if($ig_sorfs->{$ig_sorf_for_id}{'FPKM_begin'} <= $key){
                push(@window,$ig_sorf_for_id);
            }else{last;}
        }
        # Get rid of ig_sorf_for elements already in @window
        @$ig_sorf_for_FPKM = grep { $ig_sorfs->{$_}{'FPKM_begin'} > $key} @$ig_sorf_for_FPKM;
        
        # Get rid of ig_sorf_ids in @$window where ig_end < window_pos
        @window = grep { $ig_sorfs->{$_}{'FPKM_end'} >= $key} @window;
        
        # Loop over window and add read position to window_ig_sorfs
        foreach my $window_id (@window){
            push(@{$ig_sorfs->{$window_id}{'FPKM_CHX'}},$key);
        }
    }
    
    #Empty @window
    @window = ();
    
    # Loop over CHX_reverse
    foreach my $key (sort {$a <=> $b} keys %{$CHX_rev}){
        
        # Push all ig_sorf_ids to @window where ig_sorf_start < window_pos
        foreach my $ig_sorf_rev_id (@$ig_sorf_rev_FPKM){
            if($ig_sorfs->{$ig_sorf_rev_id}{'FPKM_begin'} <= $key){
                push(@window,$ig_sorf_rev_id);
            }else{last;}
        }
        # Get rid of ig_sorf_for elements already in @window
        @$ig_sorf_rev_FPKM = grep { $ig_sorfs->{$_}{'FPKM_begin'} > $key} @$ig_sorf_rev_FPKM;
        
        # Get rid of ig_sorf_ids in @$window where ig_sorf_end < window_pos
        @window = grep { $ig_sorfs->{$_}{'FPKM_end'} >= $key} @window;
        
        # Loop over window and add read position to window_ig_sorfs
        foreach my $window_id (@window){
            push(@{$ig_sorfs->{$window_id}{'FPKM_CHX'}},$key);
        }
    }
    
    #Return
    return($ig_sorfs);
}

### Get reads ###

sub get_reads {
    
    # Catch
    my $chr             =   $_[0];
    
    # Init
    my $dbh         = dbh($dsn_results,$us_results,$pw_results);
    my $CHX_for     = {};
    my $LTM_for     = {};
    my $CHX_rev     = {};
    my $LTM_rev     = {};
    
    # Get Reads
    my $query = "SELECT * FROM count_fastq1 WHERE chr = '$chr' and strand = '1'";
    my $sth = $dbh->prepare($query);
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
	return($CHX_for,$LTM_for,$CHX_rev,$LTM_rev);    
}

### get intergenic sORF start positions from TIS-calling table ###

sub get_intergenic_sORF_starts_per_chromosome {
    
    # Catch
    my $dbh         =   $_[0];
    my $id          =   $_[1];
    my $chr         =   $_[2];
    
    # Init
    my $ig_sORF_starts = {};
    
    # Get transcripts
    my $query = "SELECT intergene_id||'_'||start as intergene_start,intergene_id,chr,strand,start,annotation,upstream_gene_distance,downstream_gene_distance,start_codon,peak_shift,count from TIS_intergenic_".$id." WHERE chr = '".$chr."'";
	my $sth = $dbh->prepare($query);
	$sth->execute();
	$ig_sORF_starts = $sth->fetchall_hashref('intergene_start');
    
	# Return
	return($ig_sORF_starts);
}

### get snp information from SQLite into hash ###

sub get_SNPs_per_chromosome {
    
    # Catch
    my $dbh         = $_[0];
    my $chr         = $_[1];
    my $snp         = $_[2];
    my $CHX_lane    = $_[3];
    
    # Init
    my $SNPs = {};
    
   	# Get SNPs based on chr
	my $query = "select chr,pos,ref,alt,af from snp_".$snp."_lane".$CHX_lane." where chr = '".$chr."'";
    
	my $sth = $dbh->prepare($query);
	$sth->execute();
	$SNPs = $sth->fetchall_hashref('pos');
	
	# Return
	return($SNPs);
}

### Fetch SNPs per intergenic sORF

sub fetch_SNPs_sORF {
    
    # Catch
    my $sorf_begin  =   $_[0];
    my $sorf_end    =   $_[1];
    my %SNPs        =   %{$_[2]};
    
    # Init
    my (@slice_list,%sORF_SNPs_all,%sORF_SNPs_def,$key,$value);
    
    # Take hash slice based on sORF_begin/sORF_end
    @slice_list = ($sorf_begin..$sorf_end);
    @sORF_SNPs_all{@slice_list} = @SNPs{@slice_list};
    foreach $key (keys %sORF_SNPs_all) {
        # Only print defined keys
        if ($sORF_SNPs_all{$key}) { $sORF_SNPs_def{$key} =  $sORF_SNPs_all{$key} }
    }
    
    #Return
    return(%sORF_SNPs_def);
}

### CREATE BIN CHROMS ###

sub create_BIN_chromosomes {
    
    # Catch
    my $BIN_chrom_dir = $_[0];
    my $cores = $_[1];
    my $chr_sizes = $_[2];
    my $TMP = $_[3];
    
    # Create BIN_CHR directory
    system ("mkdir -p ".$BIN_chrom_dir);
    
    # Create binary chrom files
    ## Init multi core
    my $pm = new Parallel::ForkManager($cores);
    print "   Using ".$cores." core(s)\n   ---------------\n";
    
    foreach my $chr (keys %$chr_sizes){
        
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

### GET_INPUT_VARS ###

sub get_input_vars {
    # Catch
    my $dbh_results = $_[0];
    
    my ($query,$sth);
    
    # Get input variables
    $query = "select value from arguments where variable = \'ensembl_version\'";
    $sth = $dbh_results->prepare($query);
	$sth->execute();
	my $ensemblversion = $sth->fetch()->[0];
    
    $query = "select value from arguments where variable = \'species\'";
    $sth = $dbh_results->prepare($query);
	$sth->execute();
	my $species = $sth->fetch()->[0];
    
    $query = "select total from statistics where sample like \'%fastq1%\' and type = \'genomic\'";
    $sth = $dbh_results->prepare($query);
	$sth->execute();
	my $Mreads = $sth->fetch()->[0];
    $Mreads = $Mreads/1000000;
    
    # Return input variables
    return($ensemblversion,$species,$Mreads);
}

### Get the analysis ids that need to be processed

sub get_analysis_ids {
    
    # Catch
    my $dbh    = $_[0];
    my $ids_in = $_[1]; #Either comma separated list of identifiers or "all"
    
    
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
    
    return $idsref;
}

### DNA to AA ###

sub DNAtoAA{
    
    #Catch
    my $triplet = $_[0];
    
    #Init
    my %AA1 = (
    'TTT','F','TTC','F','TTA','L','TTG','L','TCT','S','TCC','S','TCA','S','TCG','S',
    'TAT','Y','TAC','Y','TAA','*','TAG','*','TGT','C','TGC','C','TGA','*','TGG','W',
    'CTT','L','CTC','L','CTA','L','CTG','L','CCT','P','CCC','P','CCA','P','CCG','P',
    'CAT','H','CAC','H','CAA','Q','CAG','Q','CGT','R','CGC','R','CGA','R','CGG','R',
    'ATT','I','ATC','I','ATA','I','ATG','M','ACT','T','ACC','T','ACA','T','ACG','T',
    'AAT','N','AAC','N','AAA','K','AAG','K','AGT','S','AGC','S','AGA','R','AGG','R',
    'GTT','V','GTC','V','GTA','V','GTG','V','GCT','A','GCC','A','GCA','A','GCG','A',
    'GAT','D','GAC','D','GAA','E','GAG','E','GGT','G','GGC','G','GGA','G','GGG','G',
    'AAN','X','ATN','X','ACN','T','AGN','X','ANN','X','ANA','X','ANG','X','ANC','X',
    'ANT','X','TAN','X','TTN','X','TCN','S','TGN','X','TNN','X','TNA','X','TNG','X',
    'TNC','X','TNT','X','CAN','X','CTN','L','CCN','P','CGN','R','CNN','X','CNA','X',
    'CNG','X','CNC','X','CNT','X','GAN','X','GTN','V','GCN','A','GGN','G','GNN','X',
    'GNA','X','GNG','X','GNC','X','GNT','X','NAN','X','NAA','X','NAG','X','NAC','X',
    'NAT','X','NTN','X','NTA','X','NTG','X','NTC','X','NTT','X','NGN','X','NGA','X',
    'NGG','X','NGC','X','NGT','X','NCN','X','NCA','X','NCG','X','NCC','X','NCT','X',
    'NNN','X','NNA','X','NNG','X','NNC','X','NNT','X',);
    
    my $AA = $AA1{uc($triplet)};
    
    #Return
    return($AA);
    
}

### get reverse complement sequence ###

sub revdnacomp {
    my $dna = shift;
    my $revcomp = reverse($dna);
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
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