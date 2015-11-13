#!/usr/bin/perl -w

$|=1;

use strict;
use warnings;
use DBI;
use Data::Dumper;
use Getopt::Long;
use v5.10;
use Parallel::ForkManager;
use Cwd;


##############
##Command-line
# ./2_assembly.pl  --sqliteRES sqlite_results_DB --tis_ids list_of_analysis_ids (--snp NO --dir /data/RIBO_runs/RIBO_Ingolia_GerbenM/ --tmp /data/RIBO_runs/RIBO_Ingolia_GerbenM/tmp/ --localmax 1 --mincount_aTIS 10 --R_aTIS .05 and the other mincount/R values --out_sqlite $out_sqlite)
#  --igenomes_root IGENOMES_ROOT --cores nr_of_cores are taken from the arguments table (SQLiteDB, after mapping script)

# For GALAXY
# ./2_assembly.pl


# get the command line arguments
my ($resultDB,$ensDB,$work_dir,$seqFileName,$run_name,$species,$ensemblversion,$cores,$mapper,$readlength,$readtype,$tmpfolder,$adaptorSeq,$unique,$CHX_lane,$LTM_lane,$local_max,$min_count_aTIS,$Rltm_minus_Rchx_aTIS,$min_count_5UTR,$Rltm_minus_Rchx_5UTR,$min_count_3UTR,$Rltm_minus_Rchx_3UTR,$min_count_CDS,$Rltm_minus_Rchx_CDS,$min_count_no_translation,$Rltm_minus_Rchx_no_translation,$snp,$IGENOMES_ROOT,$tis_ids,$out_sqlite);


GetOptions(
"sqliteRES=s"=>\$resultDB,                  # The sqlite DB holding all RIBO-pipeline results,                      mandatory argument
#Get sqliteENS from arguments table
#"sqliteENS=s"=>\$ensDB,                     # The sqlite DB holding the ensembl information,
"dir:s"=>\$work_dir,                        # Path to the working directory,                                        optional  argument
#Get nr_of_cores from arguments table
#"cores=i"=>\$cores,                         # Number of cores to use for Bowtie Mapping,                            mandatory argument
"tmp:s" =>\$tmpfolder,                      # Folder where temporary files are stored,                              optional  argument (default = $TMP env setting)
"localmax:i" =>\$local_max,                 # The range wherein the localmax is (e.g.: 1 means +/- one triplet)     optional  argument (default 1)
"mincount_aTIS:i" =>\$min_count_aTIS,       # The minimum count of riboseq profiles mapping to the TIS site         optional  argument (default 5)
"R_aTIS:f" => \$Rltm_minus_Rchx_aTIS,       # The Rltm - Rchx value calculated based on both CHX and LTM data       optional  argument (default .05)
"mincount_5UTR:i" =>\$min_count_5UTR,       # The minimum count of riboseq profiles mapping to the TIS site         optional  argument (default 10)
"R_5UTR:f" => \$Rltm_minus_Rchx_5UTR,       # The Rltm - Rchx value calculated based on both CHX and LTM data       optional  argument (default .05)
"mincount_3UTR:i" =>\$min_count_3UTR,       # The minimum count of riboseq profiles mapping to the TIS site         optional  argument (default 10)
"R_3UTR:f" => \$Rltm_minus_Rchx_3UTR,       # The Rltm - Rchx value calculated based on both CHX and LTM data       optional  argument (default .05)
"mincount_CDS:i" =>\$min_count_CDS,         # The minimum count of riboseq profiles mapping to the TIS site         optional  argument (default 15)
"R_CDS:f" => \$Rltm_minus_Rchx_CDS,         # The Rltm - Rchx value calculated based on both CHX and LTM data       optional  argument (default .15)
"mincount_no_translation:i" =>\$min_count_no_translation,     # The minimum count of riboseq profiles mapping to the TIS site         optional argument (default 10)
"R_no_translation:f" => \$Rltm_minus_Rchx_no_translation,     # The Rltm - Rchx value calculated based on both CHX and LTM data       optional argument (default .05)
"snp:s"=>\$snp,                             # The snp calling algorithm applied                                     optional  argument (default "NO", others can be "samtools", "samtools_dbSNP")
#Get igenomes root from ARGUMENTS table
#"igenomes_root=s" =>\$IGENOMES_ROOT,        # IGENOMES ROOT FOLDER                                                  mandatory argument
"tis_ids=s"  =>\$tis_ids,                   # list of analysis ids                                                  mandatory argument
"out_sqlite=s"=>\$out_sqlite,               # The sqlite DB holding all the output                                  optional argument (default, same as mandatory input $resultDB)

);

print "TIS_IDS= $tis_ids\n";

#open(STDERR, ">&STDOUT");
my $CWD             = getcwd;
my $HOME            = $ENV{'HOME'};
my $TMP             = ($ENV{'TMP'}) ? $ENV{'TMP'} : ($tmpfolder) ? $tmpfolder : "$CWD/tmp"; # (1) get the TMP environment variable, (2) get the $tmpfolder variable, (3) get current_working_dir/tmp
print "The following tmpfolder is used                          : $TMP\n";

#Check if tmpfolder exists, if not create it...
if (!-d "$TMP") {
    system ("mkdir ". $TMP);
}

# comment on these
if ($resultDB){
    print "SqliteDB used is                                         : $resultDB\n";
} else {
    die "\nDon't forget to pass the SQLite DB using the --sqliteDB argument!\n\n";
}
#if ($ensDB){
#    print "EnsDB used is                                            : $ensDB\n";
#} else {
#    die "\nDon't forget to pass the Ensembl custom DB using the --sqliteENS argument!\n\n";
#}
if ($out_sqlite){
    print "Output SqliteDB used is                                  : $out_sqlite\n";
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
if ($snp){
    print "The snp calling algorithm used is                        : $snp\n";
} else {
    #Choose default value
    $snp = "NO";
    print "snp information included                                 : $snp\n";
}
#Get nr_of_cores from sqliteDB-arguments table
#if ($cores){
#    print "Number of cores to use for Mapping                       : $cores\n";
#} else {
#    die "\nDon't forget to pass number of cores to use for mapping using the --cores or -c argument!\n\n";
#}
if ($local_max){
    print "The regio for local max is set to                        : $local_max\n";
} else {
    #Choose default value
    $local_max = 1;
}
####### ALL THE min_count and R values for the different categories (aTIS, 5UTR, 3UTR, CDS, and no_translation)
####### DEFAULT SETTINGS ARE:
#ID|CHX_lane|LTM_lane|ensembl_version|intergenic|species     |local_max|min_count_aTIS|R_aTis|min_count_5UTR|R_5UTR|min_count_CDS|R_CDS|min_count_3UTR|R_3UTR|min_count_no_trans|R_no_trans
#1 |lane3   |lane4   |72             |N         |Mus_musculus|1        |5             |0.05  |10            |0.05  |15           |0.15 |10            |0.05  |10                |0.05
if ($min_count_aTIS){
    print "The minimum coverage of a TIS is                         : $min_count_aTIS\n";
} else {
    #Choose default value
    $min_count_aTIS = 5;
}
if ($Rltm_minus_Rchx_aTIS){
    print "The aTIS R value used is                                 : $Rltm_minus_Rchx_aTIS\n";
} else {
    #Choose default value
    $Rltm_minus_Rchx_aTIS = 0.05;
}
if ($min_count_5UTR){
    print "The minimum coverage of a 5UTR is                        : $min_count_5UTR\n";
} else {
    #Choose default value
    $min_count_5UTR = 10;
}
if ($Rltm_minus_Rchx_5UTR){
    print "The 5UTR R value used is                                 : $Rltm_minus_Rchx_5UTR\n";
} else {
    #Choose default value
    $Rltm_minus_Rchx_5UTR = 0.05;
}
if ($min_count_3UTR){
    print "The minimum coverage of a 3UTR is                        : $min_count_3UTR\n";
} else {
    #Choose default value
    $min_count_3UTR = 10;
}
if ($Rltm_minus_Rchx_3UTR){
    print "The 3UTR R value used is                                 : $Rltm_minus_Rchx_3UTR\n";
} else {
    #Choose default value
    $Rltm_minus_Rchx_3UTR = 0.05;
}
if ($min_count_CDS){
    print "The minimum coverage of a CDS is                         : $min_count_CDS\n";
} else {
    #Choose default value
    $min_count_CDS = 15;
}
if ($Rltm_minus_Rchx_CDS){
    print "The CDS R value used is                                  : $Rltm_minus_Rchx_CDS\n";
} else {
    #Choose default value
    $Rltm_minus_Rchx_CDS = 0.15;
}
if ($min_count_no_translation){
    print "The minimum coverage of a no_translation is              : $min_count_no_translation\n";
} else {
    #Choose default value
    $min_count_no_translation = 10;
}
if ($Rltm_minus_Rchx_no_translation){
    print "The no_translation R value used is                       : $Rltm_minus_Rchx_no_translation\n";
} else {
    #Choose default value
    $Rltm_minus_Rchx_no_translation = 0.05;
}

#ADDED FOR TEST ISSUES
my $clusterPhoenix = "N";


# DB settings for results sqlite RIBOseq
my $db_results  = $resultDB;
my $dsn_results = "DBI:SQLite:dbname=$db_results";
my $us_results  = "";
my $pw_results  = "";

# Get the input variables
my $dbh_results = dbh($dsn_results,$us_results,$pw_results);
($run_name,$ensemblversion,$species,$mapper,$unique,$adaptorSeq,$readlength,$readtype,$IGENOMES_ROOT,$cores,$ensDB)=get_input_vars($dbh_results);

# DB settings for results sqlite Ensembl
my $db_ENS  = $ensDB;
my $dsn_ENS = "DBI:SQLite:dbname=$db_ENS";
my $us_ENS  = "";
my $pw_ENS  = "";


$IGENOMES_ROOT      = ($ENV{'IGENOMES_ROOT'}) ? $ENV{'IGENOMES_ROOT'} : $IGENOMES_ROOT;
print "The following igenomes folder is used                    : $IGENOMES_ROOT\n";

#Conversion for species terminology
my $spec = ($species eq "mouse") ? "Mus_musculus" : ($species eq "human") ? "Homo_sapiens" : ($species eq "arabidopsis") ? "Arabidopsis_thaliana" : ($species eq "fruitfly") ? "Drosophila_melanogaster" : "";
my $spec_short = ($species eq "mouse") ? "mmu" : ($species eq "human") ? "hsa" : ($species eq "arabidopsis") ? "ath" : ($species eq "fruitfly") ? "dme" : "";
#Old mouse assembly = NCBIM37, new one is GRCm38. Old human assembly = GRCh37. The new one is GRCh38.
my $assembly = ($species eq "mouse" && $ensemblversion >= 70 ) ? "GRCm38"
: ($species eq "mouse" && $ensemblversion < 70 ) ? "NCBIM37"
: ($species eq "human" && $ensemblversion >= 76) ? "GRCh38"
: ($species eq "human" && $ensemblversion < 76) ? "GRCh37"
: ($species eq "arabidopsis") ? "TAIR10"
: ($species eq "fruitfly") ? "BDGP5" : "";

#my $chromosome_sizes = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Annotation/Genes/ChromInfo.txt";
my $BIN_chrom_dir = ($clusterPhoenix eq "Y") ? $HOME."/igenomes_extra/".$spec."/Ensembl/".$assembly."/Sequence/Chromosomes_BIN" : $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Sequence/Chromosomes_BIN";


###################################
## ASSEMBLY OF TRANSLATION PRODUCTS
###################################

# Start time
my $start = time;

## Get chromosome sizes
print "Getting chromosome sizes  ...\n";
my $chr_sizes = get_chr_sizes();

# Create binary chromosomes if they don't exist
print "Checking/Creating binary chrom files ...\n";
if (!-d "$BIN_chrom_dir") {
    create_BIN_chromosomes($BIN_chrom_dir,$cores,$chr_sizes,$TMP);
}

# Get the analysis_id that corresponds to the TIS-calling input parameters
my $idsref = get_analysis_ids($dbh_results,$tis_ids); #$tis_ids is input variable
#my $analysis_id = get_analysis_id($dbh_results,$local_max,$min_count_aTIS,$Rltm_minus_Rchx_aTIS,$min_count_5UTR,$Rltm_minus_Rchx_5UTR,$min_count_3UTR,$Rltm_minus_Rchx_3UTR,$min_count_CDS,$Rltm_minus_Rchx_CDS,$min_count_no_translation,$Rltm_minus_Rchx_no_translation);


#print "analysis id is $analysis_id\n";


# Assembly of the translation products
print "Assembly of translation products ...\n";
my $annot = "Y";
my $fusion = "N";
my $novel = "N";
my $ins = "N";
my $del = "N";
my $transl_type = "total"; #total,splice
my $transl_RF = "1"; # 1,3,6

#FOR ALL INPUT PARAMS MENTIONED ABOVE
# CHECK IF SQLITE DB ARE AVAILABLE AND IF COMBI's ARE POSSIBLE
# IF NOT, SPIT OUTPUT TELLING THAT THIS OUTPUT IS MISSING

my $analysis_id;
#Loop over all selected analysis_ids
foreach $analysis_id (@$idsref) {
    print "Processing analysis_id $analysis_id ...\n";
    construct_trans_prod($annot,$fusion,$novel,$snp,$ins,$del,$transl_type,$transl_RF,$readtype,$work_dir,$TMP,$analysis_id);
    #Add snp info in TIS_OVERVIEW records
    update_TIS_overview($snp,$analysis_id);
}

#Move to galaxy history
system ("mv ".$resultDB." ".$out_sqlite);

# End time
my $end = time - $start;
printf("runtime assembly: %02d:%02d:%02d\n\n",int($end/3600), int(($end % 3600)/60), int($end % 60));


############
# THE SUBS #
############

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


### GET_INPUT_VARS ###

sub get_input_vars {
    # Catch
    my $dbh_results = $_[0];
    
    my ($query,$sth);
    
    # Get input variables
    $query = "select value from arguments where variable = \'run_name\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $run_name = $sth->fetch()->[0];
    
    $query = "select value from arguments where variable = \'ensembl_version\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $ensemblversion = $sth->fetch()->[0];
    
    $query = "select value from arguments where variable = \'species\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $species = $sth->fetch()->[0];
    
    $query = "select value from arguments where variable = \'mapper\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $mapper = $sth->fetch()->[0];
    
    $query = "select value from arguments where variable = \'unique\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $unique = $sth->fetch()->[0];
    
    $query = "select value from arguments where variable = \'adaptor\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $adaptorSeq = $sth->fetch()->[0];
    
    $query = "select value from arguments where variable = \'readlength\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $readlength = $sth->fetch()->[0];
    
    $query = "select value from arguments where variable = \'readtype\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $readtype = $sth->fetch()->[0];
    
    $query = "select value from arguments where variable = \'igenomes_root\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $IGENOMES_ROOT = $sth->fetch()->[0];
    
    $query = "select value from arguments where variable = \'nr_of_cores\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $nr_of_cores = $sth->fetch()->[0];
    
    $query = "select value from arguments where variable = \'ens_db\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $ensDB = $sth->fetch()->[0];
    
    # Return input variables
    return($run_name,$ensemblversion,$species,$mapper,$unique,$adaptorSeq,$readlength,$readtype,$IGENOMES_ROOT,$nr_of_cores,$ensDB);
}


### GET CHRs ###
sub get_chrs {
    
    # Catch
    my $dbh          =  $_[0];
    my $chr_file    =   $_[1];
    my $assembly    =   $_[2];
    
    # Init
    my $chrs    =   {};
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

### GET CHR SIZES FROM IGENOMES ###
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

### ASSEMBLY OF TRANS PRODS ###
sub construct_trans_prod {
    
    # Catch
    my $annot = $_[0];
    my $fusion = $_[1];
    my $novel = $_[2];
    my $snp = $_[3];
    my $ins = $_[4];
    my $del = $_[5];
    my $transl_type = $_[6];
    my $transl_RF = $_[7];
    my $readtype = $_[8];
    my $work_dir = $_[9];
    my $tmp = $_[10];
    my $analysis_id = $_[11];
    
    #print "analysis_id = $analysis_id\n";
    # Get chromosomes name,seq_region_id
    my $dbh_ENS = dbh($dsn_ENS,$us_ENS,$pw_ENS);
    my $chrs = get_chrs($dbh_ENS,get_chr_sizes(),$assembly);
    my %chrs = %{$chrs};
    $dbh_ENS->disconnect();
    
    # Init multi core
    my $pm = new Parallel::ForkManager($cores);
    print "   Using ".$cores." core(s)\n   ---------------\n";
    
    
    # Open each chromosome in seperate core
    foreach my $chr (sort keys %{$chrs}){
        
        #my @chrs = (1);
        #foreach my $chr (@chrs) {
        my $SNP_chr = {};
        
        ### Start parallel process
        $pm->start and next;
        
        ## Analyse Peaks per Transcript
        #print "\n\n   Starting translation on chromosome ".$chr."\n";
        
        ### DBH per process
        my $dbh_results = dbh($dsn_results,$us_results,$pw_results);
        my $dbh_ENS = dbh($dsn_ENS,$us_ENS,$pw_ENS);
        
        ### Output  db-file per process
        open TMP_db, "+>> ".$TMP."/".$chr."_tmp.csv" or die $!;
        
        ## Get SNPs per chr from SNP SQLite
        if ( $snp ne "NO" ) {
            $SNP_chr = get_SNPs_per_chromosome($dbh_results,$chr,$snp);
            
        }
        
        ## Get transcript_ids per chr from transcript calling
        my $transcript_starts = get_transcripts_per_chromosome($dbh_results,$chr,$analysis_id);
        #print Dumper($transcript_starts); exit;
        # Init
        my ($cnt, $transcript_start,$tr_id,$TIS,$exons,$strand,$start_codon,$dist_to_transcript_start,$dist_to_aTIS,$annotation,$aTIS_call,$peak_shift,$count,$Rltm_min_Rchx,$exon,$tr_stable_id,$e_rank,$e_start,$e_end,$e_start_phase,$e_end_phase,$e_stable_id,$seq,$tr_seq,$AA_seq,$exon_SNP,$CDS_tmp_length);
        
        #my @tr = ("440589_118628041");
        #foreach $transcript_start ( @tr){
        foreach $transcript_start ( keys %{$transcript_starts}){
            my %tr_SNPs;
            $tr_seq = '';
            $cnt++;
            
            # Split transcript and TIS
            $tr_id                          = $transcript_starts->{$transcript_start}{'transcript_id'};
            $TIS                            = $transcript_starts->{$transcript_start}{'start'};
            $strand                         = $transcript_starts->{$transcript_start}{'strand'};
            $start_codon                    = $transcript_starts->{$transcript_start}{'start_codon'};
            $dist_to_transcript_start       = $transcript_starts->{$transcript_start}{'dist_to_transcript_start'};
            $dist_to_aTIS                   = $transcript_starts->{$transcript_start}{'dist_to_aTIS'};
            $annotation                     = $transcript_starts->{$transcript_start}{'annotation'};
            $aTIS_call                      = $transcript_starts->{$transcript_start}{'aTIS_call'};
            $peak_shift                     = $transcript_starts->{$transcript_start}{'peak_shift'};
            $count                          = $transcript_starts->{$transcript_start}{'count'};
            $Rltm_min_Rchx                  = $transcript_starts->{$transcript_start}{'Rltm_min_Rchx'};
            
            
            #print "$tr_id,$TIS,$strand,$start_codon,$dist_to_transcript_start,$dist_to_aTIS,$annotation,$peak_shift,$count,$Rltm_min_Rchx\n";
            
            # Get exons of transcript and their chrom positions
            $exons = get_exons_for_tr($dbh_ENS,$tr_id);
            my $start_exon_check;
            $CDS_tmp_length = 0;
            # Hash will be needed to include positions of genomic positions of the ORF -> needed for ORF based counts
            my $tmp_orf_structure = {};
            my $orf_structure = {};
            my $orf_exon_counter = 0;
            my $stop_coord;
            
            # Get sequence starting from TIS and concatenated exon seqs
            # Loop over exons ascending for sense, descending for antisense
            # Reverse complement antisense
            # Set TIS=exon_start for sense, TIS=exon_end for antisense
            #if ($strand eq '1') {
            foreach $exon (sort {$a <=> $b} keys %{$exons}){
                
                $tr_stable_id   = $exons->{$exon}{'tr_stable_id'};
                $e_rank         = $exons->{$exon}{'rank'};
                $chr            = $exons->{$exon}{'chr'};
                $e_start        = $exons->{$exon}{'seq_region_start'};
                $e_end          = $exons->{$exon}{'seq_region_end'};
                $e_start_phase  = $exons->{$exon}{'phase'};
                $e_end_phase    = $exons->{$exon}{'end_phase'};
                $e_stable_id    = $exons->{$exon}{'e_stable_id'};
                
                #print "$tr_stable_id,$e_rank,$chr,$e_start,$e_end,$e_start_phase,$e_end_phase,$e_stable_id\n";
                
                # Find first exon and set start_pos to TIS site
                if ($TIS ~~ [$e_start..$e_end]) {
                    #print "$TIS\n";
                    $start_exon_check = 1;
                    
                    # Sense specific
                    if ($strand eq '1') { $e_start = $TIS; }
                    if ($strand eq '-1') { $e_end = $TIS; }
                    
                }
                
                # concatenation starts here...
                if ($start_exon_check) {
                    
                    #Slice exon SNPs from chrom SNP-hash and add transcript position to SNP_hash
                    my %exon_SNPs = fetch_SNPs_exon($e_start,$e_end,$SNP_chr);
                    
                    #print Dumper ('exon_SNPs_before',\%exon_SNPs);
                    
                    if (%exon_SNPs) {
                        foreach $exon_SNP (keys %exon_SNPs){
                            #print Dumper('between',$exon_SNPs{$exon_SNP});
                            
                            # Sense specific
                            #print "exon_end minus snp_pos plus CDS_tmp_length : $e_end minus $exon_SNPs{$exon_SNP}{'pos'} plus $CDS_tmp_length \n";
                            
                            # For the negative strand: take into account the length of the SNP entry (length($exon_SNPs{$exon_SNP}{'ref'}) -1)
                            $exon_SNPs{$exon_SNP}{'tr_pos'} =
                            ($strand eq '1') ? $exon_SNPs{$exon_SNP}{'pos'} - $e_start + 1 + $CDS_tmp_length
                            : $e_end - ($exon_SNPs{$exon_SNP}{'pos'} + (length($exon_SNPs{$exon_SNP}{'ref'}) -1)) + 1 + $CDS_tmp_length;
                            
                            $tr_SNPs{$exon_SNPs{$exon_SNP}{'tr_pos'}}{'pos'} =       $exon_SNPs{$exon_SNP}{'pos'};
                            $tr_SNPs{$exon_SNPs{$exon_SNP}{'tr_pos'}}{'tr_pos'} =    $exon_SNPs{$exon_SNP}{'tr_pos'};
                            $tr_SNPs{$exon_SNPs{$exon_SNP}{'tr_pos'}}{'ref'} =       ($strand eq '1') ? $exon_SNPs{$exon_SNP}{'ref'} : revdnacomp($exon_SNPs{$exon_SNP}{'ref'});
                            $tr_SNPs{$exon_SNPs{$exon_SNP}{'tr_pos'}}{'alt'} =       ($strand eq '1') ? $exon_SNPs{$exon_SNP}{'alt'} : revdnacomp($exon_SNPs{$exon_SNP}{'alt'});
                            $tr_SNPs{$exon_SNPs{$exon_SNP}{'tr_pos'}}{'chr'} =       $exon_SNPs{$exon_SNP}{'chr'};
                            $tr_SNPs{$exon_SNPs{$exon_SNP}{'tr_pos'}}{'af'} =        $exon_SNPs{$exon_SNP}{'af'};
                        }
                        
                        
                        #print Dumper ('exon_SNPs_after',\%exon_SNPs);
                    }
                    
                    #Sense specific
                    $seq = ($strand eq '1') ? get_sequence($chr,$e_start,$e_end) : revdnacomp(get_sequence($chr,$e_start,$e_end));
                    #print "$seq\n";
                    $tr_seq = $tr_seq . $seq;
                    
                    #Save ORF structure
                    $orf_exon_counter++;
                    $tmp_orf_structure->{$orf_exon_counter}->{'start'} = $e_start;
                    $tmp_orf_structure->{$orf_exon_counter}->{'end'} = $e_end;
                    $tmp_orf_structure->{$orf_exon_counter}->{'e_stable_id'} = $e_stable_id;
                    $tmp_orf_structure->{$orf_exon_counter}->{'rank'} = $e_rank;
                    $tmp_orf_structure->{$orf_exon_counter}->{'length'} = $e_end - $e_start + 1;
                    
                    # Build the temporary CDS length
                    my $l = length($seq);
                    $CDS_tmp_length = $CDS_tmp_length + ($e_end - $e_start + 1);
                    #print "length exon_seq = $l\n";
                    #print "CDS_tmp_length = $CDS_tmp_length\n";
                }
            }
            #print Dumper ('tr_SNPs',\%tr_SNPs);
            
            #Only output transcripts where TIS is in exonic regions. The transcripts missing the exon that overspans the TIS are excluded
            if ($start_exon_check) {
                #print "DNA  $tr_seq\n";
                #translate
                my $tr_seqs_all = translate($tr_seq,\%tr_SNPs);
                #print Dumper ('seqs_all',$tr_seqs_all);
                
                #Loop over all entries in AoH_tr_seqs_all
                my $out_cnt = 0;
                foreach (@$tr_seqs_all) {
                    $out_cnt++;
                    $AA_seq = ($_->{'AAseq'}) ? $_->{'AAseq'} : '';
                    $tr_seq = $_->{'seq'};
                    
                    #Only output transcripts that actually terminate with a STOP-codon (transcript forms with missing 3' exons resulting in missing STOP-codon are excluded)
                    if ($AA_seq =~ /\*$/ && $AA_seq ne '') {
                        
                        #Adapt the orf structure in function of the length of the AA seq to find the stop codon coordinate
                        if($strand eq '1'){
                            #Determine how long the translated sequence spans over the exon structure
                            my $orf_length = length($AA_seq) * 3;
                            my $i = 1;
                            for(;$orf_length>$tmp_orf_structure->{$i}->{'length'};$i++){
                                $orf_structure->{$i} = $tmp_orf_structure->{$i};
                                $orf_length = $orf_length - $tmp_orf_structure->{$i}->{'length'};
                            }
                            $orf_structure->{$i} = $tmp_orf_structure->{$i};
                            #Determine stop coordinate and adapt in last exon of the orf
                            $stop_coord = $orf_structure->{$i}->{'start'} + $orf_length - 1;
                            $orf_structure->{$i}->{'end'} = $stop_coord;
                            $orf_structure->{$i}->{'length'} = $orf_length;
                        }
                        elsif($strand eq '-1'){ #For antisense
                            my $orf_length = length($AA_seq) * 3;
                            my $i = 1;
                            for(;$orf_length>$tmp_orf_structure->{$i}->{'length'};$i++){
                                $orf_structure->{$i} = $tmp_orf_structure->{$i};
                                $orf_length = $orf_length - $tmp_orf_structure->{$i}->{'length'};
                            }
                            $orf_structure->{$i} = $tmp_orf_structure->{$i};
                            #Determine stop coordinate and adapt in last exon of the orf
                            $stop_coord = $orf_structure->{$i}->{'end'} - $orf_length + 1;
                            $orf_structure->{$i}->{'start'} = $stop_coord;
                            $orf_structure->{$i}->{'length'} = $orf_length;
                        }
                        
                        #Make underscore seperated sequences of start and stop coordinates
                        my $starts_seq="";
                        my $ends_seq="";
                        if($strand eq '1'){
                            for(my $i=1;$i<=scalar(keys(%{$orf_structure}));$i++){
                                if($i==1){
                                    $starts_seq = $orf_structure->{$i}->{'start'};
                                    $ends_seq = $orf_structure->{$i}->{'end'};
                                } elsif ($i>1) {
                                    $starts_seq = $starts_seq."_".$orf_structure->{$i}->{'start'};
                                    $ends_seq = $ends_seq."_".$orf_structure->{$i}->{'end'};
                                }
                            }
                        } elsif ($strand eq '-1'){
                            for(my $i=scalar(keys(%{$orf_structure}));$i>=1;$i--){
                                if($i==scalar(keys(%{$orf_structure}))){
                                    $starts_seq = $orf_structure->{$i}->{'start'};
                                    $ends_seq = $orf_structure->{$i}->{'end'};
                                } elsif($i<scalar(keys(%{$orf_structure})) && $i>=1){
                                    $starts_seq = $starts_seq."_".$orf_structure->{$i}->{'start'};
                                    $ends_seq = $ends_seq."_".$orf_structure->{$i}->{'end'};
                                }
                            }
                        }
                        
                        #Replace near-cognate start to cognate methionine...
                        $AA_seq = (substr($AA_seq,0,1) ne 'M') ? 'M'.substr($AA_seq,1) : $AA_seq;
                        print TMP_db $tr_stable_id.",".$chr.",".$strand.",".$TIS.",".$start_codon.",".$stop_coord.",".$starts_seq.",".$ends_seq.",".$dist_to_transcript_start.",".$dist_to_aTIS.",".$annotation.",".$aTIS_call.",".$peak_shift.",".$count.",".$Rltm_min_Rchx.",,,".$_->{'SNP_NS'}.",".$tr_seq.",".$AA_seq."\n";
                        
                        
                    }
                    else {
                        #print "$tr_stable_id,$AA_seq\n";
                    }
                }
            }
        }
        
        close(TMP_db);
        
        #print $chr." has ".$cnt." transcript_starts \n";
        ### Finish childs
        print "     * Finished translating chromosome ".$chr."\n";
        $dbh_ENS->disconnect();
        $dbh_results->disconnect();
        $pm->finish;
    }
    
    print "Waiting for all childs to finish...\n";
    $pm->wait_all_children();
    
    #Store in db
    store_in_db($dsn_results,$us_results,$pw_results,$work_dir,$chrs,$analysis_id,$snp);
    
    
}

### Store in DB ##
sub store_in_db{
    
    # Catch
    my $dsn         =   $_[0];
    my $us          =   $_[1];
    my $pw          =   $_[2];
    my $work_dir    =   $_[3];
    my $chrs        =   $_[4];
    my $analysis_id =   $_[5];
    my $snp         =   $_[6];
    
    # Init
    my $dbh     =   dbh($dsn,$us,$pw);
    
    my $snp_table_annot = ($snp eq "NO") ? "" : $snp."_";
    my $table   =   "TIS_".$analysis_id."_".$snp_table_annot."transcripts";
    
    my $query_drop = "DROP TABLE IF EXISTS `".$table."`";
    $dbh->do($query_drop);
    
    # Create table
    my $query = "CREATE TABLE IF NOT EXISTS `".$table."` (
    `tr_stable_id` varchar(128) NOT NULL default '',
    `chr` char(50) NOT NULL default '',
    `strand` int(2) NOT NULL default '',
    `start` int(10) NOT NULL default '',
    `start_codon` varchar(128) NOT NULL default '',
    `stop` int(10) NOT NULL default '',
    `starts_list` varchar(512) NOT NULL default '',
    `ends_list` varchar(512) NOT NULL default '',
    `dist_to_transcript_start` int(10) NOT NULL default '',
    `dist_to_aTIS` int(10) NOT NULL default 'NA',
    `annotation` varchar(128) NOT NULL default 'NA',
    `aTIS_call` varchar(128) NOT NULL default 'NA',
    `peak_shift` int(2) NOT NULL default '',
    `count` float default NULL,
    `Rltm_min_Rchx` decimal(11,8) NOT NULL default '0',
    `coverage` decimal(11,8) NOT NULL default '0',
    `FPKM` decimal(11,8) NOT NULL default '0',
    `SNP` varchar(256) NOT NULL default '0',
    `tr_seq` TEXT NOT NULL default '',
    `aa_seq` TEXT NOT NULL default '' )"  ;
    
    $dbh->do($query);
    
    # Store
    foreach my $chr (sort keys %{$chrs}){
        system("sqlite3 -separator , ".$resultDB." \".import ".$TMP."/".$chr."_tmp.csv ".$table."\"")== 0 or die "system failed: $?";
    }
    #Disconnect dbh
    $dbh->disconnect();
    
    # Unlink tmp csv files
    foreach my $chr (sort keys %{$chrs}){
        unlink $TMP."/".$chr."_tmp.csv";
    }
}

### Get the analysis_id that corresponds to the input paramaters for the TIS calling
sub get_analysis_id {
    
    # Catch
    my $dbh = $_[0];
    my $local_max = $_[1];
    my $min_count_aTIS = $_[2];
    my $Rltm_minus_Rchx_aTIS = $_[3];
    my $min_count_5UTR = $_[4];
    my $Rltm_minus_Rchx_5UTR = $_[5];
    my $min_count_3UTR = $_[6];
    my $Rltm_minus_Rchx_3UTR = $_[7];
    my $min_count_CDS = $_[8];
    my $Rltm_minus_Rchx_CDS = $_[9];
    my $min_count_no_translation = $_[10];
    my $Rltm_minus_Rchx_no_translation = $_[11];
    
    # Init
    my $analysis_id;
    
    #ID|CHX_lane|LTM_lane|ensembl_version|intergenic|species|local_max|min_count_aTIS|R_aTis|min_count_5UTR|R_5UTR|min_count_CDS|R_CDS|min_count_3UTR|R_3UTR|min_count_no_trans|R_no_trans
    #1|lane3|lane4|72|N|Mus_musculus|1|5|0.05|10|0.05|15|0.15|10|0.05|10|0.05
    
    my $query_analysis_id = "select max(ID) as ID from TIS_overview where local_max = ".$local_max. " and min_count_aTIS = ".$min_count_aTIS." and R_aTis = ".$Rltm_minus_Rchx_aTIS.
    " and min_count_5UTR = ".$min_count_5UTR." and R_5UTR = ".$Rltm_minus_Rchx_5UTR.
    " and min_count_3UTR = ".$min_count_3UTR." and R_3UTR= ".$Rltm_minus_Rchx_3UTR.
    " and min_count_CDS = ".$min_count_CDS." and R_CDS = ".$Rltm_minus_Rchx_CDS.
    " and min_count_no_trans = ".$min_count_no_translation." and R_no_trans = ".$Rltm_minus_Rchx_no_translation;
    
    #print "$query_analysis_id\n";
    
    my $sth = $dbh->prepare($query_analysis_id);
    $sth->execute();
    $analysis_id = $sth->fetch()->[0];
    print "Analysis ID = $analysis_id\n";
    return($analysis_id);
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

### get snp information from SQLite into hash ###
sub get_SNPs_per_chromosome {
    
    # Catch
    my $dbh = $_[0];
    my $chr_name = $_[1];
    my $snp = $_[2];
    
    # Init
    my $SNPs = {};
    my $query;
    
   	# Get chrs with seq_region_id
    if ($snp eq "samtools") {
        $query = "select chr,pos,ref,alt,case when new='m' then 0.5 else af end as af from snp_".$snp." where chr = '".$chr_name."' and af <> 'INDEL' and new<>'m' ";
    }
    elsif ($snp eq "samtools_dbSNP") {
        my @name_split = split(/_/,$snp);
        $query = "select chr,pos,ref,alt,case when new='m' then 0.5 else af end as af from snp_".$name_split[0]." where chr = '".$chr_name."' and af <> 'INDEL' ";
    }
    
    
    
    #print "$query\n";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    $SNPs = $sth->fetchall_hashref('pos');
    
    # Return
    return($SNPs);
    
}

sub fetch_SNPs_exon {
    
    # Catch
    my $e_start = $_[0];
    my $e_end   = $_[1];
    my %SNPs    = %{$_[2]};
    
    # Init
    my (@slice_list,%exon_SNPs_all,%exon_SNPs_def,$key,$value);
    
    
    ####LOOK INTO THIS TO CODE BETTER, FASTER (PERFORMANCE) TODOTODOTODO
    # Take hash slice based on e_start/e_end
    @slice_list = ($e_start..$e_end);
    #print Dumper(\@slice_list);
    @exon_SNPs_all{@slice_list} = @SNPs{@slice_list};
    #print Dumper('all',\%exon_SNPs_all);
    foreach $key (keys %exon_SNPs_all) {
        # Only print defined keys
        if ($exon_SNPs_all{$key}) { $exon_SNPs_def{$key} =  $exon_SNPs_all{$key} }
    }
    
    
    return %exon_SNPs_def;
}

### get transcripts from TIS-calling table ###
sub get_transcripts_per_chromosome {
    
    # Catch
    my $dbh = $_[0];
    my $chr_name = $_[1];
    my $analysis_id = $_[2];
    
    # Init
    my $transcript_starts = {};
    
    #print "chr_name = $chr_name\n";
    my $query = "SELECT transcript_id||'_'||start as transcript_start,transcript_id,biotype,chr,strand,start,dist_to_transcript_start,dist_to_aTIS,annotation,aTIS_call,start_codon,peak_shift,count,Rltm_min_Rchx FROM TIS_".$analysis_id." WHERE chr = '".$chr_name."'";
    ##Test for annotation='aTIS' and aTIS_call='True'
    #my $query = "SELECT transcript_id||'_'||start as transcript_start,transcript_id,biotype,chr,strand,start,dist_to_transcript_start,dist_to_aTIS,annotation,aTIS_call,start_codon,peak_shift,count,Rltm_min_Rchx FROM TIS_".$analysis_id." WHERE chr = '".$chr_name."' (";
    #print "$query\n";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    $transcript_starts = $sth->fetchall_hashref('transcript_start');
    
    # Return
    return($transcript_starts);
}

### Get all exons of transcript (Ensembl) ###
sub get_exons_for_tr {
    
    # Catch
    my $dbh = $_[0];
    my $transcript_id = $_[1];
    
    # Init
    my $exons = {};
    
    # Get exons
    my $query = "select tr.transcript_id tr_id, tr.stable_id tr_stable_id,etr.rank, r.name chr,e.seq_region_start,e.seq_region_end,e.seq_region_strand,e.phase,e.end_phase,e.stable_id e_stable_id ".
    "from transcript tr ".
    "   inner join exon_transcript etr on tr.transcript_id = etr.transcript_id ".
    "   inner join exon e on e.exon_id = etr.exon_id ".
    "   inner join seq_region r on r.seq_region_id = tr.seq_region_id ".
    " where tr.transcript_id = '". $transcript_id ."'";
    #print "$query\n";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    $exons = $sth->fetchall_hashref('rank');
    
    #Return
    return($exons);
    
}

### translation subroutine ###
sub translate {
    
    
    
    #Catch
    my $seq = $_[0];
    my %tr_SNPs = %{$_[1]};
    my (@AoH_tr_seqs,$href_tr_seq,$pos,$alt,$seq_tmp,@split_ALT,@AoH_tr_seqs_extra);
    
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
    
    
    #Array of hashes to hold resulting tr_seqs based on indel/SNPs changes
    # @AoH = (
    #  { 'SNP'  => @(trpos_ref_alt_af),
    #    'offset' => nr of extra bases due to inserts,
    #    'seq' => seq
    #  } ,
    #  { 'SNP'  => @(trpos_ref_alt_af),
    #    'offset' => nr of extra bases due to inserts,
    #    'seq' => seq
    #  }
    # )
    
    
    #Push input sequence in hash of possibilities
    push @AoH_tr_seqs, { 'SNP' => "", 'offset' => 0, 'seq' => $seq };
    
    #print Dumper('initial AoH',\@AoH_tr_seqs);
    # Only do the SNP-looping if input parameter is set to $snp <> "NO"
    if ($snp ne "NO") {
        # Loop over tr_SNPs and add in information
        # Needs to be sorted because if insertions are present, the downstream tr_pos will be influenced
        if (%tr_SNPs) {
            #print Dumper ('tr_SNPs',\%tr_SNPs); exit;
            foreach my $tr_SNP (sort {$a<=>$b}  keys %tr_SNPs) {
                #print Dumper ('tr_SNP',$tr_SNP);
                # Allele frequence equals 1 (only alt need to be retained)
                if ($tr_SNPs{$tr_SNP}{'af'} > 0.99) {
                    
                    #Only one ALTernative
                    if ($tr_SNPs{$tr_SNP}{'alt'} !~ /,/) {
                        
                        #Update info for all sequences
                        for $href_tr_seq ( @AoH_tr_seqs ) {
                            #print Dumper('AF99',$href_tr_seq);
                            #Get sequence
                            $seq_tmp = $href_tr_seq->{'seq'};
                            #print "BEF  $seq_tmp\n";
                            #Replace SNP in sequence (take into account the created offset caused by previous INDEL)
                            #print "$tr_SNPs{$tr_SNP}{'tr_pos'} + $href_tr_seq->{'offset'}\n";
                            
                            my $pos = $tr_SNPs{$tr_SNP}{'tr_pos'} + $href_tr_seq->{'offset'}  - 1;
                            
                            #pos = $pos;
                            #exit;
                            #$seq_tmp =~ s/\G($tr_SNPs{$tr_SNP}{'ref'})/$tr_SNPs{$tr_SNP}{'alt'}/;
                            $seq_tmp = substr($seq_tmp,0,$pos).$tr_SNPs{$tr_SNP}{'alt'}.substr($seq_tmp,$pos+(length($tr_SNPs{$tr_SNP}{'ref'})));
                            #print "AFT  $seq_tmp";
                            #Save changed sequence
                            my $pos_and_offset = $tr_SNPs{$tr_SNP}{'tr_pos'} + $href_tr_seq->{'offset'};
                            $href_tr_seq->{'SNP'} = $href_tr_seq->{'SNP'} . ":" . $pos_and_offset."_".$tr_SNPs{$tr_SNP}{'ref'}."_".$tr_SNPs{$tr_SNP}{'alt'}."_".$tr_SNPs{$tr_SNP}{'af'};
                            $href_tr_seq->{'offset'} = $href_tr_seq->{'offset'} + (length($tr_SNPs{$tr_SNP}{'alt'}) - length($tr_SNPs{$tr_SNP}{'ref'}));
                            $href_tr_seq->{'seq'} = $seq_tmp;
                            
                        }
                    }
                    
                    #Multiple ALTernative
                    elsif ($tr_SNPs{$tr_SNP}{'alt'} =~ /,/) {
                        
                        #Split ALT snps
                        @split_ALT = split(/,/, $tr_SNPs{$tr_SNP}{'alt'});
                        
                        #Update info for all sequences
                        for $href_tr_seq ( @AoH_tr_seqs ) {
                            
                            my $ref_SNP     = $href_tr_seq->{'SNP'};
                            my $ref_offset  = $href_tr_seq->{'offset'};
                            my $ref_seq      = $href_tr_seq->{'seq'};
                            
                            my $alt_cnt = 0;
                            foreach (@split_ALT) {
                                
                                $alt_cnt++;
                                # Update first one
                                if ($alt_cnt == 1) {
                                    
                                    #Get sequence
                                    $seq_tmp = $ref_seq;
                                    
                                    
                                    #Replace SNP in sequence (take into account the created offset caused by previous INDEL)
                                    my $pos = ($tr_SNPs{$tr_SNP}{'tr_pos'} + $ref_offset)  -1;
                                    #$seq_tmp =~ s/\G($tr_SNPs{$tr_SNP}{'ref'})/$_/;
                                    $seq_tmp = substr($seq_tmp,0,$pos).$_.substr($seq_tmp,$pos+(length($tr_SNPs{$tr_SNP}{'ref'})));
                                    #Save changed sequence
                                    my $pos_and_offset = $tr_SNPs{$tr_SNP}{'tr_pos'} + $href_tr_seq->{'offset'};
                                    $href_tr_seq->{'SNP'} = $ref_SNP . ":" . $pos_and_offset."_".$tr_SNPs{$tr_SNP}{'ref'}."_".$_."_".$tr_SNPs{$tr_SNP}{'af'};
                                    $href_tr_seq->{'offset'} = $ref_offset + (length($_) - length($tr_SNPs{$tr_SNP}{'ref'}));
                                    $href_tr_seq->{'seq'} = $seq_tmp;
                                }
                                # Push following ones in extra AoH
                                else {
                                    #Get sequence
                                    $seq_tmp = $ref_seq;
                                    
                                    
                                    #Replace SNP in sequence (take into account the created offset caused by previous INDEL)
                                    my $pos = ($tr_SNPs{$tr_SNP}{'tr_pos'} + $ref_offset)  -1;
                                    $seq_tmp = substr($seq_tmp,0,$pos).$_.substr($seq_tmp,$pos+(length($tr_SNPs{$tr_SNP}{'ref'})));
                                    #$seq_tmp =~ s/\G($tr_SNPs{$tr_SNP}{'ref'})/$_/;
                                    
                                    #Push new ones in extra AoH
                                    my $pos_and_offset = $tr_SNPs{$tr_SNP}{'tr_pos'} + $href_tr_seq->{'offset'};
                                    my $extra_SNP = $ref_SNP . ":" . $pos_and_offset."_".$tr_SNPs{$tr_SNP}{'ref'}."_".$_."_".$tr_SNPs{$tr_SNP}{'af'};
                                    my $extra_offset = $ref_offset + (length($_) - length($tr_SNPs{$tr_SNP}{'ref'}));
                                    my $extra_seq = $seq_tmp;
                                    push (@AoH_tr_seqs_extra, { 'SNP' => "$extra_SNP", 'offset' => $extra_offset, 'seq' => "$extra_seq" });
                                }
                            }
                        }
                    }
                    
                    #Add extra to regular AoH and empty extra
                    push (@AoH_tr_seqs, @AoH_tr_seqs_extra);
                    undef @AoH_tr_seqs_extra;
                    
                }
                
                # Allele frequence equals .5 (both ref and alt need te be retained)
                elsif ($tr_SNPs{$tr_SNP}{'af'} > 0.49 && $tr_SNPs{$tr_SNP}{'af'} < 0.51) {
                    
                    #Only one ALTernative
                    if ($tr_SNPs{$tr_SNP}{'alt'} !~ /,/) {
                        
                        #Keep existing one and push new in extra AoH
                        for $href_tr_seq ( @AoH_tr_seqs ) {
                            
                            #Get sequence
                            $seq_tmp = $href_tr_seq->{'seq'};
                            
                            #Replace SNP in sequence (take into account the created offset caused by previous INDEL)
                            my $pos = $tr_SNPs{$tr_SNP}{'tr_pos'} + $href_tr_seq->{'offset'}  - 1;
                            
                            #$seq_tmp =~ s/\G($tr_SNPs{$tr_SNP}{'ref'})/$tr_SNPs{$tr_SNP}{'alt'}/;
                            $seq_tmp = substr($seq_tmp,0,$pos).$tr_SNPs{$tr_SNP}{'alt'}.substr($seq_tmp,$pos+(length($tr_SNPs{$tr_SNP}{'ref'})));
                            
                            #Push new one in extra AoH
                            my $pos_and_offset = $tr_SNPs{$tr_SNP}{'tr_pos'} + $href_tr_seq->{'offset'};
                            my $extra_SNP = $href_tr_seq->{'SNP'} . ":" . $pos_and_offset."_".$tr_SNPs{$tr_SNP}{'ref'}."_".$tr_SNPs{$tr_SNP}{'alt'}."_".$tr_SNPs{$tr_SNP}{'af'};
                            my $extra_offset = $href_tr_seq->{'offset'} + (length($tr_SNPs{$tr_SNP}{'alt'}) - length($tr_SNPs{$tr_SNP}{'ref'}));
                            my $extra_seq = $seq_tmp;
                            push (@AoH_tr_seqs_extra, { 'SNP' => "$extra_SNP", 'offset' => $extra_offset, 'seq' => "$extra_seq" });
                            #print Dumper ('extra',\@AoH_tr_seqs_extra);
                            #exit;
                        }
                    }
                    
                    #Multiple ALTernative
                    elsif ($tr_SNPs{$tr_SNP}{'alt'} =~ /,/) {
                        
                        #Split ALT snps
                        @split_ALT = split(/,/, $tr_SNPs{$tr_SNP}{'alt'});
                        
                        #Update info for all sequences
                        for $href_tr_seq ( @AoH_tr_seqs ) {
                            
                            my $ref_SNP     = $href_tr_seq->{'SNP'};
                            my $ref_offset  = $href_tr_seq->{'offset'};
                            my $ref_seq      = $href_tr_seq->{'seq'};
                            
                            my $alt_cnt = 0;
                            foreach (@split_ALT) {
                                #Push all in extra AoH and keep original one
                                
                                #Get sequence
                                $seq_tmp = $ref_seq;
                                
                                
                                #Replace SNP in sequence (take into account the created offset caused by previous INDEL)
                                my $pos = ($tr_SNPs{$tr_SNP}{'tr_pos'} + $ref_offset)  -1;
                                #$seq_tmp =~ s/\G($tr_SNPs{$tr_SNP}{'ref'})/$_/;
                                $seq_tmp = substr($seq_tmp,0,$pos).$_.substr($seq_tmp,$pos+(length($tr_SNPs{$tr_SNP}{'ref'})));
                                
                                #Push new ones in extra AoH
                                my $pos_and_offset = $tr_SNPs{$tr_SNP}{'tr_pos'} + $href_tr_seq->{'offset'};
                                my $extra_SNP = $ref_SNP . ":" . $pos_and_offset."_".$tr_SNPs{$tr_SNP}{'ref'}."_".$_."_".$tr_SNPs{$tr_SNP}{'af'};
                                my $extra_offset = $ref_offset + (length($_) - length($tr_SNPs{$tr_SNP}{'ref'}));
                                my $extra_seq = $seq_tmp;
                                push @AoH_tr_seqs_extra, { 'SNP' => "$extra_SNP", 'offset' => $extra_offset, 'seq' => "$extra_seq" };
                            }
                        }
                    }
                    
                    #Add extra to regular AoH and empty extra
                    push (@AoH_tr_seqs, @AoH_tr_seqs_extra);
                    undef @AoH_tr_seqs_extra;
                    #print Dumper('normal + extra',\@AoH_tr_seqs);
                    
                    
                }
            }
        }
    }
    #print Dumper('final AoH',\@AoH_tr_seqs);
    
    for $href_tr_seq ( @AoH_tr_seqs ) {
        my $AAseq='';
        my ($i,$triplet,$AA);
        my (@snp_pos,@snp_ref,@snp_alt,@snp_af);
        my $SNP_NS = '';
        $AA='';
        #print Dumper("sequence",$href_tr_seq);
        if ($href_tr_seq->{'SNP'} ne '') {
            #Get snp positions, to check whether (non-)synonymous mutation!!
            #Only non-synonymous snps are included into the description line of proteins (with their according tr_position)
            my $SNP =$href_tr_seq->{'SNP'};
            #print "$SNP\n";
            my @split_SNP = split(/:/,$SNP);
            #print "@split_SNP\n";
            foreach my $tmp_snp (@split_SNP) {
                my @tmp_snpinfo = split(/_/,$tmp_snp);
                push (@snp_pos,$tmp_snpinfo[0]); push (@snp_ref,$tmp_snpinfo[1]); push (@snp_alt,$tmp_snpinfo[2]); push (@snp_af,$tmp_snpinfo[3]);
            }
            shift(@snp_pos); shift(@snp_ref); shift(@snp_alt); shift(@snp_af);
            #print "pos = @snp_pos\n";
            #print "ref = @snp_ref\n";
            #print "alt = @snp_alt\n";
            #print "af = @snp_af\n";
        }
        
        for($i=0;$i<=length($href_tr_seq->{'seq'})-2;$i+=3){
            my $cds_start = $i+1;
            my $cds_end  = $i+3;
            $triplet = substr($href_tr_seq->{'seq'},$i,3);
            #Stop translation if trailing bases
            if (length($triplet) < 3) { last; }
            $AA=$AA1{uc($triplet)};
            $AAseq = $AAseq . $AA;
            if ($AA eq '*') { last; }
            
            #Check the snp-inclusive triplets
            my ($AAREF,$tripletREF);
            if ($href_tr_seq->{'SNP'} ne '' && $snp_pos[0] && $snp_pos[0]) {
                #print Dumper($href_tr_seq);
                #print "$snp_pos[0], $cds_start, $cds_end\n";
                if ($snp_pos[0] >= $cds_start && $snp_pos[0] <= $cds_end) {
                    my $p = ($snp_pos[0] == $cds_start) ? 0 : ($snp_pos[0] == $cds_end) ? 2 : 1;
                    #print "$p\n";
                    my $tripletALT= $triplet;
                    substr($triplet,$p,1,$snp_ref[0]);
                    my $tripletREF = $triplet;
                    $AAREF=$AA1{uc($triplet)};
                    if ($AAREF ne $AA) {
                        my $tmp_SNP_NS = $snp_pos[0]."_".$snp_ref[0]."_".$snp_alt[0]."_".$snp_af[0];
                        $SNP_NS = $SNP_NS . $tmp_SNP_NS . ":";
                    }
                    shift(@snp_pos); shift(@snp_ref); shift(@snp_alt); shift(@snp_af);
                    
                    #print "$i, start=$cds_start,stop=$cds_end,tripletALT=$tripletALT,AAALT=$AA,tripletREF=$tripletREF,AAREF=$AAREF\n";
                }
            }
        }
        $SNP_NS =~ s/:$//;
        $href_tr_seq->{'seq'} = substr($href_tr_seq->{'seq'},0,length($AAseq)*3);
        $href_tr_seq->{'AAseq'} = $AAseq;
        $href_tr_seq->{'SNP_NS'} = $SNP_NS;
    }
    
    #print Dumper('final AoH_with_AA',\@AoH_tr_seqs);
    return \@AoH_tr_seqs;
}

### get sequence from binary chromosomes ###
sub get_sequence {
    
    #Catch
    my $chr     = $_[0];
    my $e_start = $_[1];
    my $e_end   = $_[2];
    
    # Init
    my $seq;
    
    open(IN, "< ".$BIN_chrom_dir."/".$chr.".fa");
    binmode(IN);
    
    my $buffer; # Buffer to get binary data
    my $length = $e_end - $e_start + 1; # Length of string to read
    my $offset = $e_start-1; # Offset where to start reading
    
    seek(IN,$offset,0);
    read(IN,$buffer,$length);
    
    #print "buffer:  $buffer\n";
    
    close(IN);
    
    # Return
    return($buffer);
    
    
}

### get reverse complement sequence ###
sub revdnacomp {
    my $dna = shift;
    my $revcomp = reverse($dna);
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}

### get complement sequence ###
sub dnacomp {
    my $dna = shift;
    $dna =~ tr/ACGTacgt/TGCAtgca/;
    return $dna;
}

### make combination of all members of array ###
sub combine {
    
    my ($list, $n) = @_;
    die "Insufficient list members" if $n > @$list;
    
    return map [$_], @$list if $n == 1;
    
    my @comb;
    
    for (my $i = 0; $i+$n <= @$list; $i++) {
        my $val = $list->[$i];
        my @rest = @$list[$i+1..$#$list];
        push @comb, [$val, @$_] for combine(\@rest, $n-1);
    }
    
    return @comb;
}


### update TIS_OVERVIEW table: add SNP info
sub update_TIS_overview {
    
    my ($snp, $analysis_id) = @_;
    my $dbh = dbh($dsn_results,$us_results,$pw_results);
    
    my $update = "update TIS_OVERVIEW set SNP = '".$snp."' where ID = ".$analysis_id;
    $dbh->do($update);
}