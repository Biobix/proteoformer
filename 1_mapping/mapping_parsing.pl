#!/usr/bin/perl -w

#####################################
##	PROTEOFORMER: deep proteome coverage through ribosome profiling and MS integration
##
##	Copyright (C) 2014 G. Menschaert, J.Crappé, E. Ndah, A. Koch & S. Steyaert
##
##	This program is free software: you can redistribute it and/or modify
##	it under the terms of the GNU General Public License as published by
##	the Free Software Foundation, either version 3 of the License, or
##	(at your option) any later version.
##
##	This program is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##	GNU General Public License for more details.
##
##	You should have received a copy of the GNU General Public License
##	along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
## 	For more (contact) information visit http://www.biobix.be/PROTEOFORMER
#####################################


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
##############
# ./1_mapping_parsing.pl --out_sqlite SQLite/results.db --offset standard

#For GALAXY
#1_mapping.pl --out_sqlite SQLite/results.db

# get the command line arguments
my ($work_dir,$tmpfolder,$out_sqlite,$offset_option,$offset_file_untr,$offset_file_tr,$help);

GetOptions(
"tmp:s" =>\$tmpfolder,                  	# Folder where temporary files are stored,                          			optional  argument (default = $TMP or $CWD/tmp env setting)
"work_dir:s" =>\$work_dir,              	# Working directory ,                                               			optional  argument (default = $CWD env setting)
"out_sqlite:s" =>\$out_sqlite,          	# sqlite DB output file,                                             			optional  argument (default = SQLite/results.db)
"offset:s" =>\$offset_option,                      # offset option (standard, from_file, plastid)                           optional  argument (default = standard)
#                                                       Standard: if you use the standard offset lengths used by Ingolia et al. (2009)
#                                                       From_file: from txt file with rpf_length, offset format as in plastid
#                                                       Plastid: if you have run the mapping module with plastid before
"offset_file_untr:s" =>\$offset_file_untr,  # offset input file for untreated data                                         mandatory argument if offset argument is 'from_file'
"offset_file_tr:s" => \$offset_file_tr,     # offset input file for treated data                                           mandatory argument if offset argument is 'from_file' and readtype argument in arguments table is 'ribo'
"help" => \$help                            #Help text option
);

if ($help){
    print_help_text();
    exit;
}

###########################################################################
#Check all input variable and/or get default values and set extra variables
###########################################################################

my $CWD             = getcwd;
my $HOME            = $ENV{'HOME'};
my $TMP             = ($ENV{'TMP'}) ? $ENV{'TMP'} : ($tmpfolder) ? $tmpfolder : "$CWD/tmp"; # (1) get the TMP environment variable, (2) get the $tmpfolder variable, (3) get current_working_dir/tmp
print "The following tmpfolder is used                          : $TMP\n";

#Check if tmpfolder exists, if not create it...
if (!-d "$TMP") {
    system ("mkdir ". $TMP);
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
if ($offset_option) {
    if ($offset_option eq "standard" || $offset_option eq "from_file" || $offset_option eq "plastid" || $offset_option eq "cst_5prime" || $offset_option eq "cst_3prime") {
        print "Offset source                                            : $offset_option\n";
    } else {
        die "Offset argument needs to be \" standard\", \"from_file\", \"plastid\", \"cst_5prime\" or \"cst_3prime\"!";
    }
} else {
    $offset_option = "standard";
    print "Offset source                                            : $offset_option\n";
}


#Get all other necessary arguments out of SQLite DB
my $db_sqlite_results  = $out_sqlite;
# Sqlite results
my $dsn_sqlite_results = "DBI:SQLite:dbname=$db_sqlite_results";
my $us_sqlite_results  = "";
my $pw_sqlite_results  = "";
# Get arguments vars
my ($species,$ensemblversion,$IGENOMES_ROOT,$cores,$seqFileName1,$seqFileName2,$mapper,$unique,$rpf_split,$FirstRankMultiMap,$truseq,$readtype,$out_bg_s_untr,$out_bg_as_untr,$out_bg_s_tr,$out_bg_as_tr,$out_sam_untr,$out_sam_tr,$run_name,$maxmultimap,$min_l_count,$max_l_count,$cst_prime_offset,$min_cst_prime_offset,$max_cst_prime_offset) = get_ARG_vars($dsn_sqlite_results,$us_sqlite_results,$pw_sqlite_results,$offset_option);

# Get executables
my $sqlite_loc = "sqlite3";

if ($offset_option eq "from_file"){
    if ($offset_file_untr) {
        print "Offset input file untreated                              : $offset_file_untr\n";
    } else {
        die "Do not forget the offset untreated input file if offset argument is \"from_file\"!";
    }
} else {
    $offset_file_untr = "";
}
if ($offset_option eq "from_file" && $readtype eq 'ribo'){
    if ($offset_file_tr) {
        print "Offset input file treated                                 : $offset_file_tr\n";
    } else {
        die "Do not forget the offset treated input file if offset argument is \"from_file\" and readtype is \"ribo\"!";
    }
} else {
    $offset_file_tr = "";
}


#Comment on input variables from argument table
print "The following species is used                            : $species\n";
print "The following Ensembl version is used                    : $ensemblversion\n";
print "The following igenomes folder is used			 : $IGENOMES_ROOT\n";
print "Number of cores to use for analysis			 : $cores\n";

#Conversion for species terminology
my $spec = (uc($species) eq "MOUSE") ? "Mus_musculus" 
: (uc($species) eq "RAT") ? "Rattus_norvegicus" 
: (uc($species) eq "HORSE") ? "Equus_caballus"
: (uc($species) eq "ARCTIC_SQUIRREL") ? "Urocitellus_parryii" 
: (uc($species) eq "C.ELEGANS") ? "Caenorhabditis_elegans"
: (uc($species) eq "CNECNA3") ? "Cryptococcus_neoformans_var_grubii_h99_gca_000149245" 
: (uc($species) eq "SL1344") ? "SL1344" 
: (uc($species) eq "MYC_ABS_ATCC_19977") ? "mycobacterium_abscessus_atcc_19977" 
: (uc($species) eq "HUMAN") ? "Homo_sapiens" 
: (uc($species) eq "ARABIDOPSIS") ? "Arabidopsis_thaliana"
: (uc($species) eq "EARTHMOSS") ? "Physcomitrium_patens"
: (uc($species) eq "FRUITFLY") ? "Drosophila_melanogaster" 
: (uc($species) eq "YEAST") ? "Saccharomyces_cerevisiae" 
: (uc($species) eq "ZEBRAFISH") ? "Danio_rerio" : "";
my $spec_short = (uc($species) eq "MOUSE") ? "mmu" 
: (uc($species) eq "RAT") ? "rnor" 
: (uc($species) eq "HORSE") ? "eca" 
: (uc($species) eq "ARCTIC_SQUIRREL") ? "upa"
: (uc($species) eq "CNECNA3") ? "cnecna3" 
: (uc($species) eq "SL1344") ? "sl1344" 
:  (uc($species) eq "MYC_ABS_ATCC_19977") ? "MYC_ABS_ATCC_19977" 
:(uc($species) eq "HUMAN") ? "hsa" 
: (uc($species) eq "ARABIDOPSIS") ? "ath" 
: (uc($species) eq "EARTHMOSS") ? "ppa" 
: (uc($species) eq "FRUITFLY") ? "dme" 
: (uc($species) eq "YEAST") ? "sce" 
: (uc($species) eq "ZEBRAFISH") ? "dre" : "";
#Old mouse assembly = NCBIM37, new one is GRCm38. Old human assembly = GRCh37, the new one is GRCh38
my $assembly = (uc($species) eq "MOUSE" && $ensemblversion >= 103 ) ? "GRCm39"
: (uc($species) eq "MOUSE" && $ensemblversion >= 70 && $ensemblversion < 103 ) ? "GRCm38"
: (uc($species) eq "MOUSE" && $ensemblversion < 70 ) ? "NCBIM37"
: (uc($species) eq "RAT" && $ensemblversion >=80 ) ? "Rnor_6.0"
: (uc($species) eq "RAT" && $ensemblversion < 80) ? "Rnor_5.0"
: (uc($species) eq "HORSE" && $ensemblversion > 94) ? "EquCab3.0"
: (uc($species) eq "ARCTIC_SQUIRREL" && $ensemblversion > 95) ? "ASM342692v1"
: (uc($species) eq "C.ELEGANS") ? "WBcel235"
: (uc($species) eq "HUMAN" && $ensemblversion >= 76) ? "GRCh38"
: (uc($species) eq "HUMAN" && $ensemblversion < 76) ? "GRCh37"
: (uc($species) eq "ARABIDOPSIS") ? "TAIR10"
: (uc($species) eq "EARTHMOSS") ? "Phypa_V3"
: (uc($species) eq "SL1344") ? "ASM21085v2"
: (uc($species) eq "MYC_ABS_ATCC_19977") ? "ASM6918v1"
: (uc($species) eq "ZEBRAFISH") ? "GRCz10"
: (uc($species) eq "YEAST") ? "R64-1-1"
: (uc($species) eq "CNECNA3") ? "CNA3"
: (uc($species) eq "FRUITFLY" && $ensemblversion < 79) ? "BDGP5"
: (uc($species) eq "FRUITFLY" && $ensemblversion >= 79) ? "BDGP6" : "";

my $chromosome_sizes = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Annotation/Genes/ChromInfo.txt";
print "chromosome_sizes_file=".$IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Annotation/Genes/ChromInfo.txt\n";


####################
## MAPPING PARSING #
####################
# Dependent on RIBO-seq or RNA-seq run one has to iterate 2 times (RIBO-seq, both untreated and treated, or just once (RNA-seq).
# For paired-end reads the read files are passed as comma-separated list.
my @loopfastQ;
if (uc($readtype) eq 'RIBO') {
    @loopfastQ = ($seqFileName1,$seqFileName2);
}
else {
    if (uc($readtype) =~ m/PE/) {
        my $concatfastQ = $seqFileName1.",".$seqFileName2;
        @loopfastQ = ($concatfastQ);
    } else {
        @loopfastQ = ($seqFileName1);
    }
}

# Start loop, create numeric fastq filename to give along to subroutines, next to fastq file
my $cnt=0;

foreach (@loopfastQ) {
    $cnt++;
    #next unless ($cnt == 2);
    my $fastqName = "fastq".$cnt;
    print "$fastqName, $_\n";

    if (uc($mapper) eq "BOWTIE") {
        die "\nNo mapping parsing for Bowtie available.\n\n";
    }
    elsif (uc($mapper) eq "BOWTIE2") {
        die "\nNo mapping parsing for Bowtie2 available.\n\n";
    }
    elsif (uc($mapper) eq "TOPHAT2") {

        print "Mapping parsing of TopHat2\n";
        my $start = time;
        if (uc($readtype) eq "RIBO" || $readtype eq "ribo_untr") {
            RIBO_parse_store($_,$fastqName, 'Y', $rpf_split,$offset_option,$offset_file_untr,$offset_file_tr,$out_bg_s_untr,$out_bg_as_untr,$out_bg_s_tr,$out_bg_as_tr,$out_sam_untr,$out_sam_tr,$run_name,$maxmultimap,$min_l_count,$max_l_count,$species,$cst_prime_offset,$min_cst_prime_offset,$max_cst_prime_offset); # Only A-site parsing if RIBO-seq
            if ($unique eq "N" && $FirstRankMultiMap eq "N") {RIBO_parse_store($_,$fastqName, $unique, $rpf_split,$offset_option,$offset_file_untr,$offset_file_tr,$out_bg_s_untr,$out_bg_as_untr,$out_bg_s_tr,$out_bg_as_tr,$out_sam_untr,$out_sam_tr,$run_name,$maxmultimap,$min_l_count,$max_l_count,$species,$cst_prime_offset,$min_cst_prime_offset,$max_cst_prime_offset)}
        }
        if (uc($readtype) eq "SE_POLYA") {
            RNA_parse_store($_,$fastqName, $unique, $truseq,$out_bg_s_untr,$out_bg_as_untr,$out_bg_s_tr,$out_bg_as_tr,$out_sam_untr,$out_sam_tr,$run_name,$maxmultimap);
        }
        if (uc($readtype) =~ m/PE/) {
            RNA_parse_store($_,$fastqName, $unique, $truseq,$out_bg_s_untr,$out_bg_as_untr,$out_bg_s_tr,$out_bg_as_tr,$out_sam_untr,$out_sam_tr,$run_name,$maxmultimap);
        }

        my $end = time - $start;
        printf("runtime TopHat mapping parsing: %02d:%02d:%02d\n\n",int($end/3600), int(($end % 3600)/60), int($end % 60));
    }
    elsif (uc($mapper) eq "STAR") {
        
		print "Mapping parsing of STAR\n";
        my $start = time;
        if (uc($readtype) eq "RIBO" || $readtype eq "ribo_untr") {
            RIBO_parse_store($_,$fastqName, 'Y', $rpf_split,$offset_option,$offset_file_untr,$offset_file_tr,$out_bg_s_untr,$out_bg_as_untr,$out_bg_s_tr,$out_bg_as_tr,$out_sam_untr,$out_sam_tr,$run_name,$maxmultimap,$min_l_count,$max_l_count,$species,$cst_prime_offset,$min_cst_prime_offset,$max_cst_prime_offset); # Only A-site parsing if RIBO-seq
			if ($unique eq "N" && $FirstRankMultiMap eq "N") {RIBO_parse_store($_,$fastqName, $unique, $rpf_split,$offset_option,$offset_file_untr,$offset_file_tr,$out_bg_s_untr,$out_bg_as_untr,$out_bg_s_tr,$out_bg_as_tr,$out_sam_untr,$out_sam_tr,$run_name,$maxmultimap,$min_l_count,$max_l_count,$species,$cst_prime_offset,$min_cst_prime_offset,$max_cst_prime_offset)}
        }
        if (uc($readtype) eq "SE_POLYA") {
			RNA_parse_store($_,$fastqName, $unique, $truseq,$out_bg_s_untr,$out_bg_as_untr,$out_bg_s_tr,$out_bg_as_tr,$out_sam_untr,$out_sam_tr,$run_name,$maxmultimap);
        }
        if (uc($readtype) =~ m/PE/) {
            RNA_parse_store($_,$fastqName , $unique, $truseq,$out_bg_s_untr,$out_bg_as_untr,$out_bg_s_tr,$out_bg_as_tr,$out_sam_untr,$out_sam_tr,$run_name,$maxmultimap);
        }

        my $end = time - $start;
        printf("runtime STAR mapping parsing: %02d:%02d:%02d\n\n",int($end/3600), int(($end % 3600)/60), int($end % 60));
    }

}



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

### SYSTEMERROR ###
sub systemError {
    my ($command,$returnValue,$errorMessage) = @_;
    if ($returnValue == -1){
        die "$command failed!\n$errorMessage\n\n";
    }
}


sub RIBO_parse_store {

    # Catch
    my $seqFile = $_[0];
    my $seqFileName = $_[1];
	my $uniq = $_[2];
    my $rpf_split = $_[3];
    my $offset_option = $_[4];
    my $offset_file_untr = $_[5];
    my $offset_file_tr = $_[6];
    my $out_bg_s_untr = $_[7];
    my $out_bg_as_untr = $_[8];
    my $out_bg_s_tr = $_[9];
    my $out_bg_as_tr = $_[10];
    my $out_sam_untr = $_[11];
    my $out_sam_tr = $_[12];
    my $run_name = $_[13];
    my $maxmultimap = $_[14];
    my $min_l_count = $_[15];
    my $max_l_count = $_[16];
    my $species = $_[17];
    my $cst_prime_offset = $_[18];
    my $min_cst_prime_offset = $_[19];
    my $max_cst_prime_offset = $_[20];

    my $bedgr_s = ($seqFileName  eq 'fastq1') ? $out_bg_s_untr : $out_bg_s_tr;
    my $bedgr_as = ($seqFileName  eq 'fastq1') ? $out_bg_as_untr : $out_bg_as_tr;
    my $sam = ($seqFileName  eq 'fastq1') ? $out_sam_untr : $out_sam_tr;
    my $offset_file = ($seqFileName eq 'fastq1') ? $offset_file_untr : $offset_file_tr;

    #If you only want the unique reads
    if ($unique eq 'Y' || ($unique eq 'N' && $FirstRankMultiMap eq 'Y')) {
        $seqFileName = $seqFileName;
    #For multimapping, make two sorts of count tables
    } elsif ($unique eq 'N' && $FirstRankMultiMap eq 'N') {
        if ($uniq eq 'N') {
            #the count tables with unique and multimapping reads
            $seqFileName = $seqFileName;
        } elsif ($uniq eq 'Y') {
            #the count tables with only unique reads, although the option for multimapping was selected
            $seqFileName = $seqFileName."_unique";
        }
    }

    ## Get chromosome sizes and cDNA identifiers #############
    print "Getting chromosome sizes and cDNA to chromosome mappings ...\n";
    my %chr_sizes = %{get_chr_sizes($chromosome_sizes)};

    print "Splitting genomic mapping per chromosome...\n";
    split_SAM_per_chr(\%chr_sizes,$work_dir,$seqFileName,$run_name,$sam,$uniq,$maxmultimap);

    #Create count tables

    #Init dbh
    my $dbh_count = dbh($dsn_sqlite_results,$us_sqlite_results,$pw_sqlite_results);

    # Create table if not exist
    my $query_table = "CREATE TABLE IF NOT EXISTS `count_".$seqFileName."` (
    `chr` char(50) NOT NULL default '',
    `strand` char(1) NOT NULL default '0',
    `start` int(10) NOT NULL default '0',
    `count` float default NULL)";

    # Create table if not exist
    my $query_table2 = "CREATE TABLE IF NOT EXISTS `count_".$seqFileName."_splitRPF` (
    `chr` char(50) NOT NULL default '',
    `strand` char(1) NOT NULL default '0',
    `start` int(10) NOT NULL default '0',
    `RPF` int(10) NOT NULL default '0',
    `count` float default NULL)";

    $dbh_count->do($query_table);
    $dbh_count->do($query_table2);

    # Disco
    $dbh_count->disconnect;
    
    # Construct p offset hash
    my $offset_hash = {};
    if($offset_option eq "plastid"){
        #Offset from plastid table
        #Init
        $offset_hash->{"min"} = 1000;
        $offset_hash->{"max"} = 0;
        my $offset_table;
        if ($seqFileName eq 'fastq1' || $seqFileName eq "fastq1_unique"){
            $offset_table = "p_offsets_untreated";
        } else {
            $offset_table = "p_offsets_treated";
        }
        #
        #Connect to result db
        my $dbh_plastid = dbh($dsn_sqlite_results, $us_sqlite_results, $pw_sqlite_results);
        
        #Query and parse plastid table
        my $query = "SELECT * FROM ".$offset_table.";";
        my $sth = $dbh_plastid->prepare($query);
        $sth->execute();
        while(my @row = $sth->fetchrow_array()){
            $offset_hash->{$row[0]} = $row[1];
            if($row[0]<$offset_hash->{"min"}){
                $offset_hash->{"min"} = $row[0];
            }
            if($row[0]>$offset_hash->{"max"}){
                $offset_hash->{"max"} = $row[0];
            }
        }
        
        $dbh_plastid->disconnect();
    } elsif($offset_option eq "from_file"){
        #Init
        $offset_hash->{"min"} = 1000;
        $offset_hash->{"max"} = 0;
        #Read in file
        open(my $FR, $offset_file) or die "Could not open $offset_file";
        
        #Parse
        while(my $line = <$FR>){
            if($line =~ /^(\d+)\s+(\d+)$/){
                my $length = $1;
                my $offset = $2;
                $offset_hash->{$length} = $offset;
                if($length<$offset_hash->{"min"}){
                    $offset_hash->{"min"} = $length;
                }
                if($length>$offset_hash->{"max"}){
                    $offset_hash->{"max"} = $length;
                }
            }
        }
    } elsif($offset_option eq "cst_5prime"){
        #Define constant 5prime offsets
        $offset_hash->{"min"} = $min_cst_prime_offset;
        $offset_hash->{"max"} = $max_cst_prime_offset;
        for(my $rpf = $offset_hash->{"min"}; $rpf<=$offset_hash->{"max"}; $rpf++){
            $offset_hash->{$rpf} = $cst_prime_offset;
        }
    } elsif($offset_option eq "cst_3prime"){
        #Translate constant 3 prime offsets into 5 prime-based offsets
        $offset_hash->{"min"} = $min_cst_prime_offset;
        $offset_hash->{"max"} = $max_cst_prime_offset;
        for(my $rpf = $offset_hash->{"min"}; $rpf<=$offset_hash->{"max"}; $rpf++){
            $offset_hash->{$rpf} = $rpf - $cst_prime_offset - 1;
        }
    } else {
        #Standard offset options from Ingolia paper
        if(uc($species) eq 'FRUITFLY'){
            $offset_hash->{25} = 12;
        }
        $offset_hash->{26} = 12;
        $offset_hash->{27} = 12;
        $offset_hash->{28} = 12;
        $offset_hash->{29} = 12;
        $offset_hash->{30} = 12;
        $offset_hash->{31} = 13;
        $offset_hash->{32} = 13;
        $offset_hash->{33} = 13;
        $offset_hash->{34} = 14;
        
        #Boundaries
        if(uc($species) eq 'FRUITFLY'){
            $offset_hash->{"min"} = 25;
        } else {
            $offset_hash->{"min"} = 26;
        }
        $offset_hash->{"max"} = 34;
    }
    
    #print Dumper $offset_hash;

    # Init multi core
    my $pm = new Parallel::ForkManager($cores);
    print "   Using ".$cores." core(s)\n   ---------------\n";

    foreach my $chr (keys %chr_sizes){

        ### Start parallel process
        $pm->start and next;

        ### DBH per process
        my $dbh = dbh($dsn_sqlite_results,$us_sqlite_results,$pw_sqlite_results);

        ### RIBO parsing
        my ($hits,$hits_splitRPF) = RIBO_parsing_genomic_per_chr($work_dir,$seqFileName,$run_name,$sam,$chr,$offset_hash,$min_l_count,$max_l_count);

        ### To File
        store_in_file_per_chr($hits,$hits_splitRPF,$dbh,$seqFileName,$chr,$run_name);

        ### Finish
        print "* Finished chromosome ".$chr."\n";
        $dbh->disconnect();
        $pm->finish;
    }

    # Finish all subprocesses
    $pm->wait_all_children;

    # Create indexes on riboseq tables


    my $table_name = "count_".$seqFileName;
    my $table_name_RPF = "count_".$seqFileName."_splitRPF";

    my $index1_st =  "create index if not exists ".$table_name."_chr on ".$table_name." (chr)";
    my $index2_st = "create index if not exists ".$table_name."_strand on ".$table_name." (strand)";

    my $system_cmd1 = $sqlite_loc." ".$db_sqlite_results." \"".$index1_st."\"";
    my $system_cmd2 = $sqlite_loc." ".$db_sqlite_results." \"".$index2_st."\"";

    system($system_cmd1);
    system($system_cmd2);


    ###################################
    #Start combining/generating output#
    ###################################

    ###SQLite DUMP

    # Gather all temp_chrom_csv files and dump import into SQLite DB
    my $temp_csv_all = $TMP."/genomic/".$run_name."_".$seqFileName.".csv";
    system("touch ".$temp_csv_all);

    foreach my $chr (keys %chr_sizes){
        my $temp_csv = $TMP."/genomic/".$run_name."_".$seqFileName."_".$chr."_tmp.csv";
        system ("cat ".$temp_csv." >>". $temp_csv_all);
    }
    #Remove _tmp.csv files
    system("rm -rf ".$TMP."/genomic/".$run_name."_".$seqFileName."_*_tmp.csv");

    system($sqlite_loc." -separator , ".$db_sqlite_results." \".import ".$temp_csv_all." ".$table_name."\"")== 0 or die "system failed: $?";
    system ("rm -rf ".$temp_csv_all);

    ####SPLIT FOR RPF
    # Gather all temp_chrom_csv_splitRPF files and dump import into SQLite DB
    my $temp_csv_all_RPF = $TMP."/genomic/".$run_name."_".$seqFileName."_splitRPF.csv";
    system("touch ".$temp_csv_all_RPF);

    foreach my $chr (keys %chr_sizes){
        my $temp_csv_splitRPF = $TMP."/genomic/".$run_name."_".$seqFileName."_".$chr."_tmp_splitRPF.csv";
        system ("cat ".$temp_csv_splitRPF." >>". $temp_csv_all_RPF);
    }
    #Remove _tmp.csv files
    system("rm -rf ".$TMP."/genomic/".$run_name."_".$seqFileName."_*_tmp_splitRPF.csv");

    system($sqlite_loc." -separator , ".$db_sqlite_results." \".import ".$temp_csv_all_RPF." ".$table_name_RPF."\"")== 0 or die "system failed: $?";
    system ("rm -rf ".$temp_csv_all_RPF);

    #Store mean pruned alignment length in arguments
    my $total_pruned_alignment_length = 0;
    my $total_nr_reads = 0;
    my ($mean_pruned_length,$query,$sth,@R);
    my $dbh_sqlite_results = dbh($dsn_sqlite_results,$us_sqlite_results,$pw_sqlite_results);

    if(uc($species) eq 'FRUITFLY'){

        foreach my $chr (keys %chr_sizes){

            #Get data from R_temp file
            open (RTMPin,"<".$TMP."/R_TMP_".$seqFileName."_".$chr.".txt") || die "ERROR reading R_tmp file \n";
            my $RTMPIN = <RTMPin>;
            my @R=split("\t",$RTMPIN);
            close(RTMPin);
            $total_pruned_alignment_length += $R[0];
            $total_nr_reads += $R[1];

        }
        $mean_pruned_length = $total_pruned_alignment_length/$total_nr_reads;
    }else{
        $mean_pruned_length = 1
    }

    $query = "INSERT INTO arguments (variable,value) VALUES (\'mean_length_".$seqFileName."\',\'".$mean_pruned_length."\')";
    $dbh_sqlite_results->do($query);
    $dbh_sqlite_results->disconnect();

    ###BEDGRAPH/BED output RIBOseq

    # Gather all + give header
    # BEDGRAPH /split for sense and antisense (since double entries, both (anti)sense cannot be visualized)
    my $bed_allgr_sense = $TMP."/genomic/".$run_name."_".$seqFileName."_sense.bedgraph";
    my $bed_allgr_antisense = $TMP."/genomic/".$run_name."_".$seqFileName."_antisense.bedgraph";
    open (BEDALLGRS,">".$bed_allgr_sense) || die "Cannot open the BEDGRAPH sense output file";
    open (BEDALLGRAS,">".$bed_allgr_antisense) || die "Cannot open the BEDGRAPH antisense output file";
    print BEDALLGRS "track type=bedGraph name=\"".$run_name."_".$seqFileName."_s\" description=\"".$run_name."_".$seqFileName."_s\" visibility=full color=3,189,0 priority=20\n";
    print BEDALLGRAS "track type=bedGraph name=\"".$run_name."_".$seqFileName."_as\" description=\"".$run_name."_".$seqFileName."_as\" visibility=full color=239,61,14 priority=20\n";
    close(BEDALLGRS);
    close(BEDALLGRAS);

    # BED has simple header for H2G2 upload
    my $bed_all = $TMP."/genomic/".$run_name."_".$seqFileName.".bed";
    open (BEDALL,">".$bed_all) || die "Cannot open the BED output file";
    print BEDALL "track name=\"".$run_name."_".$seqFileName."\" description=\"".$run_name."_".$seqFileName."\"\n";
    close(BEDALL);

    #Write temp files into bundled file (both BED and BEDGRAPH)
    foreach my $chr (keys %chr_sizes){
        if ($chr eq "MT") { next; } # Skip mitochondrial mappings
        my $temp_bed = $TMP."/genomic/".$run_name."_".$seqFileName."_".$chr."_tmp.bed";
        system ("cat ".$temp_bed." >>". $bed_all);
        my $temp_bedgr_s = $TMP."/genomic/".$run_name."_".$seqFileName."_".$chr."_s_tmp.bedgraph";
        my $temp_bedgr_as = $TMP."/genomic/".$run_name."_".$seqFileName."_".$chr."_as_tmp.bedgraph";

        system ("cat ".$temp_bedgr_s." >>". $bed_allgr_sense);
        system ("cat ".$temp_bedgr_as." >>". $bed_allgr_antisense);
    }
    
    ####RPF LENGTH SPECIFIC BEDGRAPH OUTPUT RIBOseq
    if($rpf_split eq "Y"){
        for (my $rpf_length=26; $rpf_length<=34; $rpf_length++){
            my $bedgr_rpf_s = $work_dir."/output/rpf_specific_bedgraph/".$run_name."_".$seqFileName."_rpf".$rpf_length."_s.bedgraph";
            my $bedgr_rpf_as = $work_dir."/output/rpf_specific_bedgraph/".$run_name."_".$seqFileName."_rpf".$rpf_length."_as.bedgraph";
            open (BEDRPFGRS,">".$bedgr_rpf_s) || die "Cannot open the RPF specific BEDGRAPH sense output file for RPF length ".$rpf_length;
            open (BEDRPFGRAS,">".$bedgr_rpf_as) || die "Cannot open the RPF specific BEDGRAPH antisense output file for RPF length ".$rpf_length;
            #Write headers
            print BEDRPFGRS "track type=bedGraph name=\"".$run_name."_".$seqFileName."_".$rpf_length."_s\" description=\"".$run_name."_".$seqFileName."_".$rpf_length."_s\" visibility=full color=3,189,0 priority=20\n";
            print BEDRPFGRAS "track type=bedGraph name=\"".$run_name."_".$seqFileName."_".$rpf_length."_as\" description=\"".$run_name."_".$seqFileName."_".$rpf_length."_as\" visibility=full color=239,61,14 priority=20\n";
            
            #Fetch all counts for that specific rpf length (strand specific)
            my $dbh_sqlite_results = dbh($dsn_sqlite_results,$us_sqlite_results,$pw_sqlite_results);
            #Sense
            my $query_fetch_rpf_s = "SELECT chr||'_'||start, chr, start, count FROM ".$table_name_RPF." WHERE strand='1' AND RPF='".$rpf_length."';";
            my $sth_fetch_rpf_s= $dbh_sqlite_results->prepare($query_fetch_rpf_s);
            $sth_fetch_rpf_s->execute();
            my $RPF_counts_s = $sth_fetch_rpf_s->fetchall_hashref("chr||'_'||start");
            my @keys_RPF_counts_s = keys %{$RPF_counts_s};
            @keys_RPF_counts_s = sort {lc($a) cmp lc($b)} @keys_RPF_counts_s;
            foreach my $chr_start (@keys_RPF_counts_s){
                my $RPF_start_pos = $RPF_counts_s->{$chr_start}->{"start"};
                my $RPF_start_pos_0based = $RPF_start_pos - 1;
                if($RPF_counts_s->{$chr_start}->{"chr"} eq "MT"){
                    print BEDRPFGRS "chrM\t".$RPF_start_pos_0based."\t".$RPF_start_pos."\t".$RPF_counts_s->{$chr_start}->{"count"}."\n";
                } else {
                    print BEDRPFGRS "chr".$RPF_counts_s->{$chr_start}->{"chr"}."\t".$RPF_start_pos_0based."\t".$RPF_start_pos."\t".$RPF_counts_s->{$chr_start}->{"count"}."\n";
                }
            }
            #Antisense
            my $query_fetch_rpf_as = "SELECT chr||'_'||start, chr, start, count FROM ".$table_name_RPF." WHERE strand='-1' AND RPF='".$rpf_length."';";
            my $sth_fetch_rpf_as= $dbh_sqlite_results->prepare($query_fetch_rpf_as);
            $sth_fetch_rpf_as->execute();
            my $RPF_counts_as = $sth_fetch_rpf_as->fetchall_hashref("chr||'_'||start");
            my @keys_RPF_counts_as = keys %{$RPF_counts_as};
            @keys_RPF_counts_as = sort {lc($a) cmp lc($b)} @keys_RPF_counts_as;
            foreach my $chr_start (@keys_RPF_counts_as){
                my $RPF_start_pos = $RPF_counts_as->{$chr_start}->{"start"};
                my $RPF_start_pos_0based = $RPF_start_pos - 1;
                if($RPF_counts_as->{$chr_start}->{"chr"} eq "MT"){
                    print BEDRPFGRAS "chrM\t".$RPF_start_pos_0based."\t".$RPF_start_pos."\t".$RPF_counts_as->{$chr_start}->{"count"}."\n";
                } else {
                    print BEDRPFGRAS "chr".$RPF_counts_as->{$chr_start}->{"chr"}."\t".$RPF_start_pos_0based."\t".$RPF_start_pos."\t".$RPF_counts_as->{$chr_start}->{"count"}."\n";
                }
            }
            
            #Disconnect
            $dbh_sqlite_results->disconnect();
            
            close(BEDALLGRS);
            close(BEDALLGRAS);
        }
    }
    
    #Remove _tmp.bed/_tmp.bedgraph and sorted sam files and move bundled files into output folder
    #    system("mv ". $bed_allgr_sense." ".$work_dir."/output/".$run_name."_".$seqFileName."_s.bedgraph");
    #    system("mv ". $bed_allgr_antisense." ".$work_dir."/output/".$run_name."_".$seqFileName."_as.bedgraph");
    system("mv ". $bed_allgr_sense." ".$bedgr_s);
    system("mv ". $bed_allgr_antisense." ".$bedgr_as);
    system("mv ". $bed_all." ".$work_dir."/output/".$run_name."_".$seqFileName.".bed");
    system("rm -rf ".$TMP."/genomic/".$run_name."_".$seqFileName."_*_tmp.bed");
    system("rm -rf ".$TMP."/genomic/".$run_name."_".$seqFileName."_*_s_tmp.bedgraph");
    system("rm -rf ".$TMP."/genomic/".$run_name."_".$seqFileName."_*_as_tmp.bedgraph");
    system("rm -rf ".$TMP."/genomic/".$seqFileName."_*");
}

### Generate bedgraph and bam for SE RNA-seq ###
### Also works for PE RNA-seq (see: http://cancan.cshl.edu/labmembers/gordon/files/bedtools_genome_cov1.png ) ###
### BED file creation still needs to be included (for H2G2 browser environment) ###


sub RNA_parse_store {

    # Catch
    my $seqFile = $_[0];
    my $seqFileName = $_[1];
    my $uniq = $_[2];
    my $truseq = $_[3];
    my $out_bg_s_untr = $_[4];
    my $out_bg_as_untr = $_[5];
    my $out_bg_s_tr = $_[6];
    my $out_bg_as_tr = $_[7];
    my $out_sam_untr = $_[8];
    my $out_sam_tr = $_[9];
    my $run_name = $_[10];
    my $maxmultimap = $_[11];
    
    my $bedgr_s = ($seqFileName  eq 'fastq1') ? $out_bg_s_untr : $out_bg_s_tr;
    my $bedgr_as = ($seqFileName  eq 'fastq1') ? $out_bg_as_untr : $out_bg_as_tr;
    my $sam = ($seqFileName  eq 'fastq1') ? $out_sam_untr : $out_sam_tr;

    #my $sorted_bam = $work_dir."/".$mapper."/".$seqFileName."/Aligned.sorted.bam";
    my $sorted_bam = $work_dir."/".$mapper."/".$seqFileName."/Aligned.sortedByCoord.out.bam";

	#if ($unique eq 'Y') {
	#	$sorted_bam = $work_dir."/".$mapper."/".$seqFileName."/Aligned.sorted.unique.bam";
	#}

    ## Get chromosome sizes and cDNA identifiers
    print "Getting chromosome sizes and cDNA to chromosome mappings ... \n";
    my %chr_sizes = %{get_chr_sizes($chromosome_sizes)};

    print "Splitting genomic mapping per chromosome...\n";
    split_SAM_per_chr(\%chr_sizes, $work_dir, $seqFileName, $run_name, $sam, $uniq, $maxmultimap);

    #Create count tables

    #Init dbh
    my $dbh_count = dbh($dsn_sqlite_results, $us_sqlite_results, $pw_sqlite_results);


    #Create table if not exist
    my $query_table = "CREATE TABLE IF NOT EXISTS `count_".$seqFileName."` (
    `chr` char(50) NOT NULL default '',
    `strand` char(1) NOT NULL default '0',
    `start` int(10) NOT NULL default '0',
    `count` float default NULL)";


    $dbh_count->do($query_table);

    #Disconnect
    $dbh_count->disconnect;

    #Init multi core
    my $pm = new Parallel::ForkManager($cores);
    print "   Using ".$cores." core(s)\n  --------------\n";

    foreach my $chr (keys %chr_sizes){

        ### Start parallel process
        $pm->start and next;

        ### DBH per process
        my $dbh = dbh($dsn_sqlite_results, $us_sqlite_results, $pw_sqlite_results);

        ### RNA parsing
        my $hits = RNA_parsing_genomic_per_chr($work_dir,$seqFileName,$run_name,$sam,$chr);

        ### To file
        store_in_file_per_chr_RNA($hits,$dbh,$seqFileName,$chr,$run_name, $truseq);

        ### Finish
        print "* Finished chromosome ".$chr."\n";
        $dbh->disconnect();
        $pm->finish;
    }

    #Finish all subprocesses
    $pm->wait_all_children;

    # Creat indexes on count table
    my $table_name = "count_".$seqFileName;
    my $index1_st = "create index if not exists ".$table_name."_chr on ".$table_name." (chr)";
    my $index2_st = "create index if not exists ".$table_name."_strand on ".$table_name." (strand)";
    my $system_cmd1 = $sqlite_loc." ".$db_sqlite_results." \"".$index1_st."\"";
    my $system_cmd2 = $sqlite_loc." ".$db_sqlite_results." \"".$index2_st."\"";
    system($system_cmd1);
    system($system_cmd2);

    ####### Combine and generate all output ########
    ##SQLite Dump

    #Gather all temp_chrom_csv files in one file
    my $temp_csv_all = $TMP."/genomic/".$run_name."_".$seqFileName."_".$seqFileName.".csv";
    system("touch ".$temp_csv_all);

    foreach my $chr (keys %chr_sizes){
        my $temp_csv = $TMP."/genomic/".$run_name."_".$seqFileName."_".$chr."_tmp.csv";
        system("cat ".$temp_csv." >> ".$temp_csv_all);
    }

    #Remove chr tmp csv files
    system("rm -rf ".$TMP."/genomic/".$run_name."_".$seqFileName."_*_tmp.csv");

    #Dump into sqlite
    system($sqlite_loc." -separator , ".$db_sqlite_results." \".import ".$temp_csv_all." ".$table_name."\"")==0 or die "system failed: $?";
    system("rm -rf ".$temp_csv_all);

    # BEDGRAPH /split for sense and antisense (since double entries, both (anti)sense cannot be visualized)
    my $bed_allgr_sense = $TMP."/genomic/".$run_name."_".$seqFileName."_sense.bedgraph";
    my $bed_allgr_antisense = $TMP."/genomic/".$run_name."_".$seqFileName."_antisense.bedgraph";
	system("touch ".$bedgr_s);
	system("touch ".$bedgr_as);

	# Convert BAM to BEDGRAPH files for each strand (strand assignment depends on RNAseq protocol)
    if($truseq eq "Y"){
        my $command_bedgraphS = "bedtools genomecov -bg -strand - -split -ibam ".$sorted_bam." -g ".$chromosome_sizes." >".$bedgr_s." 2>&1";
        my $command_bedgraphAS = "bedtools genomecov -bg -strand + -split -ibam ".$sorted_bam." -g ".$chromosome_sizes." >".$bedgr_as." 2>&1";
        system($command_bedgraphS);
        system($command_bedgraphAS);
    } else {
        my $command_bedgraphS = "bedtools genomecov -bg -strand + -split -ibam ".$sorted_bam." -g ".$chromosome_sizes." >".$bedgr_s." 2>&1";
        my $command_bedgraphAS = "bedtools genomecov -bg -strand - -split -ibam ".$sorted_bam." -g ".$chromosome_sizes." >".$bedgr_as." 2>&1";
        system($command_bedgraphS);
        system($command_bedgraphAS);
    }

    #Adapt BEDGRAPH files to UCSC conventions
    system("mv ".$bedgr_s." ".$TMP."/genomic/bedgraph_old_sense.bedgraph");
    system("mv ".$bedgr_as." ".$TMP."/genomic/bedgraph_old_antisense.bedgraph");
    open(BEDALLGRS,">".$bedgr_s) || die "Cannot open the BEDGRAPH sense output file";
    open(BEDALLGRAS,">".$bedgr_as) || die "Cannot open the BEDGRAPH antisense output file";
    print BEDALLGRS "track type=bedGraph name=\"".$run_name."_".$seqFileName."_s\" description=\"".$run_name."_".$seqFileName."_s\" visibility=full color=3,189,0 priority=20\n";
    print BEDALLGRAS "track type=bedGraph name=\"".$run_name."_".$seqFileName."_as\" description=\"".$run_name."_".$seqFileName."_as\" visibility=full color=239,61,14 priority=20\n";
    open(OLDBEDGRS,"<".$TMP."/genomic/bedgraph_old_sense.bedgraph");
    while(<OLDBEDGRS>){
        if($_ =~ m/^MT/){
            print BEDALLGRS "chrM";
        } else {
            print BEDALLGRS "chr".$_;
        }
    }
    close(OLDBEDGRS);
    system("rm -rf ".$TMP."/genomic/bedgraph_old_sense.bedgraph");
    close(BEDALLGRS);
    open(OLDBEDGRAS,"<".$TMP."/genomic/bedgraph_old_antisense.bedgraph");
    while(<OLDBEDGRAS>){
        if($_ =~ m/^MT/){
            print BEDALLGRAS "chrM";;
        } else {
            print BEDALLGRAS "chr".$_;
        }
    }
    close(OLDBEDGRAS);
    system("rm -rf ".$TMP."/genomic/bedgraph_old_antisense.bedgraph");
    close(BEDALLGRAS);
    
    system("rm -rf ".$TMP."/genomic/");
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

### SPLIT SAM PER CHR ###
sub split_SAM_per_chr {

    # Catch
    my %chr_sizes = %{$_[0]};
    my $work_dir = $_[1];
    my $seqFileName   = $_[2];
    my $run_name       = $_[3];
    my $sam = $_[4];
    my $uni = $_[5];
    my $maxmultimap = $_[6];

    my @splitsam = split(/\//, $sam );
    my $samFileName = $splitsam[$#splitsam];

    my $directory = $work_dir."/".$mapper."/".$seqFileName."/";
    my ($chr,@mapping_store,$file_in_loc,$file_in,$file_out);

    #Create chromosome sub directory in temp
    system("mkdir -p ".$TMP."/genomic/");
    system("rm -f ".$TMP."/genomic/".$samFileName."_*");   # Delete existing

    # Touch per chr
    foreach $chr (keys %chr_sizes){
        system("touch ".$TMP."/genomic/".$samFileName."_".$chr);   # Touch new
    }

    ## Split files into chromosomes
    $file_in_loc = $sam;
    system ("mv ". $file_in_loc ." ".$TMP."/genomic/");
    $file_in = $TMP."/genomic/".$samFileName;

    # Open
    open (I,"<".$file_in) || die "Cannot open ".$file_in." file\n";


    #For unsorted SAM file (genomic location)
    my $prev_chr="0";

    while(my $line=<I>){

        #Skip annotation lines
        if ($line =~ m/^@/) { next; }

        #Process alignment line
        @mapping_store = split(/\t/,$line);
        $chr = $mapping_store[2];

        # Unique vs. (Unique+Multiple) alignment selection
        # NH:i:1 means that only 1 alignment is present
        # HI:i:xx means that this is the xx-st ranked (for Tophat ranking starts with 0, for STAR ranking starts with 1)
        if ($uni eq "Y") {
            next unless (($mapping_store[4] == 255 && uc($mapper) eq "STAR") || ($line =~ m/NH:i:1\D/ && uc($mapper) eq "TOPHAT2"));
        }
        elsif ($uni eq "N") {
            #If multiple: best scoring or random (if equally scoring) is chosen
            if ($FirstRankMultiMap eq "Y") {
                next unless (($mapping_store[12] eq "HI:i:1" && uc($mapper) eq "STAR") || (($line =~ m/HI:i:0/ || $line =~ m/NH:i:1\D/) && uc($mapper) eq "TOPHAT2"));
            }
            #Keep all mappings, also MultipleMapping locations are available (alternative to pseudogenes mapping) GM:07-10-2013
            #Note that we only retain the up until <16 multiple locations (to avoid including TopHat2 peak @ 16)
            #For STAR the maxMultiMap is not included in the output (e.g. if maxmultimap is set to 16, up untill 15 is included)
            #For TopHat2 the maxMultiMap is included in the output (e.g. if maxmultimap is set to 16, up untill 16 is included)
            #In order to have the same actual limit, the maxmultimap is discarded (see record below)
            next unless ( $line !~ m/NH:i:$maxmultimap/ );

        }

        # Write off
        if ($prev_chr ne $chr) {
            if ($prev_chr ne "0") { close(A);}
            $file_out = $TMP."/genomic/".$samFileName;
            open (A,">>".$file_out."_".$chr) || die "Cannot open the sep file";
            print A $line;
        }
        elsif ($prev_chr eq $chr) {
            print A $line;
        }
        $prev_chr = $chr;

    }

    # Close
    close(A);
    close(I);

    # Move back from TMP folder to original location
    system ("mv ". $file_in ." ".$sam);
}


### RIBO PARSE PER CHR ###
sub RIBO_parsing_genomic_per_chr {

    #Catch
    my $work_dir = $_[0];
    my $seqFileName = $_[1];
    my $run_name = $_[2];
    my $sam = $_[3];
    my $chr = $_[4];
    my $offset_hash = $_[5];
    my $min_l_count = $_[6];
    my $max_l_count = $_[7];

    my @splitsam = split(/\//, $sam );
    my $samFileName = $splitsam[$#splitsam];


    #Initialize
    my $directory = $work_dir."/".$mapper."/".$seqFileName."/";
    my $hits_genomic = {};
    my $hits_genomic_splitRPF = {};
    my $plus_count = 0; my $min_count = 0; my $lineCount = 0;
    my ($genmatchL,$offset,$start,$intron_total,$extra_for_min_strand,$pruned_alignmentL,$prunedalignment);
    my $lendistribution;
    my $read_mapped_length = 0;
    my $nr_reads = 0;

    open (LD,">".$TMP."/LD_".$seqFileName."_".$chr.".txt");
    open (I,"<".$TMP."/genomic/".$samFileName."_".$chr) || die "Cannot open ".$samFileName." file\n";
    while(my $line=<I>){

        $lineCount++;

        #Process alignment line
        my @mapping_store = split(/\t/,$line);

        #Get strand specifics
        # Sam flag is bitwise. (0x10 SEQ being reverse complemented)
        # 0x10 = 16 in decimal. -> negative strand.
        my $strand = ($mapping_store[1] & 16) ? "-": "+";
        my $CIGAR = $mapping_store[5];

        #Parse CIGAR to obtain offset,genomic matching length and total covered intronic region before reaching the offset
        if(uc($species) eq 'FRUITFLY'){
            ($offset,$genmatchL,$intron_total,$extra_for_min_strand,$pruned_alignmentL,$prunedalignment) = parse_dme_RIBO_CIGAR($CIGAR,$strand);
            $lendistribution->{$genmatchL}++;


            if ($pruned_alignmentL > 0){
                $read_mapped_length = $read_mapped_length + $pruned_alignmentL;
                $nr_reads++;
                if ($strand eq "+") { $plus_count++;} elsif ($strand eq "-") { $min_count++; }
                foreach my $n (keys %{$prunedalignment}){
                    $start = ($strand eq "+") ? $mapping_store[3] + $prunedalignment->{$n}{'intron_total'} + $n -1: ($strand eq "-") ? $mapping_store[3] -$n - $prunedalignment->{$n}{'intron_total'} + $extra_for_min_strand : "";
                    if ( $genmatchL >= $min_l_count && $genmatchL <= $max_l_count) {
                        if ( exists $hits_genomic->{$chr}->{$start}->{$strand} ){
                            $hits_genomic->{$chr}->{$start}->{$strand} = $hits_genomic->{$chr}->{$start}->{$strand} + (1/$pruned_alignmentL);
                        }else {
                            $hits_genomic->{$chr}->{$start}->{$strand} = 0;
                            $hits_genomic->{$chr}->{$start}->{$strand} = $hits_genomic->{$chr}->{$start}->{$strand} + (1/$pruned_alignmentL);
                        }
                    }
                }
            }
        }else{
            ($offset,$genmatchL,$intron_total,$extra_for_min_strand) = parse_RIBO_CIGAR($CIGAR,$strand,$offset_hash);
            $lendistribution->{$genmatchL}++;
            #Determine genomic position based on CIGAR string output and mapping position and direction
            $start = ($strand eq "+") ? $mapping_store[3] + $offset + $intron_total : ($strand eq "-") ? $mapping_store[3] - $offset - $intron_total + $extra_for_min_strand -1 : "";

            $hits_genomic_splitRPF->{$chr}->{$start}->{$genmatchL}->{$strand}++;

            if ( $genmatchL >= $min_l_count && $genmatchL <= $max_l_count) {
                $hits_genomic->{$chr}->{$start}->{$strand}++;
                if ($strand eq "+") { $plus_count++;} elsif ($strand eq "-") { $min_count++; }
            }
        }
    }

    #Create R_temp file and write read alignment statistics to file
    open (RTMP,">".$TMP."/R_TMP_".$seqFileName."_".$chr.".txt") || die "ERROR opening R_tmp file to write \n";
    print RTMP "$read_mapped_length\t$nr_reads";

    my $cnttot = $plus_count + $min_count;
    for my $key ( sort { $a <=> $b } keys %$lendistribution ) {
        print LD "$key\t$lendistribution->{$key}\n";
    }

    close(LD);
    close(I);
    close(RTMP);
    return($hits_genomic,$hits_genomic_splitRPF);
}

### RNA PARSE PER CHR ###
sub RNA_parsing_genomic_per_chr {

    #Catch
    my $work_dir = $_[0];
    my $seqFileName = $_[1];
    my $run_name = $_[2];
    my $sam = $_[3];
    my $chr = $_[4];

    my @splitsam = split(/\//, $sam);
    my $samFileName = $splitsam[$#splitsam];

    #Initialize
    my $lineCount=0; my $plus_count = 0; my $min_count = 0;
    my $hits_genomic = {};
    my ($genmatchL,$offset,$intron_total,$extra_for_min_strand);
    my $start;
    my $lendistribution;
    my $read_mapped_length=0; #needed if you want to expand this RNA module to be expanded for fruitfly)
    my $nr_reads=0;

    open(LD,">".$TMP."/LD_".$seqFileName."_".$chr.".txt");
    open(I,"<".$TMP."/genomic/".$samFileName."_".$chr) || die "Cannot open ".$samFileName." file for chromosome ".$chr."\n";
    while(my $line=<I>){

        $lineCount++;

        #Proccess alignment line
        my @mapping_store = split(/\t/,$line);

        #Get strand specifics
        #The sam flag (second tab in alignment line) is bitwise. (0x10 SEQ being reverse complemented)
        # 0x10 = 16 in decimal. -> negative strand
        my $strand = ($mapping_store[1] & 16) ? "-": "+"; # & is a bitwise operation
        my $CIGAR = $mapping_store[5];

        #Parse CIGAR to obtain offset, genomic matching length and total covered intronic region before reaching the offset (Not adapted for FRUITFLY yet!). For RNA, offset will be set to 26 nt as in RiboTaper paper (Calviello et al. 2015).
        ($offset,$genmatchL,$intron_total,$extra_for_min_strand) = parse_RNA_CIGAR($CIGAR,$strand);
        $lendistribution->{$genmatchL}++;

        #Determine the genomic position based on CIGAR string output, mapping position and strand direction
        $start = ($strand eq "+") ? $mapping_store[3] + $offset + $intron_total : ($strand eq "-") ? $mapping_store[3] - $offset - $intron_total + $extra_for_min_strand + 1 : "";

        #Save in reads hash
        $hits_genomic->{$chr}->{$start}->{$strand}++;
        if ($strand eq "+") {$plus_count++;} elsif ($strand eq "-") { $min_count ++; }

    }

    #Create R_temp file and write read alignment statistics to file (for fruitfly, although not used for RNA so far)
    open (RTMP,">".$TMP."/R_TMP_".$seqFileName."_".$chr.".txt") || die "ERROR opening R_tmp file to write \n";
    print RTMP "$read_mapped_length\t$nr_reads";

    my $cnttot = $plus_count + $min_count;
    for my $key ( sort { $a <=> $b } keys %$lendistribution ) {
        print LD "$key\t$lendistribution->{$key}\n";
    }

    close(LD);
    close(I);
    close(RTMP);

    return $hits_genomic;
}

### STORE IN FILE PER CHR ###
sub store_in_file_per_chr {

    # Catch
    my $hits = $_[0];
    my $hits_splitRPF = $_[1];
    my $dbh  = $_[2];
    my $seqFileName = $_[3];
    my $chromosome = $_[4];
    my $run_name = $_[5];

    my $directory = $work_dir."/".$mapper."/".$seqFileName."/";

    #Init temporary csv-file/bed-file/bedgraph-file
    my $temp_csv = $TMP."/genomic/".$run_name."_".$seqFileName."_".$chromosome."_tmp.csv";
    my $temp_csv_splitRPF = $TMP."/genomic/".$run_name."_".$seqFileName."_".$chromosome."_tmp_splitRPF.csv";
    my $temp_bed = $TMP."/genomic/".$run_name."_".$seqFileName."_".$chromosome."_tmp.bed";
    my $temp_bedgr_s = $TMP."/genomic/".$run_name."_".$seqFileName."_".$chromosome."_s_tmp.bedgraph";
    my $temp_bedgr_as = $TMP."/genomic/".$run_name."_".$seqFileName."_".$chromosome."_as_tmp.bedgraph";

    open TMP, "+>>".$temp_csv or die $!;
    open TMP_SPLITRPF, "+>>".$temp_csv_splitRPF or die $!;
    open TMPBED, "+>>".$temp_bed or die $!;
    open TMPBEDGRS, "+>>".$temp_bedgr_s or die $!;
    open TMPBEDGRAS, "+>>".$temp_bedgr_as or die $!;


    # Store
    foreach my $start (sort {$a <=> $b} keys %{$hits->{$chromosome}}){

        # Vars
        my $size = 1; #$bins->{$start}->{"o"} - $bins->{$start}->{"a"} + 1;
        my $plus_count = ($hits->{$chromosome}->{$start}->{'+'}) ? $hits->{$chromosome}->{$start}->{'+'}/$size : 0;
        my $min_count =  ($hits->{$chromosome}->{$start}->{'-'}) ? $hits->{$chromosome}->{$start}->{'-'}/$size : 0;
        my $start_pos = $start;
        #Convert to 0-based (BED=0-based instead of SAM=1-based)
        my $start_pos_Obased = $start_pos -1;
        my $sign;

        my $strand;
        # To db
        if ($min_count != 0) {
            $strand = "-1";
            $sign ="-";
            $min_count = sprintf("%.3f", $min_count);

            print TMP $chromosome.",".$strand.",".$start_pos.",".$min_count."\n";
            print TMPBED "chr$chromosome\t$start_pos_Obased\t$start_pos\t \t$min_count\t$sign\t0\t0\t239,34,5\t\t\t\t\t\n";
            print TMPBEDGRAS "chr$chromosome\t$start_pos_Obased\t$start_pos\t$sign$min_count\n";

        }
        if ($plus_count != 0) {
            $strand = "1";
            $sign ="+";
            $plus_count = sprintf("%.3f", $plus_count);

           print TMP $chromosome.",".$strand.",".$start_pos.",".$plus_count."\n";
            print TMPBED "chr$chromosome\t$start_pos_Obased\t$start_pos\t \t$plus_count\t$sign\t0\t0\t23,170,35\t\t\t\t\t\n";
            print TMPBEDGRS "chr$chromosome\t$start_pos_Obased\t$start_pos\t$plus_count\n";


        }
    }

    foreach my $start (sort {$a <=> $b} keys %{$hits_splitRPF->{$chromosome}}){
        foreach my $RPF (sort {$a <=> $b} keys %{$hits_splitRPF->{$chromosome}->{$start}}){
            # Vars
            my $size = 1; #$bins->{$start}->{"o"} - $bins->{$start}->{"a"} + 1;
            my $plus_count = ($hits_splitRPF->{$chromosome}->{$start}->{$RPF}->{'+'}) ? $hits_splitRPF->{$chromosome}->{$start}->{$RPF}->{'+'}/$size : 0;
            my $min_count =  ($hits_splitRPF->{$chromosome}->{$start}->{$RPF}->{'-'}) ? $hits_splitRPF->{$chromosome}->{$start}->{$RPF}->{'-'}/$size : 0;
            #print "$plus_count,$min_count\n"; exit;
            my $start_pos = $start;
            #Convert to 0-based (BED=0-based instead of SAM=1-based)
            my $start_pos_Obased = $start_pos -1;
            my $sign;

            my $strand;
            # To db
            if ($min_count != 0) {
                $strand = "-1";
                $sign ="-";
                $min_count = sprintf("%.3f", $min_count);

                print TMP_SPLITRPF $chromosome.",".$strand.",".$start_pos.",".$RPF.",".$min_count."\n";

            }
            if ($plus_count != 0) {
                $strand = "1";
                $sign ="+";
                $plus_count = sprintf("%.3f", $plus_count);

                print TMP_SPLITRPF $chromosome.",".$strand.",".$start_pos.",".$RPF.",".$plus_count."\n";

            }
        }
    }


    close(TMP);
    close(TMP_SPLITRPF);
    close(TMPBED);
    close(TMPBEDGRS);
    close(TMPBEDGRAS);
}

### STORE IN FILE PER CHR (RNA), only for count table ###
sub store_in_file_per_chr_RNA {

    #Catch
    my $hits = $_[0];
    my $dbh = $_[1];
    my $seqFileName = $_[2];
    my $chromosome = $_[3];
    my $run_name = $_[4];
    my $truseq = $_[5];

    my $directory = $work_dir."/".$mapper."/".$seqFileName."/";

    #Init temporary csv-file
    my $temp_csv = $TMP."/genomic/".$run_name."_".$seqFileName."_".$chromosome."_tmp.csv";

    open TMP, "+>>".$temp_csv or die $!;

    # Store
    foreach my $start (sort {$a <=> $b} keys %{$hits->{$chromosome}}){

        #Vars
        my $size = 1;
        my $plus_count = ($hits->{$chromosome}->{$start}->{'+'}) ? $hits->{$chromosome}->{$start}->{'+'}/$size : 0;
        my $min_count = ($hits->{$chromosome}->{$start}->{'-'}) ? $hits->{$chromosome}->{$start}->{'-'}/$size : 0;
        my $start_pos = $start;
        my $strand;

        #Write to files
        if ($min_count != 0){
            if ($truseq eq 'Y'){
                $strand="1";
            } else {
                $strand="-1";
            }
            $min_count = sprintf("%.3f", $min_count);

            print TMP $chromosome.",".$strand.",".$start_pos.",".$min_count."\n";
        }
        if ($plus_count != 0){
            if ($truseq eq "Y"){
                $strand="-1";
            } else {
                $strand="1";
            }
            $plus_count = sprintf("%.3f", $plus_count);

            print TMP $chromosome.",".$strand.",".$start_pos.",".$plus_count."\n";
        }
    }

    close(TMP);

}

#Parse dme RIBO_CIGARS to obtain pruned alignment,read mapping length and total intronic length for each pruned alignment position
sub parse_dme_RIBO_CIGAR {

    #Catch
    my $CIGAR = $_[0];
    my $strand = $_[1];

    my $CIGAR_SPLIT = splitCigar($CIGAR);
    my $CIGAR_SPLIT_STR = [];
    @$CIGAR_SPLIT_STR = ($strand eq "-") ? reverse @$CIGAR_SPLIT : @$CIGAR_SPLIT;
    my $op_total = @$CIGAR_SPLIT_STR;

    my $genmatchL = 0;
    my $op_count = 0;
    #To keep track of total length of genomic + intron (negative strand, reverse position)
    my $extra_for_min_strand = 0;
    #Loop over operation to get total mapping length to calculate the A-site offset
    # and to get total extra length for min_strand (i.e. S(not 5adapt nor 1stTRIM), N (splicing), M, D,I)
    foreach my $operation (@$CIGAR_SPLIT_STR) {
        my $op_length = $operation->[0];
        my $op_type = $operation->[1];
        $op_count++;

        if($op_type =~ /^S$/) {
            #Trim leading substitution if only 1 substitution @ 5'
            if ($op_count == 1 && $op_length == 1) {
                next;
            }
            #Clip trailing adaptor substitution
            elsif ($op_count == $op_total) {
                next;
            }
            #Other substitutions are added to RIBO-read genomic-match length
            #And also added to the total matching count
            else {
                $genmatchL = $genmatchL + $op_length;
                $extra_for_min_strand = $extra_for_min_strand + $op_length;
            }
        }
        #Sum matching operations until the offset is reached, then change status to "Y"
        elsif($op_type =~ /^M$/) {
            $genmatchL = $genmatchL + $op_length;
            $extra_for_min_strand = $extra_for_min_strand + $op_length;
        }
        #Insertions elongate the readL and the insertion size is added to the total matching count
        elsif($op_type =~ /^I$/) {
            $genmatchL = $genmatchL + $op_length;
            $extra_for_min_strand = $extra_for_min_strand + $op_length;
        }
        #Splice intronic regions are added to the extra_for_min_strand
        elsif($op_type =~ /^N$/) {
            $extra_for_min_strand = $extra_for_min_strand + $op_length;
        }
    }
    #print "total mapped sequence = $genmatchL\n";
    my $offset = 12;
    my $match_count_total = 0;
    $op_count = 0;

    #Create hash for pruned alignment.
    my $prunedalignmentL = $genmatchL - (2 * $offset);
    my $prunedalignment = {};
    my $pruned_alignment_position;
    my $intron_total = 0;

    #Return if genmatchL too short
    if ($prunedalignmentL <= 0) {
        return ($offset,$genmatchL,0,$extra_for_min_strand,$prunedalignmentL,$prunedalignment);
    }

    #Run over each pruned alignment position
    for(my $i=1;$i<=$prunedalignmentL;$i++){

        my $offset_covered = "N";
        $intron_total = 0;
        $match_count_total = 0;
        $op_count = 0;
        $pruned_alignment_position = $offset + $i;

        #Loop over operations to caculate the total intron length
        foreach my $operation (@$CIGAR_SPLIT_STR) {
            my $op_length = $operation->[0];
            my $op_type = $operation->[1];
            $op_count++;

            if($op_type =~ /^S$/) {
                #Trim leading substitution if only 1 substitution @ 5'
                if ($op_count == 1 && $op_length == 1) {
                    next;
                }
                #Clip trailing adaptor substitution
                elsif ($op_count == $op_total) {
                    next;
                }
                #Other substitutions are added to RIBO-read genomic-match length
                #And also added to the total matching count
                else {
                    $match_count_total = $match_count_total + $op_length;
                    if ($match_count_total >= $pruned_alignment_position) {
                        $offset_covered = "Y";
                        last;
                    }
                }
            }
            #Sum matching operations until the offset is reached, then change status to "Y"
            elsif($op_type =~ /^M$/) {
                $match_count_total = $match_count_total + $op_length;
                if ($match_count_total >= $pruned_alignment_position) {
                    $offset_covered = "Y";
                    last;
                }
            }
            #Sum intronic region lengths untill the offset has been covered by the matching operations
            elsif($op_type =~ /^N$/ && $offset_covered eq "N") {
                $intron_total = $intron_total + $op_length;
            }
            #Deletion are not counted for the readL
            elsif($op_type =~ /^D$/) {
                next;
            }
            #Insertions elongate the readL and the insertion size is added to the total matching count
            elsif($op_type =~ /^I$/) {
                $match_count_total = $match_count_total + $op_length;
                if ($match_count_total >= $pruned_alignment_position) {
                    $offset_covered = "Y";
                    last;
                }
            }
        }

        #Save in prunedalignment_hash
        $prunedalignment->{$pruned_alignment_position}{'offset_covered'} = $offset_covered;
        $prunedalignment->{$pruned_alignment_position}{'intron_total'} = $intron_total;
    }

    return($offset,$genmatchL,$intron_total,$extra_for_min_strand,$prunedalignmentL,$prunedalignment);
}

#Parse RIBO_CIGARS to obtain offset,genomic read mapping length and total intronic length before offset is reached
sub parse_RIBO_CIGAR {

    #Catch
    my $CIGAR = $_[0];
    my $strand = $_[1];
    my $offset_hash = $_[2];

    my $CIGAR_SPLIT = splitCigar($CIGAR);
    my $CIGAR_SPLIT_STR = [];
    @$CIGAR_SPLIT_STR = ($strand eq "-") ? reverse @$CIGAR_SPLIT : @$CIGAR_SPLIT;
    my $op_total = @$CIGAR_SPLIT_STR;

    my $genmatchL = 0;
    my $op_count = 0;
    #To keep track of total length of genomic + intron (negative strand, reverse position)
    my $extra_for_min_strand = 0;
    #Loop over operation to get total mapping length to calculate the A-site offset
    # and to get total extra length for min_strand (i.e. S(not 5adapt nor 1stTRIM), N (splicing), M, D,I)
    foreach my $operation (@$CIGAR_SPLIT_STR) {
        my $op_length = $operation->[0];
        my $op_type = $operation->[1];
        $op_count++;

        if($op_type =~ /^S$/) {
            #Trim leading substitution if only 1 substitution @ 5'
            if ($op_count == 1 && $op_length == 1) {
                next;
            }
            #Clip trailing adaptor substitution
            elsif ($op_count == $op_total) {
                next;
            }
            #Other substitutions are added to RIBO-read genomic-match length
            #And also added to the total matching count
            else {
                $genmatchL = $genmatchL + $op_length;
                $extra_for_min_strand = $extra_for_min_strand + $op_length;
            }
        }
        #Sum matching operations until the offset is reached, then change status to "Y"
        elsif($op_type =~ /^M$/) {
            $genmatchL = $genmatchL + $op_length;
            $extra_for_min_strand = $extra_for_min_strand + $op_length;
        }
        #Insertions elongate the readL and the insertion size is added to the total matching count
        elsif($op_type =~ /^I$/) {
            $genmatchL = $genmatchL + $op_length;
            $extra_for_min_strand = $extra_for_min_strand + $op_length;
        }
        #Splice intronic regions are added to the extra_for_min_strand
        elsif($op_type =~ /^N$/) {
            $extra_for_min_strand = $extra_for_min_strand + $op_length;
        }
    }
    #print "total mapped sequence = $genmatchL\n";
    my $offset = get_offset($genmatchL, $offset_hash);
    my $match_count_total = 0;
    $op_count = 0;
    my $offset_covered = "N";
    my $intron_total = 0;
    #Loop over operations to caculate the total intron length
    foreach my $operation (@$CIGAR_SPLIT_STR) {
        my $op_length = $operation->[0];
        my $op_type = $operation->[1];
        $op_count++;
        #print "$op_type,$op_length\n";
        if($op_type =~ /^S$/) {
            #Trim leading substitution if only 1 substitution @ 5'
            if ($op_count == 1 && $op_length == 1) {
                next;
            }
            #Clip trailing adaptor substitution
            elsif ($op_count == $op_total) {
                next;
            }
            #Other substitutions are added to RIBO-read genomic-match length
            #And also added to the total matching count
            else {
                $match_count_total = $match_count_total + $op_length;
                if ($match_count_total >= $offset) {
                    $offset_covered = "Y";
                    last;
                }
            }
        }
        #Sum matching operations until the offset is reached, then change status to "Y"
        elsif($op_type =~ /^M$/) {
            $match_count_total = $match_count_total + $op_length;
            if ($match_count_total >= $offset) {
                $offset_covered = "Y";
                last;
            }
        }
        #Sum intronic region lengths untill the offset has been covered by the matching operations
        elsif($op_type =~ /^N$/ && $offset_covered eq "N") {
            $intron_total = $intron_total + $op_length;

        }
        #Deletion are not counted for the readL
        elsif($op_type =~ /^D$/) {
            next;
            #$genmatchL = $genmatchL - $op_length;
        }
        #Insertions elongate the readL and the insertion size is added to the total matching count
        elsif($op_type =~ /^I$/) {
            $match_count_total = $match_count_total + $op_length;
            if ($match_count_total >= $offset) {
                $offset_covered = "Y";
                last;
            }
        }
        #print "$match_count_total,$offset_covered,$offset\n";
    }

    return($offset,$genmatchL,$intron_total,$extra_for_min_strand)
}

### Parse RNA_CIGARs to obtain offset, genomic read mapping length and total intronic length before offset is reached
sub parse_RNA_CIGAR {

    #Catch
    my $CIGAR = $_[0];
    my $strand = $_[1];

    my $CIGAR_SPLIT = splitCigar($CIGAR);
    my $CIGAR_SPLIT_STR = [];
    @$CIGAR_SPLIT_STR = ($strand eq "-") ? reverse @$CIGAR_SPLIT : @$CIGAR_SPLIT;
    my $op_total = @$CIGAR_SPLIT_STR;

    #Init
    my $genmatchL = 0;
    my $op_count = 0;
    my $extra_for_min_strand = 0; #To keep track of total length of genomic + intron (negative strand, reverse position)

    #Loop over operations to get total mapping length and to get the total
    #extra length for min_strand (i.e. S(not 5adapt nor 1st TRIM), N (splicing), M, D, I)
    foreach my $operation (@$CIGAR_SPLIT_STR) {
        my $op_length = $operation->[0];
        my $op_type = $operation->[1];
        $op_count++;

        if($op_type =~ /^S$/) {
            #Trim leading substitution if only 1 substitution @ 5'
            if ($op_count ==1 && $op_length == 1){
                next;
            }
            #Clip trailing adaptor substitution
            elsif ($op_count == $op_total) {
                next;
            }
            #Other substitutions are added to RNAread genomic-match length
            #And also added to the total matching count
            else{
                $genmatchL = $genmatchL + $op_length;
                $extra_for_min_strand = $extra_for_min_strand + $op_length;
            }
        }
        #Sum matching operations
        elsif($op_type =~ /^M$/) {
            $genmatchL = $genmatchL + $op_length;
            $extra_for_min_strand = $extra_for_min_strand + $op_length;
        }
        #Insertions elongate the readL and the insertion size is added to the total matching count
        elsif($op_type =~ /^I$/) {
            $genmatchL = $genmatchL + $op_length;
            $extra_for_min_strand = $extra_for_min_strand + $op_length;
        }
        #Splice intronic regions are added to the extra length needed for the min strand
        elsif($op_type =~ /^N$/) {
            $extra_for_min_strand = $extra_for_min_strand + $op_length;
        }
    }

    #Init
    my $offset = 26; #For RNA, reads are pinpointed at nucleotide 26 as in RiboTaper paper (Calviello et al. 2015).
    my $match_count_total = 0;
    $op_count = 0; #Reset to zero
    my $offset_covered = "N";
    my $intron_total = 0;

    #Loop over operations to calculate the total intron length
    foreach my $operation (@$CIGAR_SPLIT_STR) {
        my $op_length = $operation->[0];
        my $op_type = $operation->[1];
        $op_count++;

        if ($op_type =~ /^S$/) {
            #Trim leading substitution if only 1 substitution @5'
            if ($op_count == 1 && $op_length == 1) {
                next;
            }
            #Clip trailing adaptor substituion
            elsif ($op_count == $op_total) {
                next;
            }
            #Other substitutions are added to RNAread genomic-match length
            #And also adde to the total matching count
            else {
                $match_count_total = $match_count_total = $op_length;
                if ($match_count_total >= $offset) {
                    $offset_covered = "Y";
                    last;
                }
            }
        }
        #Sum matching operations until the offset is reached, then change the status to "Y"
        elsif($op_type =~ /^M$/) {
            $match_count_total = $match_count_total + $op_length;
            if ($match_count_total >= $offset) {
                $offset_covered = "Y";
                last;
            }
        }
        #Sum intronic region lengths until the offset has been covered by matching operations
        elsif($op_type =~ /^N$/ && $offset_covered eq "N") {
            $intron_total = $intron_total + $op_length;
        }
        #Deletions are not counted for the readL
        elsif($op_type =~ /^D$/) {
            next;
        }
        #Insertions elongate the readL and the insertion size is added to the total matching count
        elsif($op_type =~ /^I$/) {
            $match_count_total = $match_count_total + $op_length;
            if ($match_count_total >= $offset) {
                $offset_covered="Y";
                last;
            }
        }
    }

    return ($offset,$genmatchL,$intron_total,$extra_for_min_strand);
}


### CALCULATE A-SITE OFFSET ###
sub get_offset {

    #Catch
    my $len = $_[0];
    my $offset_hash = $_[1];
    
    #Init
    my $offset;

    if($len<$offset_hash->{"min"}){
        $offset = $offset_hash->{$offset_hash->{"min"}};
    } elsif($len>$offset_hash->{"max"}){
        $offset = $offset_hash->{$offset_hash->{"max"}};
    } else {
        $offset = $offset_hash->{$len};
    }

    return($offset);
}


# Given a Cigar string return a double array
# such that each sub array contiains the [ length_of_operation, Cigar_operation]
sub splitCigar {
    my $cigar_string = shift;
    my @returnable;
    my (@matches) = ($cigar_string =~ /(\d+\w)/g);
    foreach (@matches) {
        my @operation = ($_ =~ /(\d+)(\w)/);
        push @returnable, \@operation;
    }
    return \@returnable;
}

##Get paramters from arguments table

sub get_ARG_vars{
    # Catch
    my $dsn = $_[0];
    my $user = $_[1];
    my $pw = $_[2];
    my $offset_option = $_[3];
    
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
    
    $query = "select value from arguments where variable = \'nr_of_cores\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $cores = $sth->fetch()->[0];
    $sth->finish();
    
    $query = "select value from arguments where variable = \'seqFileName1\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $seqFileName1 = $sth->fetch()->[0];
    $sth->finish();
    
    $query = "select value from arguments where variable = \'seqFileName2\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $seqFileName2 = $sth->fetch()->[0];
    $sth->finish();
    
    $query = "select value from arguments where variable = \'mapper\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $mapper = $sth->fetch()->[0];
    $sth->finish();
    
    $query = "select value from arguments where variable = \'unique\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $unique = $sth->fetch()->[0];
    $sth->finish();
    
    $query = "select value from arguments where variable = \'rpf_split\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $rpf_split = $sth->fetch()->[0];
    $sth->finish();
    
    $query = "select value from arguments where variable = \'firstRankMultiMap\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $FirstRankMultiMap = $sth->fetch()->[0];
    $sth->finish();
    
    $query = "select value from arguments where variable = \'truseq\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $truseq = $sth->fetch()->[0];
    $sth->finish();
    
    $query = "select value from arguments where variable = \'readtype\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $readtype = $sth->fetch()->[0];
    $sth->finish();
    
    $query = "select value from arguments where variable = \'out_bg_s_untr\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $out_bg_s_untr = $sth->fetch()->[0];
    $sth->finish();
    
    $query = "select value from arguments where variable = \'out_bg_as_untr\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $out_bg_as_untr = $sth->fetch()->[0];
    $sth->finish();
    
    $query = "select value from arguments where variable = \'out_bg_s_tr\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $out_bg_s_tr = $sth->fetch()->[0];
    $sth->finish();
    
    $query = "select value from arguments where variable = \'out_bg_as_tr\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $out_bg_as_tr = $sth->fetch()->[0];
    $sth->finish();
    
    $query = "select value from arguments where variable = \'out_sam_untr\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $out_sam_untr = $sth->fetch()->[0];
    $sth->finish();
    
    $query = "select value from arguments where variable = \'out_sam_tr\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $out_sam_tr = $sth->fetch()->[0];
    $sth->finish();
    
    $query = "select value from arguments where variable = \'run_name\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $run_name = $sth->fetch()->[0];
    $sth->finish();
    
    $query = "select value from arguments where variable = \'maxmultimap\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $maxmultimap = $sth->fetch()->[0];
    $sth->finish();
    
    $query = "select value from arguments where variable = \'min_l_parsing\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $min_l_count = $sth->fetch()->[0];
    $sth->finish();
    
    $query = "select value from arguments where variable = \'max_l_parsing\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $max_l_count = $sth->fetch()->[0];
    $sth->finish();
    
    #Init
    my $cst_prime_offset = 0;
    my $min_cst_prime_offset = 0;
    my $max_cst_prime_offset = 0;
    if($offset_option eq "cst_5prime" || $offset_option eq "cst_3prime"){
        $query = "select value from arguments where variable = \'cst_prime_offset\'";
        $sth = $dbh_results->prepare($query);
        $sth->execute();
        $cst_prime_offset = $sth->fetch()->[0];
        $sth->finish();
        
        $query = "select value from arguments where variable = \'min_cst_prime_offset\'";
        $sth = $dbh_results->prepare($query);
        $sth->execute();
        $min_cst_prime_offset = $sth->fetch()->[0];
        $sth->finish();
        
        $query = "select value from arguments where variable = \'max_cst_prime_offset\'";
        $sth = $dbh_results->prepare($query);
        $sth->execute();
        $max_cst_prime_offset = $sth->fetch()->[0];
        $sth->finish();
    }

    $dbh_results -> disconnect();
    
    # Return ARG variables
    return($species,$version,$IGENOMES_ROOT,$cores,$seqFileName1,$seqFileName2,$mapper,$unique,$rpf_split,$FirstRankMultiMap,$truseq,$readtype,$out_bg_s_untr,$out_bg_as_untr,$out_bg_s_tr,$out_bg_as_tr,$out_sam_untr,$out_sam_tr,$run_name,$maxmultimap,$min_l_count,$max_l_count,$cst_prime_offset,$min_cst_prime_offset,$max_cst_prime_offset);
}

##Help text##
sub print_help_text {
    
    my $help_string = "\n\nParsing of mapping results (PROTEOFORMER)
    
Alignment results will in this tool be parsed. Based on the P-site offset option, P-site offsets will be gathered from a list of standard offsets or from the results calculated by Plastid. With these offsets, the alignments in the SAM file will be pinpointed on a certain base position. Counts per base positiion will then be saved to the SQLite table (both a table with total counts per position as also a table of counts per position and per RPF length). Furthermore, BED and BEDGRAPH files will be generated for visual inspection of the counts in genome browser tools like UCSC.

Example:
    perl mapping_parsing.pl --out_sqlite SQLite/results.db --offset plastid
    
Input arguments:
        --work_dir                      The working directory (default: CWD env setting)
        --tmp                           The folder where temporary files are stored (default: work_dir/tmp)
        --out_sqlite                    SQLite results DB (default: work_dir/SQLite/results.db)
        --offset                        Offset option (standard, from_file or plastid) (default: standard)
                                            Standard: if you use the standard offset lengths used by Ingolia et al. (2009)
                                            From_file: from tab separated txt file with rpf_length-tab-offset format
                                            Plastid: if you have run the mapping_plastid module first
                                            cst_5prime: use constant offsets relative to the 5prime read end
                                            cst_3prime: use constant offsets relative to the 3prime read end
        --offset_file_untr              Offset input file for untreated data (mandatory if offset=from_file)
        --offset_file_tr                Offset input file for treated data (mandatory if offset=from_file)
        --help                          Print help message

";
    
    print $help_string;
}
