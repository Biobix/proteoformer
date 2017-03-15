$|=1;

####### Mapping Quality Control

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

# nohup perl ./mappingQC.pl --samfile STAR/fastq1/untreat.sam --treated 'untreated' --cores 20 --result_db SQLite/results.db --ens_db SQLite/ENS_mmu_82.db --offset plastid > nohup_mappingqc.txt &

my($work_dir,$sam,$treated,$cores,$resultdb,$tmpfolder,$maxmultimap,$ens_db,$offset_option,$offset_file,$offset_img,$output_folder,$tool_dir,$html,$zip);

GetOptions(
"work_dir:s" => \$work_dir,             # The working directory                                         Optional argument (default: CWD)
"samfile=s"=>\$sam,                     # The samfile to do the analysis on                             Mandatory argument
"treated=s"=>\$treated ,                # Wheter the samfile is from the treated or untreated sample    Optional argument (untreated/treated, default untreated)
"cores=i"=>\$cores,                     # The amount of cores to use                                    Optional argument (default: 5)
"result_db=s"=>\$resultdb,              # The result db with mapping results                            Mandatory argument
"tmp:s"=>\$tmpfolder,                   # The tmp folder                                                Optional argument (default: CWD/tmp)
"maxmultimap=i"=>\$maxmultimap,         # The maximum multimapped positions for parsing                 Optional argument (default: 16)
"ens_db=s"=>\$ens_db,                   # The Ensembl db for annotation                                 Mandatory argument
"offset:s" =>\$offset_option,           # The offset source for parsing alignments                      Optional argument (default: standard)
"offset_file:s" =>\$offset_file,        # The offsets input file                                        Mandatory if offset option equals 'from_file'
"offset_img=s" =>\$offset_img,          # The offsets image from plastid                                Mandatory if offset option equals 'plastid'
"output_folder:s" => \$output_folder,   # The output folder for storing output images                   Optional argument (default: CWD/mappingQC_(un)treated/)
"tool_dir:s" => \$tool_dir,             # The directory with all necessary tools                        Optional argument (default: CWD/mqc_tools/)
"html:s" => \$html,                     # The output html file name                                     Optional argument (default: CWD/mappingqc_out.html)
"zip:s" => \$zip                        # The output zip file name of the output folder                 Optional argument (default CWD/mappingQC_(un)treated.zip )
);


###########################################################################
#Check all input variable and/or get default values and set extra variables
###########################################################################

my $CWD             = getcwd;

# comment on these
if ($work_dir){
    print "Working directory                                        : $work_dir\n";
} else {
    $work_dir = $CWD;
    print "Working directory                                        : $work_dir\n";
}
my $TMP             = ($ENV{'TMP'}) ? $ENV{'TMP'} : ($tmpfolder) ? $tmpfolder : "$CWD/tmp"; # (1) get the TMP environment variable, (2) get the $tmpfolder variable, (3) get current_working_dir/tmp
print "The following tmpfolder is used                          : $TMP\n";

#Check if tmpfolder exists, if not create it...
if (!-d "$TMP") {
    system ("mkdir ". $TMP);
}
if ($sam){
    print "the input sam file                                                : $sam\n";
} else {
    die "\nDon't forget to pass the sam file!\n\n";
}
if ($treated){
    if ($treated eq "treated" || $treated eq "untreated"){
        print "Sample treatment                                  : $treated\n";
    } else {
        print "ERROR: treated argument should be 'untreated' or 'treated'!\n";
        die;
    }
} else {
    $treated = "untreated";
    print "Sample treatment                                  : $treated\n";
}
if ($resultdb){
    print "the results DB                               : $resultdb\n";
} else {
    die "\nDon't forget to pass the results DB!\n\n";
}
if ($maxmultimap){
    print "Maximun number of loci for reads to be acceptable        : $maxmultimap\n";
} else {
    $maxmultimap = 16;
    print "Maximun number of loci for reads to be acceptable        : $maxmultimap\n";
}
if ($ens_db){
    print "the Ensembl DB                               : $ens_db\n";
} else {
    die "\nDon't forget to pass the Ensembl DB!\n\n";
}
if ($offset_option) {
    if ($offset_option eq "standard" || $offset_option eq "from_file" || $offset_option eq "plastid") {
        print "Offset source                                            : $offset_option\n";
    } else {
        die "Offset argument needs to be \" standard\", \"from_file\" or \"plastid\"!";
    }
} else {
    $offset_option = "standard";
    print "Offset source                                            : $offset_option\n";
}
if ($offset_option eq "from_file"){
    if ($offset_file) {
        print "Offset input file                                        : $offset_file\n";
    } else {
        die "Do not forget the offset input file if offset argument is \"from_file\"!";
    }
} else {
    $offset_file = "";
}
if ($offset_option eq "plastid"){
    if ($offset_img){
        print "Offset image                                              : $offset_img\n";
    } else {
        die "Do not forget the offset image if offset argument is \"plastid\"!;
    }
}

if ($output_folder){
    print "The output folder is set to     : $output_folder\n";
} else {
    $output_folder = $work_dir."/mappingQC_".$treated."/";
}
if ($tool_dir){
    print "The tool directory is set to:    $tool_dir\n";
} else {
    $tool_dir = $work_dir."/mqc_tools/";
    print "The tool directory is set to     : $tool_dir\n";
}
if ($html){
    print "Output html file name                       : $html\n";
} else {
    $html = $work_dir."/mappingqc_out.html";
    print "Output html file name                       : $html\n";
}
if ($zip){
    print "Output zip file name                         : $zip\n";
} else {
    $zip = $work_dir."/mappingQC_".$treated.".zip";
    print "Output zip file name                         : $zip\n";
}

my $dsn_sqlite_results = "DBI:SQLite:dbname=$resultdb";
my $us_sqlite_results  = "";
my $pw_sqlite_results  = "";

# Get arguments vars
my ($species,$version,$IGENOMES_ROOT) = get_ARG_vars($resultdb,$us_sqlite_results,$pw_sqlite_results);

# Igenomes
print "The following igenomes folder is used			: $IGENOMES_ROOT\n";

# Cores
if ($cores) {
    print "Number of cores to use for analysis			: $cores\n";
} else {
    $cores = 5;
    print "Number of cores to use for analysis			: $cores\n";
}

#Conversion for species terminology
my $spec = ($species eq "mouse") ? "Mus_musculus" : ($species eq "human") ? "Homo_sapiens" : ($species eq "arabidopsis") ? "Arabidopsis_thaliana" : ($species eq "fruitfly") ? "Drosophila_melanogaster" : "";
my $spec_short = ($species eq "mouse") ? "mmu" : ($species eq "human") ? "hsa" : ($species eq "arabidopsis") ? "ath" : ($species eq "fruitfly") ? "dme" : "";
#Old mouse assembly = NCBIM37, new one is GRCm38. Old human assembly = GRCh37, the new one is GRCh38
my $assembly = (uc($species) eq "MOUSE" && $version >= 70 ) ? "GRCm38"
: (uc($species) eq "MOUSE" && $version < 70 ) ? "NCBIM37"
: (uc($species) eq "HUMAN" && $version >= 76) ? "GRCh38"
: (uc($species) eq "HUMAN" && $version < 76) ? "GRCh37"
: (uc($species) eq "ARABIDOPSIS") ? "TAIR10"
: (uc($species) eq "FRUITFLY" && $version < 79) ? "BDGP5"
: (uc($species) eq "FRUITFLY" && $version >= 79) ? "BDGP6" : "";

# Get chromosomes and correct coord_system_id
print "Get chromosomes and coord_system_id...\n";
my $chromosome_sizes; my $coord_system_id; my @ch;

$chromosome_sizes = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Annotation/Genes/ChromInfo.txt";
## Get chromosome sizes and cDNA identifiers #############
print "Getting chromosome sizes and cDNA to chromosome mappings ...\n";
my %chr_sizes = %{get_chr_sizes($chromosome_sizes)};
$coord_system_id = get_coord_system_id($ens_db,$assembly);

my $chrom_dir = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Sequence/Chromosomes";
my $BIN_chrom_dir = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Sequence/Chromosomes_BIN";

#Make tmp folder specific for mappingqc
if (! -e $TMP."/mappingqc_".$treated){
    system("mkdir ".$TMP."/mappingqc_".$treated);
}

####### FOR TESTING on Y chromsome #######
#my %chr_sizesY;
#$chr_sizesY{'Y'} = $chr_sizes{'Y'};
#%chr_sizes = %chr_sizesY;
##########################################

########
# MAIN #
########

# Start time
my $start = time;

print "Get chromosomes... \n";
## Get chromosomes based on seq_region_id ##
# Sqlite Ensembl
my $db_ENS  = $ens_db;
my $dsn_ENS = "DBI:SQLite:dbname=$db_ENS";
my $us_ENS  = "";
my $pw_ENS  = "";
my $chrs = get_chrs($dsn_ENS,$us_ENS,$pw_ENS,\%chr_sizes,$assembly);

# Create binary chromosomes if they don't exist
print "Checking/Creating binary chrom files ...\n";
if (!-d "$BIN_chrom_dir") {
    create_BIN_chromosomes($BIN_chrom_dir,$cores,$chrs,$work_dir,$TMP);
}

my @splitsam = split(/\//, $sam );
my $samFileName = $splitsam[$#splitsam];
@splitsam = split(/\./,$samFileName);
$samFileName = $splitsam[0];
my $samfilechr1 = $TMP."/mappingqc/".$samFileName."_1.sam";

if (-e $samfilechr1){
    print "Splitted sam files already exist...\n";
} else {
    print "Splitting genomic mapping per chromosome...\n";
    split_SAM_per_chr(\%chr_sizes,$work_dir,$sam,$treated);
}

# Construct p offset hash
my $offset_hash = {};
if($offset_option eq "plastid"){
    #Offset from plastid table
    #Init
    $offset_hash->{"min"} = 1000;
    $offset_hash->{"max"} = 0;
    #
    #Connect to result db
    my $dbh_plastid = dbh($dsn_sqlite_results, $us_sqlite_results, $pw_sqlite_results);
    
    #Query and parse plastid table
    my $query = "SELECT * FROM p_offsets_".$treated.";";
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

#Write offsets to csv for output html file
offsets_to_csv($offset_hash, $treated, $TMP);

# Init multi core
my $pm = new Parallel::ForkManager($cores);
print "   Using ".$cores." core(s)\n   ---------------\n";

foreach my $chr (keys %chr_sizes){
    
    ### Start parallel process
    $pm->start and next;
    
    ### DBH per process
    my $dbh = dbh($dsn_sqlite_results,$us_sqlite_results,$pw_sqlite_results);
    
    ### RIBO parsing
    RIBO_parsing_genomic_per_chr($work_dir,$sam,$chr,$ens_db,$coord_system_id, $offset_hash, $treated);
    
    ### Finish
    print "* Finished chromosome ".$chr."\n";
    $dbh->disconnect();
    $pm->finish;
}

# Finish all subprocesses
$pm->wait_all_children;
print "\n";

print "Prepare data for plotting modules\n";


## RPF PHASE TABLE ##
print "\tRPF phase table\n";
#Count rpf phase table of all chromosomes together
my $temp_csv_rpf_phase = $TMP."/mappingqc_".$treated."/rpf_phase.csv";
system("touch ".$temp_csv_rpf_phase);

my $rpf_phase = {};
for (my $i=22;$i<=34;$i++){
    for (my $j=0;$j<=2;$j++){
        $rpf_phase->{$i}->{$j} = 0;
    }
}
foreach my $chr (keys %chr_sizes){
    my $infile = $TMP."/mappingqc_".$treated."/rpf_phase_".$chr.".csv";
    open(IN,"<".$infile) or die $!;
    while(my $line = <IN>){
        my @linesplit = split(',',$line);
        $rpf_phase->{$linesplit[0]}->{0} = $rpf_phase->{$linesplit[0]}->{0} + $linesplit[1];
        $rpf_phase->{$linesplit[0]}->{1} = $rpf_phase->{$linesplit[0]}->{1} + $linesplit[2];
        $rpf_phase->{$linesplit[0]}->{2} = $rpf_phase->{$linesplit[0]}->{2} + $linesplit[3];
    }
    close(IN);
    system("rm -rf ".$infile);
}

#Write rpf phase table to temp csv
open(OUT_PHASE, "+>>".$temp_csv_rpf_phase);
foreach my $rpf (keys %{$rpf_phase}){
    print OUT_PHASE $rpf.",".$rpf_phase->{$rpf}->{0}.",".$rpf_phase->{$rpf}->{1}.",".$rpf_phase->{$rpf}->{2}."\n";
}
close(OUT_PHASE);

## PHASE RELATIVE POSITION DISTRIBUTION
print "\tPhase - relative position distribution\n";
#Cat phase-position tmp files
my $temp_csv_all_pos = $TMP."/mappingqc_".$treated."/pos_table_all.csv";
system("touch ".$temp_csv_all_pos);

foreach my $chr (keys %chr_sizes){
    my $temp_csv_chr_pos = $TMP."/mappingqc_".$treated."/phase_position_".$chr.".csv";
    system("cat ".$temp_csv_chr_pos." >> ".$temp_csv_all_pos);
    system("rm -rf ".$temp_csv_chr_pos);
}

## TRIPLET IDENTITY PHASE FILE
print "\tTriplet identity distributions\n";
#Read in chr tmp files
my $triplet_phase = {};
foreach my $chr (keys %chr_sizes){
    my $infile = $TMP."/mappingqc_".$treated."/triplet_phase_".$chr.".csv";
    open(IN, "< ".$infile) or die $!;
    while(my $line = <IN>){
        chomp($line);
        my @linesplit = split(',',$line);
        my $triplet = $linesplit[0];
        my $phase = $linesplit[1];
        my $count = $linesplit[2];
        if(exists $triplet_phase->{$triplet}->{$phase}){
            $triplet_phase->{$triplet}->{$phase} = $triplet_phase->{$triplet}->{$phase} + $count;
        } else {
            $triplet_phase->{$triplet}->{$phase} = $count;
        }
    }
    close(IN);
    system("rm -rf ".$infile);
}

#Write total file for triplet identity
my $temp_total_triplet = $TMP."/mappingqc_".$treated."/total_triplet.csv";
open(OUT_TOTAL_TRIPLET, "+>> ".$temp_total_triplet);
foreach my $triplet (keys %{$triplet_phase}){
    foreach my $phase (keys %{$triplet_phase->{$triplet}}){
        print OUT_TOTAL_TRIPLET $triplet.",".$phase.",".$triplet_phase->{$triplet}->{$phase}."\n";
    }
}
close(OUT_TOTAL_TRIPLET);


## Save some results to results SQLite DB
print "Store certain results in results DB\n";
#Init dbh
my $dbh_results = dbh($dsn_sqlite_results,$us_sqlite_results,$pw_sqlite_results);

#Create table rpf phase
my $query_create_table = "CREATE TABLE IF NOT EXISTS `rpf_phase_".$treated."` (
`RPF` int(10) NOT NULL default '0',
`phase0` int(10) NOT NULL default '0',
`phase1` int(10) NOT NULL default '0',
`phase2` int(10) NOT NULL default '0')";
$dbh_results->do($query_create_table);

#Dump into sqlite
system("sqlite3 -separator , ".$resultdb." \".import ".$temp_csv_rpf_phase." rpf_phase_".$treated."\"")==0 or die "System failed: $?";
system("rm -rf ".$temp_csv_rpf_phase);

#Create table phase_triplet
$query_create_table = "CREATE TABLE IF NOT EXISTS `phase_triplet_".$treated."` (
`triplet` varchar(5) NOT NULL default '',
`phase` varchar(2) NOT NULL default '',
`count` int(10) NOT NULL default '0')";
$dbh_results->do($query_create_table);

#Dump into sqlite
system("sqlite3 -separator , ".$resultdb." \".import ".$temp_total_triplet." phase_triplet_".$treated."\"")==0 or die "System failed: $?";
system("rm -rf ".$temp_total_triplet);

#Disconnect
$dbh_results->disconnect;



#Make output folder
system("mkdir ".$output_folder);
my $treated_meta_gene;

print "\n\n";
print "# Run gene distribution module #\n";
#Run Gene Distribution script
system("perl ".$tool_dir."/gene_distribution.pl --cores ".$cores." --tooldir ".$tool_dir." --in_sqlite ".$resultdb." --output_folder ".$output_folder." --treated ".$treated." --ens_db ".$ens_db);

print "\n\n";
print "# Run metagenic classification module #\n";
# Metagenic classification
system("perl ".$tool_dir."/metagenic_classification.pl --tooldir ".$tool_dir." --cores ".$cores." --in_sqlite ".$resultdb." --ens_db ".$ens_db." --output_folder ".$output_folder." --treated ".$treated);
 
print "\n\n";
print "# Run plot generation and output HTML file creation module #\n";
#Run plot generation Python file
my $python_command = "python ".$tool_dir."/mappingQC.py -r ".$resultdb." -s ".$sam." -t ".$treated." -o ".$output_folder." -p \"".$offset_option."\" -e ".$ens_db." -h ".$html." -z ".$zip;
if ($offset_option eq "plastid"){
    $python_command = $python_command." -i ".$offset_img;
}
system($python_command);

#Clean up
#system("rm -f ".$TMP."/mappingqc_".$treated."/".$samFileName."_*");

# End time
print "   DONE! \n";
my $end = time - $start;
printf("Runtime: %02d:%02d:%02d\n\n",int($end/3600), int(($end % 3600)/60), int($end % 60));



############
# THE SUBS #
############

### RIBO PARSE PER CHR ###
sub RIBO_parsing_genomic_per_chr {
    
    #Catch
    my $work_dir = $_[0];
    my $sam = $_[1];
    my $chr = $_[2];
    my $ens_db = $_[3];
    my $coord_system_id = $_[4];
    my $offset_hash = $_[5];
    my $treated = $_[6];
    
    my @splitsam = split(/\//, $sam );
    my $samFileName = $splitsam[$#splitsam];
    @splitsam = split(/\./,$samFileName);
    $samFileName = $splitsam[0];
    
    #Construct phase library
    my ($phase_lib, $triplet_lib) = construct_phase_lib($chr, $ens_db, $coord_system_id);
    
    #Initialize
    my ($genmatchL,$offset,$start,$intron_total,$extra_for_min_strand);
    my $phase_count_RPF = {};
    for (my $i=22;$i<=34;$i++){
        for (my $j=0;$j<=2;$j++){
            $phase_count_RPF->{$i}->{$j} = 0;
        }
    }
    my $phase_count_file = $TMP."/mappingqc_".$treated."/rpf_phase_".$chr.".csv";
    my $phase_count_triplet = {};
    my $triplet_count_file = $TMP."/mappingqc_".$treated."/triplet_phase_".$chr.".csv";
    my $chr_sam_file = $TMP."/mappingqc_".$treated."/".$samFileName."_".$chr.".sam";
    my $pos_file = $TMP."/mappingqc_".$treated."/phase_position_".$chr.".csv";
    
    open (I,"<".$chr_sam_file) || die "Cannot open ".$chr_sam_file."\n";
    open (OUT_POS, "+>>".$pos_file) or die $!;
    
    while(my $line=<I>){
        #Process alignment line
        my @mapping_store = split(/\t/,$line);
        
        #Get strand specifics
        # Sam flag is bitwise. (0x10 SEQ being reverse complemented)
        # 0x10 = 16 in decimal. -> negative strand.
        my $strand = ($mapping_store[1] & 16) ? "-": "+";
        my $CIGAR = $mapping_store[5];
        #Rewrite strand info
        my $strandAlt;
        if($strand eq "+"){
            $strandAlt = 1;
        } else {
            $strandAlt = -1;
        }
        
        #Parse CIGAR to obtain offset,genomic matching length and total covered intronic region before reaching the offset
        ($offset,$genmatchL,$intron_total,$extra_for_min_strand) = parse_RIBO_CIGAR($CIGAR,$strand,$offset_hash);
        #Determine genomic position based on CIGAR string output and mapping position and direction
        $start = ($strand eq "+") ? $mapping_store[3] + $offset + $intron_total : ($strand eq "-") ? $mapping_store[3] - $offset - $intron_total + $extra_for_min_strand -1 : "";
        
        if($genmatchL>=$offset_hash->{"min"} && $genmatchL<=$offset_hash->{"max"}){
            if(exists $phase_lib->{$strandAlt}->{$start}->{"phase"}){
                #Add for RPF-splitted phase distribution
                $phase_count_RPF->{$genmatchL}->{$phase_lib->{$strandAlt}->{$start}->{"phase"}}++;
                #Print to tmp chr phase-position distribbution
                if(exists $phase_lib->{$strandAlt}->{$start}->{"transcriptomic_pos"}){
                    my $read_phase = $phase_lib->{$strandAlt}->{$start}->{"phase"};
                    print OUT_POS $read_phase.",".$phase_lib->{$strandAlt}->{$start}->{"transcriptomic_pos"}."\n";
                }
            }
            if(exists $triplet_lib->{$strandAlt}->{$start}){
                my $read_phase = $phase_lib->{$strandAlt}->{$start}->{"phase"};
                my $read_triplet = $triplet_lib->{$strandAlt}->{$start};
                if(length($read_triplet)==3){
                    if(exists $phase_count_triplet->{$read_triplet}){
                        $phase_count_triplet->{$read_triplet}->{$read_phase}++;
                    } else {
                        #Initialize if triplet is first time seen
                        $phase_count_triplet->{$read_triplet}->{0} = 0;
                        $phase_count_triplet->{$read_triplet}->{1} = 0;
                        $phase_count_triplet->{$read_triplet}->{2} = 0;
                        $phase_count_triplet->{$read_triplet}->{$read_phase}++;
                    }
                }
            }
        }
    }
    
    #Stop reading out of input files
    close(I);
    
    #Write to chromosomal rpf phase tmp file
    open (OUT_RPF_PHASE, "+>>".$phase_count_file) or die $!;
    
    for (my $i=$offset_hash->{'min'};$i<=$offset_hash->{'max'};$i++){
        print OUT_RPF_PHASE $i.",".$phase_count_RPF->{$i}->{0}.",".$phase_count_RPF->{$i}->{1}.",".$phase_count_RPF->{$i}->{2}."\n";
    }
    close(OUT_RPF_PHASE);
    
    #Write to chromosomal triplet phase tmp file
    open (OUT_TRIPLET_PHASE, "+>>".$triplet_count_file) or die $!;
    
    foreach my $triplet (keys %{$phase_count_triplet}){
        foreach my $phase (keys %{$phase_count_triplet->{$triplet}}){
            print OUT_TRIPLET_PHASE $triplet.",".$phase.",".$phase_count_triplet->{$triplet}->{$phase}."\n";
        }
    }
    
    return;
}

#Construct phase lib out of ensembl info
sub construct_phase_lib{
    
    #Catch
    my $chr = $_[0];
    my $eDB = $_[1];
    my $coord_system_id = $_[2];
    
    #Init
    my $phase_lib = {};
    my $triplet_lib = {};
    my $dsn_sqlite_ens = "DBI:SQLite:dbname=$eDB";
    my $us_sqlite_ens  = "";
    my $pw_sqlite_ens  = "";
    
    #Init dbh
    my $dbh_ens = dbh($dsn_sqlite_ens,$us_sqlite_ens,$pw_sqlite_ens);
    
    #Get seq_region_id
    my $seq_region_id = get_seq_region_id($dbh_ens, $chr, $coord_system_id);
    
    #Get transcripts (canonical protein-coding)
    my $transcripts = get_can_transcripts($dbh_ens, $seq_region_id);
    
    #Get exon structure of each transcript
    foreach my $transcript (@$transcripts){
        my($exon_struct,$strand,$max_tr_rank) = get_exon_struct_transcript($dbh_ens, $transcript, $chr);
        ($phase_lib, $triplet_lib) = add_transcript_to_phase_lib($phase_lib, $triplet_lib, $exon_struct, $strand, $max_tr_rank);
    }
    
    #Disconnect
    $dbh_ens->disconnect;
    
    return ($phase_lib, $triplet_lib);
}

#Add transcript to phase lib
sub add_transcript_to_phase_lib{
    
    #Catch
    my $phase_lib = $_[0];
    my $triplet_lib = $_[1];
    my $exon_struct = $_[2];
    my $strand = $_[3];
    my $max_tr_rank = $_[4];
    
    my $cur_phase = 0;
    my $cur_transcriptomic_pos = 1;
    my $cur_triplet = "";
    if($strand eq '1'){
        for(my $i=1;$i<=$max_tr_rank;$i++){
            my $position = $exon_struct->{$i}->{'start'};
            while($position<=$exon_struct->{$i}->{'end'}){
                $phase_lib->{$strand}->{$position}->{"phase"} = $cur_phase; #Phase lib
                #Get new triplet for phase 0
                if($cur_phase==0){
                    $cur_triplet = substr($exon_struct->{'sequence'},$cur_transcriptomic_pos-1,3);
                }
                $triplet_lib->{$strand}->{$position} = $cur_triplet;#Add triplet to triplet lib
                $phase_lib->{$strand}->{$position}->{"transcriptomic_pos"} = $cur_transcriptomic_pos; #For rel position calculation
                #Get ready for next position
                $cur_phase++;
                if($cur_phase == 3){
                    $cur_phase = 0;
                }
                $cur_transcriptomic_pos++;
                $position++;
            }
        }
        #Save the last position as the sequence length
        my $seq_length = $cur_transcriptomic_pos;
        #Then use this length to calculate normalized positions in the sequence
        for(my $i=1;$i<=$max_tr_rank;$i++){
            my $position2 = $exon_struct->{$i}->{'start'};
            while($position2<=$exon_struct->{$i}->{'end'}){
                $phase_lib->{$strand}->{$position2}->{"transcriptomic_pos"} = $phase_lib->{$strand}->{$position2}->{"transcriptomic_pos"} / $seq_length;
                $position2++;
            }
        }
    } elsif ($strand eq '-1'){
        for(my $i=1; $i<=$max_tr_rank; $i++){
            my $position = $exon_struct->{$i}->{'start'}; #strand relative
            while($position>=$exon_struct->{$i}->{'end'}){
                $phase_lib->{$strand}->{$position}->{"phase"} = $cur_phase;
                if($cur_phase==0){
                    $cur_triplet = substr($exon_struct->{'sequence'},$cur_transcriptomic_pos-1,3);
                }
                $triplet_lib->{$strand}->{$position} = $cur_triplet;
                $phase_lib->{$strand}->{$position}->{"transcriptomic_pos"} = $cur_transcriptomic_pos;
                #Set for next position
                $cur_phase++;
                if($cur_phase == 3){
                    $cur_phase = 0;
                }
                
                $cur_transcriptomic_pos++;
                $position--;
            }
        }
        #Save the last position as the sequence length
        my $seq_length = $cur_transcriptomic_pos;
        #Calculate normalized positions
        for(my $i=1; $i<=$max_tr_rank; $i++){
            my $position = $exon_struct->{$i}->{'start'}; #strand relative
            while($position>=$exon_struct->{$i}->{'end'}){
                $phase_lib->{$strand}->{$position}->{"transcriptomic_pos"} = $phase_lib->{$strand}->{$position}->{"transcriptomic_pos"} / $seq_length;
                $position--;
            }
        }
    }
    
    return $phase_lib, $triplet_lib;
}

sub get_exon_struct_transcript{
    
    #Catch
    my $dbh = $_[0];
    my $transcript = $_[1];
    my $chr = $_[2];
    
    #Init
    my $exon_struct = {};
    
    #Select transcript attributes
    my $query = "SELECT tr.start_exon_id, tr.seq_start, tr.end_exon_id, tr.seq_end, t.seq_region_strand FROM transcript as t JOIN translation as tr ON t.canonical_translation_id=tr.translation_id WHERE t.transcript_id='$transcript';";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    
    my $transcript_atts = $sth->fetchrow_hashref();
    
    #Get ranked exon structure of transcript
    $query = "SELECT et.rank, e.exon_id, e.seq_region_start, e.seq_region_end, e.phase, e.end_phase FROM exon_transcript as et JOIN exon as e ON et.exon_id=e.exon_id WHERE et.transcript_id='$transcript';";
    $sth = $dbh->prepare($query);
    $sth->execute();
    
    my $ranked_exons = $sth->fetchall_hashref('rank');
    
    #Get max rank
    $query = "SELECT MAX(rank) FROM exon_transcript WHERE transcript_id='$transcript';";
    $sth = $dbh->prepare($query);
    $sth->execute();
    my @result = $sth->fetchrow_array();
    
    my $max_rank = $result[0];
    my $tr_rank = 1;
    my $start_exon_passed = 'False';
    my $max_tr_rank=0;
    
    if($transcript_atts->{'seq_region_strand'}==1){
        for(my $rank=1;$rank<=$max_rank;$rank++){
            #Search for the translation start exon
            if($ranked_exons->{$rank}->{'exon_id'} eq $transcript_atts->{'start_exon_id'}){
                #Correct for translation start site
                $exon_struct->{$tr_rank}->{'start'} = $ranked_exons->{$rank}->{'seq_region_start'} + $transcript_atts->{'seq_start'} - 1;
                #Check if start exon is also the end exon
                if($ranked_exons->{$rank}->{'exon_id'} eq $transcript_atts->{'end_exon_id'}){
                    $exon_struct->{$tr_rank}->{'end'} = $ranked_exons->{$rank}->{'seq_region_start'} + $transcript_atts->{'seq_end'} - 1;
                    $max_tr_rank=$tr_rank;
                } else {
                    $exon_struct->{$tr_rank}->{'end'} = $ranked_exons->{$rank}->{'seq_region_end'};
                }
                $exon_struct->{$tr_rank}->{'seq'} = get_sequence($chr,$exon_struct->{$tr_rank}->{'start'},$exon_struct->{$tr_rank}->{'end'});
                #print "exon: ".$exon_struct->{$tr_rank}->{'seq'}."\n";
                #print "\tstart: ".$exon_struct->{$tr_rank}->{'start'}."\tend: ".$exon_struct->{$tr_rank}->{'end'}."\n";
                $start_exon_passed = 'True';
                $tr_rank++;
            } elsif ($start_exon_passed eq 'True'){
                #Check if in last translated exon
                if($ranked_exons->{$rank}->{'exon_id'} eq $transcript_atts->{'end_exon_id'}){
                    $exon_struct->{$tr_rank}->{'start'} = $ranked_exons->{$rank}->{'seq_region_start'};
                    $exon_struct->{$tr_rank}->{'end'} = $ranked_exons->{$rank}->{'seq_region_start'} + $transcript_atts->{'seq_end'} - 1;
                    $exon_struct->{$tr_rank}->{'seq'} = get_sequence($chr,$exon_struct->{$tr_rank}->{'start'},$exon_struct->{$tr_rank}->{'end'});
                    #print "exon: ".$exon_struct->{$tr_rank}->{'seq'}."\n";
                    #print "\tstart: ".$exon_struct->{$tr_rank}->{'start'}."\tend: ".$exon_struct->{$tr_rank}->{'end'}."\n";
                    $max_tr_rank=$tr_rank;
                    last;
                } else {
                    $exon_struct->{$tr_rank}->{'start'} = $ranked_exons->{$rank}->{'seq_region_start'};
                    $exon_struct->{$tr_rank}->{'end'} = $ranked_exons->{$rank}->{'seq_region_end'};
                    $exon_struct->{$tr_rank}->{'seq'} = get_sequence($chr,$exon_struct->{$tr_rank}->{'start'},$exon_struct->{$tr_rank}->{'end'});
                    #print "exon: ".$exon_struct->{$tr_rank}->{'seq'}."\n";
                    #print "\tstart: ".$exon_struct->{$tr_rank}->{'start'}."\tend: ".$exon_struct->{$tr_rank}->{'end'}."\n";
                    $tr_rank++;
                }
            }
        }
        $exon_struct->{'sequence'} = "";
        #Get total translated sequence by concatenation
        for(my $rank = 1; $rank<=$max_tr_rank;$rank++){
            $exon_struct->{'sequence'} = $exon_struct->{'sequence'}.$exon_struct->{$rank}->{'seq'};
        }
        #print "Max tr rank: ".$max_tr_rank."\n";
        #print "Total: ".$exon_struct->{'sequence'}."\n\n";
    } elsif($transcript_atts->{'seq_region_strand'}==-1){
        for(my $rank=1;$rank<=$max_rank;$rank++){
            #Search for the translation start exon
            if($ranked_exons->{$rank}->{'exon_id'} eq $transcript_atts->{'start_exon_id'}){
                $exon_struct->{$tr_rank}->{'start'} = $ranked_exons->{$rank}->{'seq_region_end'} - $transcript_atts->{'seq_start'} + 1;
                #Check if start exon is also the end exon
                if($ranked_exons->{$rank}->{'exon_id'} eq $transcript_atts->{'end_exon_id'}){
                    $exon_struct->{$tr_rank}->{'end'} = $ranked_exons->{$rank}->{'seq_region_end'} - $transcript_atts->{'seq_end'} + 1;
                    $max_tr_rank=$tr_rank;
                } else {
                    $exon_struct->{$tr_rank}->{'end'} = $ranked_exons->{$rank}->{'seq_region_start'};
                }
                $exon_struct->{$tr_rank}->{'seq'} = revdnacomp(get_sequence($chr,$exon_struct->{$tr_rank}->{'end'},$exon_struct->{$tr_rank}->{'start'}));
                #print "exon: ".$exon_struct->{$tr_rank}->{'seq'}."\n";
                #print "\tstart: ".$exon_struct->{$tr_rank}->{'start'}."\tend: ".$exon_struct->{$tr_rank}->{'end'}."\n";
                $start_exon_passed = 'True';
                $tr_rank++;
            } elsif ($start_exon_passed eq 'True'){
                #Check if in last translated exon
                if($ranked_exons->{$rank}->{'exon_id'} eq $transcript_atts->{'end_exon_id'}){
                    $exon_struct->{$tr_rank}->{'start'} = $ranked_exons->{$rank}->{'seq_region_end'};
                    $exon_struct->{$tr_rank}->{'end'} = $ranked_exons->{$rank}->{'seq_region_end'} - $transcript_atts->{'seq_end'} + 1;
                    $exon_struct->{$tr_rank}->{'seq'} = revdnacomp(get_sequence($chr,$exon_struct->{$tr_rank}->{'end'},$exon_struct->{$tr_rank}->{'start'}));
                    #print "exon: ".$exon_struct->{$tr_rank}->{'seq'}."\n";
                    #print "\tstart: ".$exon_struct->{$tr_rank}->{'start'}."\tend: ".$exon_struct->{$tr_rank}->{'end'}."\n";
                    $max_tr_rank=$tr_rank;
                    last;
                } else {
                    $exon_struct->{$tr_rank}->{'start'} = $ranked_exons->{$rank}->{'seq_region_end'};
                    $exon_struct->{$tr_rank}->{'end'} = $ranked_exons->{$rank}->{'seq_region_start'};
                    $exon_struct->{$tr_rank}->{'seq'} = revdnacomp(get_sequence($chr,$exon_struct->{$tr_rank}->{'end'},$exon_struct->{$tr_rank}->{'start'}));
                    #print "exon: ".$exon_struct->{$tr_rank}->{'seq'}."\n";
                    #print "\tstart: ".$exon_struct->{$tr_rank}->{'start'}."\tend: ".$exon_struct->{$tr_rank}->{'end'}."\n";
                    $tr_rank++;
                }
            }
        }
        $exon_struct->{'sequence'} = "";
        #Get total translated sequence by concatenation
        for(my $rank = 1; $rank<=$max_tr_rank;$rank++){
            $exon_struct->{'sequence'} = $exon_struct->{'sequence'}.$exon_struct->{$rank}->{'seq'};
        }
        #print "Total: ".$exon_struct->{'sequence'}."\n\n";
    }
    
    my $strand = $transcript_atts->{'seq_region_strand'};
    
    return $exon_struct, $strand, $max_tr_rank;
}

#Get transcripts
sub get_can_transcripts{
    
    #Catch
    my $dbh = $_[0];
    my $seq_region_id = $_[1];
    
    #Init
    my $transcripts = [];
    
    #Take all protein-coding genes and take the canonical transcript ID's from them for the chromosome (seq_region_id)
    my $query = "SELECT t.transcript_id FROM transcript as t JOIN gene as g ON g.canonical_transcript_id=t.transcript_id WHERE g.seq_region_id='$seq_region_id' AND g.biotype='protein_coding';";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    
    while (my @row = $sth->fetchrow_array()){
        push @{$transcripts}, $row[0];
    }
    
    return $transcripts;
}

#Parse RIBO_CIGARS to obtain offset,genomic read mapping length and total intronic length before offset is reached
sub parse_RIBO_CIGAR {
    
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

### SPLIT SAM PER CHR ###
sub split_SAM_per_chr {
    
    # Catch
    my %chr_sizes = %{$_[0]};
    my $work_dir = $_[1];
    my $sam = $_[2];
    my $treated = $_[3];
    
    my @splitsam = split(/\//, $sam );
    my $samFileName = $splitsam[$#splitsam];
    @splitsam = split(/\./,$samFileName);
    $samFileName = $splitsam[0];
    
    my ($chr,@mapping_store,$file_in_loc,$file_in,$file_out);
    
    #Remove eventual existing samfiles
    system("rm -f ".$TMP."/mappingqc_".$treated."/".$samFileName."_*.sam");   # Delete existing
    
    # Touch per chr
    foreach $chr (keys %chr_sizes){
        system("touch ".$TMP."/mappingqc_".$treated."/".$samFileName."_".$chr.".sam");   # Touch new
    }
    
    ## Split files into chromosomes
    
    # Open
    open (I,"<".$sam) || die "Cannot open ".$sam." file\n";
    
    
    #For unsorted SAM file (genomic location)
    my $prev_chr="0";
    
    while(my $line=<I>){
        
        #Skip annotation lines
        if ($line =~ m/^@/) { next; }
        
        #Process alignment line
        @mapping_store = split(/\t/,$line);
        $chr = $mapping_store[2];
        
        #Keep all mappings, also MultipleMapping locations are available (alternative to pseudogenes mapping) GM:07-10-2013
        #Note that we only retain the up until <16 multiple locations (to avoid including TopHat2 peak @ 16)
        #For STAR the maxMultiMap is not included in the output (e.g. if maxmultimap is set to 16, up untill 15 is included)
        #For TopHat2 the maxMultiMap is included in the output (e.g. if maxmultimap is set to 16, up untill 16 is included)
        #In order to have the same actual limit, the maxmultimap is discarded (see record below)
        next unless ( $line !~ m/NH:i:$maxmultimap/ );
        
        # Write off
        if ($prev_chr ne $chr) {
            if ($prev_chr ne "0") { close(A);}
            $file_out = $TMP."/mappingqc_".$treated."/".$samFileName;
            open (A,">>".$file_out."_".$chr.".sam") || die "Cannot open the sep file";
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

## Get coord system id ##
sub get_coord_system_id{
    # Catch
    my $db_ensembl = $_[0];
    my $assembly = $_[1];
    
    my $user = "";
    my $pw = "";
    
    # Connect to ensembl sqlite database
    my $dbh  = DBI->connect('DBI:SQLite:'.$db_ensembl,$user,$pw,
    { RaiseError => 1},) || die "Database connection not made: $DBI::errstr";
    
    # Get correct coord_system_id
    my $query = "SELECT coord_system_id FROM coord_system WHERE name = 'chromosome' AND version = '$assembly'";
    my $execute = $dbh->prepare($query);
    $execute->execute();
    
    my $coord_system_id;
    while(my @result = $execute->fetchrow_array()){
        $coord_system_id = $result[0];
    }
    
    $execute->finish();
    
    # Disconnect
    $dbh->disconnect();
    
    # Return
    return($coord_system_id);
} # Close sub

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

### Get sequence from chromosomal sequence fasta files ###
sub get_sequence{
    #Catch
    my $chr = $_[0];
    my $start = $_[1];
    my $end = $_[2];
    
    #Init
    my $seq;
    
    #Open chromosomal sequence fasta file (BINARY)
    open(IN, "< ".$BIN_chrom_dir."/".$chr.".fa");
    binmode(IN);
    
    #translate start and end
    my $offset = $start - 1;
    my $length = $end - $start + 1;
    
    #Define reading start position
    seek(IN, $offset, 0);
    #Read in sequence
    read(IN, $seq, $length);
    
    close(IN);
    
    return $seq;
}

### get reverse complement sequence ###
sub revdnacomp {
    #Catch
    my $dna = $_[0];
    
    my $revcomp = reverse($dna);
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}

## Get arguments out of arguments table
sub get_ARG_vars{
    # Catch
    my $db_ribo = $_[0];
    my $user = $_[1];
    my $pw = $_[2];
    
    my ($query,$sth);
    
    # Connect to db
    my $dbh_results = DBI->connect('DBI:SQLite:'.$db_ribo,$user,$pw,
    { RaiseError => 1},) || die "Database connection not made: $DBI::errstr";
    
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
    
    $dbh_results -> disconnect();
    
    # Return ARG variables
    return($species,$version,$IGENOMES_ROOT);
} # Close sub

## GET seq_region_id ##
sub get_seq_region_id{
    
    #Catch
    my $dbh = $_[0];
    my $chr = $_[1];
    my $coord_system_id = $_[2];
    
    #Init
    my $seq_region_id;
    
    my $query = "SELECT seq_region_id FROM seq_region WHERE coord_system_id = '$coord_system_id' AND name = '$chr';";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    $seq_region_id = $sth->fetch()->[0];
    $sth->finish();
    
    return $seq_region_id;
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

##Write offsets to tmp folder for html file
sub offsets_to_csv {
    
    #Catch
    my $offset_hash = $_[0];
    my $treated = $_[1];
    my $TMP = $_[2];
    
    my $outfile = $TMP."/mappingqc_".$treated."/mappingqc_offsets_".$treated.".csv";
    open(my $fw, '>', $outfile);
    
    for(my $key=$offset_hash->{'min'}; $key<=$offset_hash->{'max'}; $key++){
        my $length = $key;
        my $offset = $offset_hash->{$key};
        my $line = $length.",".$offset."\n";
        print $fw $line;
    }
    
    close($fw);
}
