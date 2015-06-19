#!/usr/bin/perl -w
use strict;
use warnings;
use DBI;
use DBD::SQLite;
use Data::Dumper;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Getopt::Long;
use v5.10;
use Parallel::ForkManager;
use Storable 'dclone';
use Cwd;

##############
##Command-line
# ./PhyloCSF_sORFs_input2.pl --sqliteDB SQLite/results.db --cores 25 --maf_root /storage2/MAF/ --tis_ids 1

# get the command line arguments
my ($resultDB,$cores,$maf_root,$tis_ids,$out_sqlite,$tmpfolder,$work_dir);

GetOptions(
"sqliteDB=s"=>\$resultDB,                  # The sqlite DB holding all RIBO-pipeline results,                      mandatory argument
"cores=i"=>\$cores,                         # Number of cores to use for Bowtie Mapping,                            mandatory argument
"maf_root=s" =>\$maf_root,                  # MAF ROOT FOLDER                                                       mandatory argument
"tis_ids=s"  =>\$tis_ids,                   # list of analysis ids                                                  mandatory argument
"out_sqlite=s"=>\$out_sqlite                # The sqlite DB holding all the output                                  mandatory argument
);

my $CWD             = getcwd;
my $TMP             = ($ENV{'TMP'}) ? $ENV{'TMP'} : ($tmpfolder) ? $tmpfolder : "$CWD/tmp"; # (1) get the TMP environment variable, (2) get the $tmpfolder variable, (3) get

print "The following tmpfolder is used                          : $TMP\n";
print "The following MAF folder is used                         : $maf_root\n";

#Check if tmpfolder exists, if not create it...
if (!-d "$TMP") {
    system ("mkdir ". $TMP);
}

#comment on these
if ($resultDB){
    print "SqliteDB used is                                         : $resultDB\n";
} else {
    die "\nDon't forget to pass the SQLite DB using the --sqliteDB argument!\n\n";
}
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
if ($cores){
    print "Number of cores to use for assembling                    : $cores\n";
} else {
    die "\nDon't forget to pass number of cores to use for mapping using the --cores or -c argument!\n\n";
}

# Create output files for command line script
if (!defined($out_sqlite))       {$out_sqlite          = $resultDB;}

# DB settings
# Sqlite Riboseq
my $db_results  = $resultDB;
my $dsn_results = "DBI:SQLite:dbname=$db_results";
my $us_results  = "";
my $pw_results  = "";

#Get the input variables (Check nog welke eruit mogen)
my $dbh_results = dbh($dsn_results,$us_results,$pw_results);
my ($ensemblversion,$species,$ens_db,$igenomes_root) = get_arguments($dsn_results,$us_results,$pw_results);

print "The following Ensembl db folder is used                  : $ens_db\n";
print "The igenomes_root folder used is                         : $igenomes_root\n";

# Sqlite Ensembl
my $db_ENS  = $ens_db;
my $dsn_ENS = "DBI:SQLite:dbname=$db_ENS";
my $us_ENS  = "";
my $pw_ENS  = "";

#Conversion for species terminology
my $spec = ($species eq "mouse") ? "Mus_musculus" : ($species eq "human") ? "Homo_sapiens" : ($species eq "arabidopsis") ? "Arabidopsis_thaliana" : ($species eq "fruitfly") ? "Drosophila_melanogaster" : "";
my $spec_short = ($species eq "mouse") ? "mmu" : ($species eq "human") ? "hsa" : ($species eq "arabidopsis") ? "ath" : ($species eq "fruitfly") ? "dme" : "";
#Old mouse assembly = NCBIM37, new one is GRCm38
my $assembly = ($species eq "mouse" && $ensemblversion >= 70 ) ? "GRCm38"
: ($species eq "mouse" && $ensemblversion < 70 ) ? "NCBIM37"
: ($species eq "human") ? "GRCh37"
: ($species eq "arabidopsis") ? "TAIR10"
: ($species eq "fruitfly") ? "BDGP5" : "";

# my $chrom_file = $igenomes_root."/".$spec."/Ensembl/".$assembly."/Annotation/Genes/ChromInfo.txt";

my $BIN_chrom_dir = $igenomes_root."/".$spec."/Ensembl/".$assembly."/Sequence/Chromosomes_BIN";

############################################
##                                        ##
## Creating MAF based PhyloCSF input data ##
##                                        ##
############################################

# Start time
my $start = time;

# Get the analysis_id that corresponds to the TIS-calling input parameters
my $idsref = get_analysis_ids($dbh_results,$tis_ids);

## Get chromosome sizes

print "Getting chromosome sizes  ...\n";

my $chr_sizes = get_chr_sizes();

print "\nGet chromosomes... \n";
## Get chromosomes based on seq_region_id ##
my $chrs = get_chrs($dsn_ENS,$us_ENS,$pw_ENS,get_chr_sizes(),$assembly,$species);

my $analysis_id;
#Loop over all selected analysis_ids
print "Creating MAF based PhyloCSF output...\n";

foreach $analysis_id (@$idsref) {
    
    print "Processing analysis_id $analysis_id ...\n";
    create_maf_sORF_alignments($chrs,$TMP,$analysis_id);
    
    ## Store in DB
    print "   Storing all sORFs and alignments in DB \n";
    store_in_db($dsn_results,$us_results,$pw_results,$analysis_id,$work_dir,$chrs,$TMP,$species);
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

## Create MAF based input for PhyloCSF

sub create_maf_sORF_alignments{

    #Catch
    my $chrs        =   $_[0];
    my $tmp         =   $_[1];
    my $analysis_id =   $_[2];
    
    # Init multi core
    my $pm = new Parallel::ForkManager($cores);
    
    ## Loop over all chromosomes
    foreach my $chr (sort keys %{$chrs}){
        
        ### Start parallel process
        $pm->start and next;
        
        ### DBH per process
        my $dbh = dbh($dsn_results,$us_results,$pw_results);
        
        ### Output fasta and db_csv file per process
        open TMP_db, "+>>".$TMP."/".$analysis_id."_".$chr."_sORF_maf.csv" or die $!;
        
        ### Get sORFs per chromosome from transcript table
        my $sorfs = get_sorfs_per_chromosome($dbh,$analysis_id,$chr);
        
        ### Match maf blocks to sORFs
        $sorfs = match_maf_blocks_to_sorfs($maf_root,$chr,$spec,$sorfs);

        ### Get MAF blocks to Hash
        # my $MAF = MafToHash($chr,$maf_root,$spec);
        
        # print Dumper($sorfs);
        
        ### Run over sORFs, create alignments
        foreach my $sorf_id (sort keys %{$sorfs}){
        
            ## Get Maf blocks overlapping sorf
            #my $sORF_MAF = SorfMafBlocks($sorfs,$sorf_id,$MAF);

            #Check for empty hash maf blocks
            if (! exists ($sorfs->{$sorf_id}{'aln'})){next;}
            
            ## Form sORF MAF alignment
            my $sORF_align = sORF_align($sorfs,$sorf_id,$species);
            
            ## Store in CSV
            store_in_csv($sORF_align,$species,$sorfs,$sorf_id);
        }
    
        #Finish chromosome analysis
        print "     * Finished analyzing chromosome ".$chr."\n";
        
        ### Finish child
        $dbh->disconnect();
        close TMP_db;
        $pm->finish;
    }
    
    #Wait all Children
    $pm->wait_all_children();
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
    my $species     =   $_[7];
    
    # Init
    my $dbh     =   dbh($dsn,$us,$pw);
    my $table   =   "TIS_sORFs_".$id."_transcripts_maf";
    my $query;
    
    print "species: ".$species."\n";
    
    # Create table
    if ($species eq "fruitfly"){
        $query = "CREATE TABLE IF NOT EXISTS `".$table."` (
        `sorf_id` varchar(128) NOT NULL default '',
        `chr` char(50) NOT NULL default '',
        `strand` int(2) NOT NULL default '',
        `sorf_begin` int(10) NOT NULL default '',
        `sorf_end` int(10) NOT NULL default '',
        `dmel` TEXT NOT NULL default '',
        `dsim` TEXT NOT NULL default '',
        `dsec` TEXT NOT NULL default '',
        `dyak` TEXT NOT NULL default '',
        `dere` TEXT NOT NULL default '',
        `dana` TEXT NOT NULL default '',
        `dpse` TEXT NOT NULL default '',
        `dper` TEXT NOT NULL default '',
        `dwil` TEXT NOT NULL default '',
        `dvir` TEXT NOT NULL default '',
        `dmoj` TEXT NOT NULL default '',
        `dgri` TEXT NOT NULL default '' )"  ;
    }elsif ($species eq "human"){
        $query = "CREATE TABLE IF NOT EXISTS `".$table."` (
        `sorf_id` varchar(128) NOT NULL default '',
        `chr` char(50) NOT NULL default '',
        `strand` int(2) NOT NULL default '',
        `sorf_begin` int(10) NOT NULL default '',
        `sorf_end` int(10) NOT NULL default '',
        `Human` TEXT NOT NULL default '',
        `Chimp` TEXT NOT NULL default '',
        `Rhesus` TEXT NOT NULL default '',
        `Bushbaby` TEXT NOT NULL default '',
        `TreeShrew` TEXT NOT NULL default '',
        `Mouse` TEXT NOT NULL default '',
        `Rat` TEXT NOT NULL default '',
        `Guinea_Pig` TEXT NOT NULL default '',
        `Squirrel` TEXT NOT NULL default '',
        `Rabbit` TEXT NOT NULL default '',
        `Pika` TEXT NOT NULL default '',
        `Alpaca` TEXT NOT NULL default '',
        `Dolphin` TEXT NOT NULL default '',
        `Cow` TEXT NOT NULL default '',
        `Horse` TEXT NOT NULL default '',
        `Cat` TEXT NOT NULL default '',
        `Dog` TEXT NOT NULL default '',
        `Microbat` TEXT NOT NULL default '',
        `Megabat` TEXT NOT NULL default '',
        `Hedgehog` TEXT NOT NULL default '',
        `Shrew` TEXT NOT NULL default '',
        `Elephant` TEXT NOT NULL default '',
        `Tenrec` TEXT NOT NULL default '',
        `Armadillo` TEXT NOT NULL default '')"  ;
    }elsif ($species eq "mouse"){
        $query = "CREATE TABLE IF NOT EXISTS `".$table."` (
        `sorf_id` varchar(128) NOT NULL default '',
        `chr` char(50) NOT NULL default '',
        `strand` int(2) NOT NULL default '',
        `sorf_begin` int(10) NOT NULL default '',
        `sorf_end` int(10) NOT NULL default '',
        `Mouse` TEXT NOT NULL default '',
        `Guinea_Pig` TEXT NOT NULL default '',
        `Kangaroo_rat` TEXT NOT NULL default '',
        `Pika` TEXT NOT NULL default '',
        `Rabbit` TEXT NOT NULL default '',
        `Rat` TEXT NOT NULL default '',
        `Squirrel` TEXT NOT NULL default '',
        `TreeShrew` TEXT NOT NULL default '',
        `Human` TEXT NOT NULL default '',
        `Mouse_lemur` TEXT NOT NULL default '',
        `Bushbaby` TEXT NOT NULL default '',
        `Chimp` TEXT NOT NULL default '',
        `Rhesus` TEXT NOT NULL default '',
        `Tarsier` TEXT NOT NULL default '',
        `Cow` TEXT NOT NULL default '',
        `Dog` TEXT NOT NULL default '',
        `Sloth` TEXT NOT NULL default '',
        `Armadillo` TEXT NOT NULL default '',
        `Tenrec` TEXT NOT NULL default '',
        `Horse` TEXT NOT NULL default '',
        `Hedgehog` TEXT NOT NULL default '',
        `Cat` TEXT NOT NULL default '',
        `Elephant` TEXT NOT NULL default '',
        `Microbat` TEXT NOT NULL default '',
        `Rock_hyrax` TEXT NOT NULL default '',
        `Megabat` TEXT NOT NULL default '',
        `Shrew` TEXT NOT NULL default '',
        `Dolphin` TEXT NOT NULL default '',
        `Alpaca` TEXT NOT NULL default '' )"  ;
    }
    $dbh->do($query);
    
    # Store
    foreach my $chr (sort keys %{$chrs}){
        system("sqlite3 -separator , ".$resultDB." \".import ".$TMP."/".$id."_".$chr."_sORF_maf.csv ".$table."\"")== 0 or die "system failed: $?";
    }
    
    #Disconnect dbh
    $dbh->disconnect();
    
    # Unlink tmp csv files
    foreach my $chr (sort keys %{$chrs}){
        unlink $TMP."/".$id."_".$chr."_sORF_maf.csv";
    }
    
}

### Store in CSV ###

sub store_in_csv{

    #Catch
    my $sORF_align  =   $_[0];
    my $species     =   $_[1];
    my $sorfs       =   $_[2];
    my $sorf_id     =   $_[3];
    
    #Init
    my $length      =   ($species eq "fruitfly") ? length($sORF_align->{'dm3'}{'DNA_seq'}) : ($species eq "human") ? length($sORF_align->{'hg19'}{'DNA_seq'}) : ($species eq "mouse") ? length($sORF_align->{'mm10'}{'DNA_seq'}) : 0;
    my $DNA_seq;
    
    my $PhyloCSF_species;
    if ($species eq "fruitfly"){
        $PhyloCSF_species =  ['dm3','droSim1','droSec1','droYak2','droEre2','droAna3','dp4','droPer1','droWil1','droVir3','droMoj3','droGri2'];
    }elsif ($species eq "human"){
        $PhyloCSF_species =  ['hg19','panTro4','rheMac3','otoGar3','tupChi1','speTri2','mm10','rn5','cavPor3','oryCun2','ochPri3','vicPac2','turTru2','bosTau7','equCab2','felCat5','canFam3','pteVam1','myoLuc2','eriEur2','sorAra2','loxAfr3','echTel2','dasNov2'];
    }elsif ($species eq "mouse"){
        $PhyloCSF_species =  ['mm10','cavPor3','dipOrd1','ochPri2','oryCun2','rn5','speTri2','tupBel1','hg19','micMur1','otoGar3','panTro4','rheMac3','tarSyr1','bosTau7','canFam3','choHof1','dasNov3','echTel1','equCab2','eriEur1','felCat5','loxAfr3','myoLuc2','proCap1','pteVam1','sorAra1','turTru2','vicPac1'];
    }
    
    #Create tmp_csv
    my $tmp_csv     =   $sorf_id.",".$sorfs->{$sorf_id}{'chr'}.",".$sorfs->{$sorf_id}{'strand'}.",".$sorfs->{$sorf_id}{'sorf_begin'}.",".$sorfs->{$sorf_id}{'sorf_end'}.",";

    foreach (@{$PhyloCSF_species}){
        if (exists $sORF_align->{$_}{'DNA_seq'}){
            if (length($sORF_align->{$_}{'DNA_seq'}) == $length){
                if ($sorfs->{$sorf_id}{'strand'} eq "1"){
                    $tmp_csv    .=  $sORF_align->{$_}{'DNA_seq'}.",";
                }elsif ($sorfs->{$sorf_id}{'strand'} eq "-1"){
                    $DNA_seq    =   reverse($sORF_align->{$_}{'DNA_seq'});
                    $DNA_seq    =~  tr/ACGTacgt/TGCATGCA/;
                    $tmp_csv    .=  $DNA_seq.",";
                }
            }else{
                $tmp_csv    .=  ",";
            }
        }else{
            $tmp_csv    .=  ",";
        }
    }
    
    # Get rid of latest ","
    chop($tmp_csv);
    
    # Write to CSV
    print TMP_db $tmp_csv."\n";
    
}

### Form DNA alignment for sORF (forward strand based) ###

sub sORF_align {
    
    # Catch
    my $sorfs       =   $_[0];
    my $sorf_id     =   $_[1];
    my $species     =   $_[2];
    
    # Init
    my $sorf_begin  =   $sorfs->{$sorf_id}{'sorf_begin'};
    my $sorf_end    =   $sorfs->{$sorf_id}{'sorf_end'};
    my $sORF_MAF    =   $sorfs->{$sorf_id}->{'aln'};
    
    # print Dumper($sORF_MAF);
    
    my %alignment   =   ();
    
    if ($species eq "fruitfly") {
        %alignment = (
        'dm3'       =>  {'PhyloCSF_name' => 'dmel',},
        'droSim1'   =>  {'PhyloCSF_name' => 'dsim',},
        'droSec1'   =>  {'PhyloCSF_name' => 'dsec',},
        'droYak2'   =>  {'PhyloCSF_name' => 'dyak',},
        'droEre2'   =>  {'PhyloCSF_name' => 'dere',},
        'droAna3'   =>  {'PhyloCSF_name' => 'dana',},
        'dp4'       =>  {'PhyloCSF_name' => 'dpse',},
        'droPer1'   =>  {'PhyloCSF_name' => 'dper',},
        'droWil1'   =>  {'PhyloCSF_name' => 'dwil',},
        'droVir3'   =>  {'PhyloCSF_name' => 'dvir',},
        'droMoj3'   =>  {'PhyloCSF_name' => 'dmoj',},
        'droGri2'   =>  {'PhyloCSF_name' => 'dgri',},
        'anoGam1'   =>  {'' => '',},
        'apiMel3'   =>  {'' => '',},
        'triCas2'   =>  {'' => '',},
        );
    }elsif ($species eq "human"){
        %alignment = (
        'hg19'          =>  {'PhyloCSF_name' => 'Human',},
        'panTro4'       =>  {'PhyloCSF_name' => 'Chimp',},
        'gorGor3'       =>  {'' => '',},
        'ponAbe2'       =>  {'' => '',},
        'nomLeu3'       =>  {'' => '',},
        'rheMac3'       =>  {'PhyloCSF_name' => 'Rhesus',},
        'macFas5'       =>  {'' => '',},
        'papHam1'       =>  {'' => '',},
        'chlSab1'       =>  {'' => '',},
        'calJac3'       =>  {'' => '',},
        'saiBol1'       =>  {'' => '',},
        'otoGar3'       =>  {'PhyloCSF_name' => 'Bushbaby',},
        'tupChi1'       =>  {'PhyloCSF_name' => 'TreeShrew',},
        'speTri2'       =>  {'PhyloCSF_name' => 'Squirrel',},
        'jacJac1'       =>  {'' => '',},
        'micOch1'       =>  {'' => '',},
        'criGri1'       =>  {'' => '',},
        'mesAur1'       =>  {'' => '',},
        'mm10'          =>  {'PhyloCSF_name' => 'Mouse',},
        'rn5'           =>  {'PhyloCSF_name' => 'Rat',},
        'hetGla2'       =>  {'' => '',},
        'cavPor3'       =>  {'PhyloCSF_name' => 'Guinea_Pig',},
        'chiLan1'       =>  {'' => '',},
        'octDeg1'       =>  {'' => '',},
        'oryCun2'       =>  {'PhyloCSF_name' => 'Rabbit',},
        'ochPri3'       =>  {'PhyloCSF_name' => 'Pika',},
        'susScr3'       =>  {'' => '',},
        'vicPac2'       =>  {'PhyloCSF_name' => 'Alpaca',},
        'camFer1'       =>  {'' => '',},
        'turTru2'       =>  {'PhyloCSF_name' => 'Dolphin',},
        'orcOrc1'       =>  {'' => '',},
        'panHod1'       =>  {'' => '',},
        'bosTau7'       =>  {'PhyloCSF_name' => 'Cow',},
        'oviAri3'       =>  {'' => '',},
        'capHir1'       =>  {'' => '',},
        'equCab2'       =>  {'PhyloCSF_name' => 'Horse',},
        'cerSim1'       =>  {'' => '',},
        'felCat5'       =>  {'PhyloCSF_name' => 'Cat',},
        'canFam3'       =>  {'PhyloCSF_name' => 'Dog',},
        'musFur1'       =>  {'' => '',},
        'ailMel1'       =>  {'' => '',},
        'odoRosDiv1'    =>  {'' => '',},
        'lepWed1'       =>  {'' => '',},
        'pteAle1'       =>  {'' => '',},
        'pteVam1'       =>  {'PhyloCSF_name' => 'Megabat',},
        'myoDav1'       =>  {'' => '',},
        'myoLuc2'       =>  {'PhyloCSF_name' => 'Microbat',},
        'eptFus1'       =>  {'' => '',},
        'eriEur2'       =>  {'PhyloCSF_name' => 'Hedgehog',},
        'sorAra2'       =>  {'PhyloCSF_name' => 'Shrew',},
        'conCri1'       =>  {'' => '',},
        'loxAfr3'       =>  {'PhyloCSF_name' => 'Elephant',},
        'eleEdw1'       =>  {'' => '',},
        'triMan1'       =>  {'' => '',},
        'chrAsi1'       =>  {'' => '',},
        'echTel2'       =>  {'PhyloCSF_name' => 'Tenrec',},
        'oryAfe1'       =>  {'' => '',},
        'dasNov3'       =>  {'PhyloCSF_name' => 'Armadillo',},
        'monDom5'       =>  {'' => '',},
        'sarHar1'       =>  {'' => '',},
        'macEug2'       =>  {'' => '',},
        'ornAna1'       =>  {'' => '',},
        'falChe1'       =>  {'' => '',},
        'falPer1'       =>  {'' => '',},
        'ficAlb2'       =>  {'' => '',},
        'zonAlb1'       =>  {'' => '',},
        'geoFor1'       =>  {'' => '',},
        'taeGut2'       =>  {'' => '',},
        'pseHum1'       =>  {'' => '',},
        'melUnd1'       =>  {'' => '',},
        'amaVit1'       =>  {'' => '',},
        'araMac1'       =>  {'' => '',},
        'colLiv1'       =>  {'' => '',},
        'anaPla1'       =>  {'' => '',},
        'galCal4'       =>  {'' => '',},
        'melGal1'       =>  {'' => '',},
        'allMis1'       =>  {'' => '',},
        'cheMyd1'       =>  {'' => '',},
        'chrPic1'       =>  {'' => '',},
        'pelSin1'       =>  {'' => '',},
        'apaSpi1'       =>  {'' => '',},
        'anoCar2'       =>  {'' => '',},
        'xenTro7'       =>  {'' => '',},
        'latCha1'       =>  {'' => '',},
        'tetNig2'       =>  {'' => '',},
        'fr3'           =>  {'' => '',},
        'takFla1'       =>  {'' => '',},
        'oreNil2'       =>  {'' => '',},
        'neoBri1'       =>  {'' => '',},
        'hapBur1'       =>  {'' => '',},
        'mayZeb1'       =>  {'' => '',},
        'punNye1'       =>  {'' => '',},
        'oryLat2'       =>  {'' => '',},
        'xipMac1'       =>  {'' => '',},
        'gasAcu1'       =>  {'' => '',},
        'gadMor1'       =>  {'' => '',},
        'danRer7'       =>  {'' => '',},
        'astMex1'       =>  {'' => '',},
        'lepOcu1'       =>  {'' => '',},
        'petMar2'       =>  {'' => '',},
         );
    }elsif ($species eq "mouse"){
        %alignment = (
        'mm10'          =>  {'PhyloCSF_name' => 'Mouse',},
        'cavPor3'       =>  {'PhyloCSF_name' => 'Guinea_Pig',},
        'dipOrd1'       =>  {'PhyloCSF_name' => 'Kangaroo_rat',},
        'hetGla2'       =>  {'' => '',},
        'ochPri2'       =>  {'PhyloCSF_name' => 'Pika',},
        'oryCun2'       =>  {'PhyloCSF_name' => 'Rabbit',},
        'rn5'           =>  {'PhyloCSF_name' => 'Rat',},
        'speTri2'       =>  {'PhyloCSF_name' => 'Squirrel',},
        'tupBel1'       =>  {'PhyloCSF_name' => 'TreeShrew',},
        'calJac3'       =>  {'' => '',},
        'gorGor3'       =>  {'' => '',},
        'hg19'          =>  {'PhyloCSF_name' => 'Human',},
        'micMur1'       =>  {'PhyloCSF_name' => 'Mouse_lemur',},
        'nomLeu2'       =>  {'' => '',},
        'otoGar3'       =>  {'PhyloCSF_name' => 'Bushbaby',},
        'panTro4'       =>  {'PhyloCSF_name' => 'Chimp',},
        'papHam1'       =>  {'' => '',},
        'ponAbe2'       =>  {'' => '',},
        'rheMac3'       =>  {'PhyloCSF_name' => 'Rhesus',},
        'saiBol1'       =>  {'' => '',},
        'tarSyr1'       =>  {'PhyloCSF_name' => 'Tarsier',},
        'ailMel1'       =>  {'' => '',},
        'bosTau7'       =>  {'PhyloCSF_name' => 'Cow',},
        'canFam3'       =>  {'PhyloCSF_name' => 'Dog',},
        'choHof1'       =>  {'PhyloCSF_name' => 'Sloth',},
        'dasNov3'       =>  {'PhyloCSF_name' => 'Armadillo',},
        'echTel1'       =>  {'PhyloCSF_name' => 'Tenrec',},
        'equCab2'       =>  {'PhyloCSF_name' => 'Horse',},
        'eriEur1'       =>  {'PhyloCSF_name' => 'Hedgehog',},
        'felCat5'       =>  {'PhyloCSF_name' => 'Cat',},
        'loxAfr3'       =>  {'PhyloCSF_name' => 'Elephant',},
        'myoLuc2'       =>  {'PhyloCSF_name' => 'Microbat',},
        'oviAri1'       =>  {'' => '',},
        'proCap1'       =>  {'PhyloCSF_name' => 'Rock_hyrax',},
        'pteVam1'       =>  {'PhyloCSF_name' => 'Megabat',},
        'sorAra1'       =>  {'PhyloCSF_name' => 'Shrew',},
        'susScr3'       =>  {'' => '',},
        'triMan1'       =>  {'' => '',},
        'turTru2'       =>  {'PhyloCSF_name' => 'Dolphin',},
        'vicPac1'       =>  {'PhyloCSF_name' => 'Alpaca',},
        'anoCar2'       =>  {'' => '',},
        'chrPic1'       =>  {'' => '',},
        'danRer7'       =>  {'' => '',},
        'fr3'           =>  {'' => '',},
        'gadMor1'       =>  {'' => '',},
        'galGal4'       =>  {'' => '',},
        'gasAcu1'       =>  {'' => '',},
        'latCha1'       =>  {'' => '',},
        'macEug2'       =>  {'' => '',},
        'melGal1'       =>  {'' => '',},
        'melUnd1'       =>  {'' => '',},
        'monDom5'       =>  {'' => '',},
        'oreNil2'       =>  {'' => '',},
        'ornAna1'       =>  {'' => '',},
        'oryLat2'       =>  {'' => '',},
        'petMar1'       =>  {'' => '',},
        'sarHar1'       =>  {'' => '',},
        'taeGut1'       =>  {'' => '',},
        'tetNig2'       =>  {'' => '',},
        'xenTro3'       =>  {'' => '',},
         );
    }
    
    # Sort maf_blocks by ASC begin position before parsing
    foreach my $key (sort { $sORF_MAF->{$a}{'begin'} <=> $sORF_MAF->{$b}{'begin'} } keys %{$sORF_MAF}){
        
        # Parsing data on sequence block
        my $aln             =   $sORF_MAF->{$key}{'mafblock'};
        my $num_sequences   =   $aln->num_sequences;
        my $seq             =   $aln->get_seq_by_pos(1);
        my $sequence        =   $seq->seq();
        my @char            =   split(//, $sequence);
        my $a               =   @char;
        my $begin           =   $seq->{start};
        my $end             =   $seq->{end};
        
        # Alignment block spans transcript begin position
        if ($begin <= $sorf_begin && $sorf_begin <= $end){
            my $offset = $sorf_begin - $begin;
            my $length = $sorf_end - $sorf_begin;
            my $offset_and_length = $offset + $length +1;
            my $count = 0;
            my $count_align = 0;
            my $alignment_begin_position = 0;
            my $alignment_end_position = 0;
            
            # Transcript begin position != alignment begin position because of multiple '-'
            foreach my $base (@char){
                
                # Get alignment position (bases and -'s)
                $count_align++;
                if ($base =~ /[ACTGactg]/){
                    
                    # Get alignment postion (only bases)
                    $count++;
                    if($count eq $offset){
                        $alignment_begin_position = $count_align;
                    }
                    elsif($count eq $offset_and_length){
                        $alignment_end_position = $count_align;
                    }
                }
            }
            
            # Transcript lays within one aligment block
            if ($begin <= $sorf_end && $sorf_end <= $end){
                
                # Get length and sequences of all species in aligment block and put in Hash
                for (my $i=1; $i<$num_sequences+1; $i++){
                    
                    my $alignment_sequence              =   $aln->get_seq_by_pos($i)->seq;
                    my $alignment_length                =   $alignment_end_position - $alignment_begin_position;
                    $alignment_sequence                 =   substr ($alignment_sequence,$alignment_begin_position,$alignment_length);
                    $alignment_sequence                 =~  tr/acgt/ACGT/;
                    my $seq_id                          =   $aln->get_seq_by_pos($i);
                    $$seq_id{ 'display_id' }            =~  m/^(\S*)\.\S*/;
                    my $seq_key                         =   $1;
                    
                    $alignment{$seq_key}->{'DNA_seq'}   =   $alignment_sequence;
                    
                    
                }
            }
            # Transcript spans over more than one alignment block
            elsif ($begin <= $sorf_end && $end <= $sorf_end){
                for (my $i=1; $i<$num_sequences+1; $i++){
                    
                    my $alignment_sequence              =   $aln->get_seq_by_pos($i)->seq;
                    my $alignment_length                =   $end - $begin - $alignment_begin_position;
                    $alignment_sequence                 =   substr ($alignment_sequence,$alignment_begin_position);
                    $alignment_sequence                 =~  tr/acgt/ACGT/;
                    my $seq_id                          =   $aln->get_seq_by_pos($i);
                    $$seq_id{ 'display_id' }            =~  m/^(\S*)\.\S*/;
                    my $seq_key                         =   $1;
                    
                    $alignment{$seq_key}->{'DNA_seq'}   =   $alignment_sequence;
                    
                }
            }
        }
        
        # Alignment block within transcript
        elsif ($sorf_begin <= $begin  && $end <= $sorf_end){
            for (my $i=1; $i<$num_sequences+1; $i++){
                
                my $alignment_sequence                  =   $aln->get_seq_by_pos($i)->seq;
                $alignment_sequence                     =~  tr/acgt/ACGT/;
                my $alignment_length                    =   $end - $begin;
                my $seq_id                              =   $aln->get_seq_by_pos($i);
                $$seq_id{ 'display_id' }                =~  m/^(\S*)\.\S*/;
                my $seq_key                             =   $1;
                
                $alignment{$seq_key}->{'DNA_seq'}       .=  $alignment_sequence;

            }
        }
        
        #Transcript spans over more as 1 alignment block, end of alignment
        elsif($begin <= $sorf_end && $sorf_end <= $end){
            
            my $offset                  =   $sorf_end - $begin + 1;
            my $count                   =   0;
            my $count_align             =   0;
            my $alignment_end_position  =   0;
            
            # Transcript begin position != alignment begin position because of multiple '-'
            foreach my $base (@char){
                
                # Get alignment position (bases and -'s)
                $count_align++;
                if ($base =~ /[ACTGactg]/){
                    
                    # Get alignment postion (only bases)
                    $count++;
                    if($count eq $offset){
                        $alignment_end_position = $count_align;
                    }
                }
            }
            for (my $i=1; $i<$num_sequences+1; $i++){
                
                my $alignment_sequence                  =   $aln->get_seq_by_pos($i)->seq;
                my $alignment_length                    =   $sorf_end - $begin;
                $alignment_sequence                     =   substr ($alignment_sequence,0,$alignment_end_position);
                $alignment_sequence                     =~  tr/acgt/ACGT/;
                my $seq_id                              =   $aln->get_seq_by_pos($i);
                $$seq_id{ 'display_id' }                =~  m/^(\S*)\.\S*/;
                my $seq_key                             =   $1;
                
                $alignment{$seq_key}->{'DNA_seq'}       .=  $alignment_sequence;
                
            }
        }
    }
    
    # Return
    return(\%alignment);
    
}

### Get MAF blocks overlapping sorf ###

sub SorfMafBlocks{
    
    # Catch
    my $sorfs   =   $_[0];
    my $sorf_id =   $_[1];
    my %MAF     =   %{$_[2]};
    
    #Init
    my %sORF_MAF    =   ();
    my $sorf_begin  =   $sorfs->{$sorf_id}{'sorf_begin'};
    my $sorf_end    =   $sorfs->{$sorf_id}{'sorf_end'};
    
    #Grep needed MAF blocks
    my @keys = grep { $MAF{$_}{'end'} > $sorf_begin -1 }  keys %MAF;
    @keys = grep {$_ < $sorf_end +1} @keys;
    @sORF_MAF{@keys} = @MAF{@keys};

    # Return
    return(dclone(\%sORF_MAF));
    
}

### Get MAF blocks from MAF file ###

sub MafToHash{
    
    # Input
    my $chr         =   $_[0];
    my $maf_root    =   $_[1];
    my $spec        =   $_[2];
    
    my $maf = {};
    
    # File handling
    #my $str = Bio::AlignIO->new(-file => $maf_root."/".$spec."/chr".$chr.".maf", format => "maf");
    my $str = Bio::AlignIO->new(-file => $maf_root."/".$spec."/tmp_chr".$chr.".maf", format => "maf");
    
    # Get seqIO object MAF blocks out of MAF file
    while ( my $aln = $str->next_aln()){
        my $seq = $aln->get_seq_by_pos(1);
        my $begin = $seq->{start};
        my $end = $seq->{end};
        
        # Read in blocks with start position as key, seqIO block as value
        $maf->{$begin}{'end'}       =   $end;
        $maf->{$begin}{'begin'}     =   $begin;
        $maf->{$begin}{'mafblock'}  =   $aln;
    }
    #print Dumper(%{$maf});
    
    # Return
    return($maf);
}

### match maf blocks to sorfs ###

sub match_maf_blocks_to_sorfs {

    #Catch
    my $maf_root    =   $_[0];
    my $chr         =   $_[1];
    my $spec        =   $_[2];
    my $sorfs       =   $_[3];
    
    #Init
    my @window = ();
    my ($sorf_array,$sorf_id);
    my $count;
    
    #if ($chr eq 'M'){$count = 0;}
    
    #Create array with sorf ids
    foreach $sorf_id (sort { $sorfs->{$a}{'sorf_begin'} <=> $sorfs->{$b}{'sorf_begin'} } keys %{$sorfs}){
        #if ($chr eq 'M'){$count++;}
        push (@$sorf_array,$sorf_id);
        #if ($chr eq 'M'){print $count."\t".$sorf_id."\n";}
    }
    
    #Loop over alignment blocks in maf file
    #my $str = Bio::AlignIO->new(-file => $maf_root."/".$spec."/chr".$chr.".maf", format => "maf");
    my $str = Bio::AlignIO->new(-file => $maf_root."/".$spec."/tmp_chr".$chr.".maf", format => "maf");

    # Get seqIO object MAF blocks out of MAF file
    while ( my $aln = $str->next_aln()){
        my $seq = $aln->get_seq_by_pos(1);
        my $begin = $seq->{start};
        my $end = $seq->{end};
        
        #print "MAF block between: ".$begin." and ".$end."\n";
        
        # Push all sorf_ids to @window where sorf_begin < window_pos
        foreach $sorf_id (@$sorf_array){
            if($sorfs->{$sorf_id}{'sorf_begin'} <= $end){
                push(@window,$sorf_id);
            }else{last;}
        }
        # Get rid of sorf_array elements already in @window
        @$sorf_array = grep { $sorfs->{$_}{'sorf_begin'} > $begin} @$sorf_array;
        
        # Get rid of sorf_ids in @$window where sorf_end < window_pos
        @window = grep { $sorfs->{$_}{'sorf_end'} >= $begin} @window;
        
        # Loop over window and add read position to window_genes
        foreach my $window_id (@window){
            $sorfs->{$window_id}->{'aln'}->{$begin}{'mafblock'}     =   $aln;
            $sorfs->{$window_id}->{'aln'}->{$begin}{'begin'}        =   $begin;
            $sorfs->{$window_id}->{'aln'}->{$begin}{'end'}          =   $end;
        }
        #if ($chr eq 'M'){print Dumper(@window);}
        #if ($chr eq 'M'){print Dumper($aln);}
    }
    
    #return
    return($sorfs);
}

### get sORFs per chromosome ###

sub get_sorfs_per_chromosome {
    
    # Catch
    my $dbh         =   $_[0];
    my $id          =   $_[1];
    my $chr         =   $_[2];
    
    # Init
    my $sorfs = {};
    
    # Get genes
    my $query = "SELECT sorf_id,chr,strand,sorf_begin,sorf_end from TIS_sORFs_".$id."_transcripts WHERE chr = '".$chr."'";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    $sorfs = $sth->fetchall_hashref('sorf_id');
    
    # Return
    return($sorfs);
}

### Get the analysis ids that need to be processed

sub get_analysis_ids {
    
    # Catch
    my $dbh    = $_[0];
    my $ids_in = $_[1]; #Either comma separated list of identifiers or "all"
    
    
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

### GET CHR SIZES FROM IGENOMES ###

sub get_chr_sizes {

    

    # Catch

    my %chr_sizes;

    my $filename = $igenomes_root."/".$spec."/Ensembl/".$assembly."/Annotation/Genes/ChromInfo.txt";

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

        $filename = $igenomes_root."/".$spec."/Ensembl/".$assembly."/Sequence/WholeGenomeFasta/GenomeSize.xml";

        open (Q,"<".$filename) || die "Cannot open chr sizes input\n";

        while (<Q>){

            if ($_ =~ /contigName=\"(.*)\".*totalBases=\"(\d+)\"/) {

                $chr_sizes{$1} = $2;

            }

        }

    }

    return(\%chr_sizes);

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
