#!/usr/bin/perl -w

# This script makes a gene distribution of the ribo-seq reads
# Output: Tabular table of the gene distibution together with three plots for visualization
###########

################
# COMMAND LINE:
# ./gene_distribution.pl --treated untreated (or 'treated') (--work_dir getcwd --in_sqlite SQLite/results.db --output_file run_name_genedistribution)
# **example: ./gene_distribution.pl --treated untreated (--work_dir /data/RIBO_runs/RIBO_Ingolia_GerbenM/ --in_sqlite SQLite/results.db --output_file mESC_GA_ens72_STAR_untreated_genedistribution)
# GALAXY
# gene_distribution.pl --tooldir "${__tool_directory__}" --in_sqlite "${sqlite_db}" --treated "${treatment}" --out_table "${out_table}" --out_pdf1 "${out_pdf1}" --out_pdf2 "${out_pdf2}" --out_pdf3 "${out_pdf3}"
################

use strict;
use warnings;
use DBI;
use Data::Dumper;
use Parallel::ForkManager;
use Getopt::Long;
use v5.10;
use Cwd;

################
## PARAMETERS & PREPARATION
################

# Parameters, get the command line arguments
my ($tooldir,$work_dir,$in_sqlite,$output_file,$treated,$out_table,$out_pdf1,$out_pdf2,$out_pdf3);
GetOptions(
#"cores=i"=>\$cores,                     # Number of cores to use 										mandatory argument
"tooldir:s"=>\$tooldir,                 # The directory the tool currently resides in (included for Galaxy) optional argument (default = $CWD env setting: for script version)
"work_dir:s" =>\$work_dir,              # Working directory												optional  argument (default = $CWD env setting)
"in_sqlite:s" =>\$in_sqlite,            # Results db (relative of CWD)                              	optional  argument (default = SQLite/results.db)
#"ens_db:s" =>\$ens_db,                  # ENSEMBL db (relative of CWD)                              	optional  argument (default = SQLite/name_build_from_arguments.db)
"output_file:s" =>\$output_file,		# Base name output files										optional  argument (default = 'run_name'_annotation)
#"igenomes_root:s" =>\$IGENOMES_ROOT,    # Igenomes root folder                                      	mandatory argument
"treated:s" => \$treated,				# Define data('no':no treatment/CHX or 'yes':treated(LTM/Harr))	mandatory argument
"out_table:s" =>\$out_table,			# Define output variables galaxy								optional argument (default = made with $output_file)
"out_pdf1:s" =>\$out_pdf1,				# Define output variables galaxy								optional argument (default = made with $output_file)
"out_pdf2:s" =>\$out_pdf2,				# Define output variables galaxy								optional argument (default = made with $output_file)
"out_pdf3:s" =>\$out_pdf3				# Define output variables galaxy								optional argument (default = made with $output_file)
);

###########################################################################
#Check all input variable and/or get default values and set extra variables
###########################################################################
# Input variables
#if($cores){
#    print "Number of cores						: $cores\n";
#}else{
#    die "\nDon't forget to pass number of cores to use for mapping using the --cores argument!\n\n";
#}
my $CWD = getcwd;
if($work_dir){
    print "Working directory					: $work_dir\n";
}else{
    $work_dir = $CWD;
    print "Working directory					: $CWD\n";
}
if($in_sqlite){
    print "Name of SQLite input db					: $in_sqlite\n";
}else{
    $in_sqlite = $work_dir."/SQLite/results.db";
    print "Name of SQLite input db					: $in_sqlite\n";
}
if ($tooldir){
    print "The tooldir is set to    : $tooldir\n";
} else {
    #Choose default value for tooldir
    $tooldir = $CWD;
    print "The tooldir is set to   : $tooldir\n";
}


# Get arguments vars
my $user = "";
my $pw = "";
my $db_ribo= $in_sqlite;
my ($species,$version,$run_name,$mapper,$IGENOMES_ROOT,$cores,$ens_db,$unique) = get_ARG_vars($db_ribo,$user,$pw);


# Create/define SQLite DBs/tables for the ribo-run
my $table_ribo;
if($treated eq 'untreated'){
    $table_ribo = ($unique eq 'N') ? 'count_fastq1_unique' : 'count_fastq1';
}elsif($treated eq 'treated'){
	$table_ribo = ($unique eq 'N') ? 'count_fastq2_unique' : 'count_fastq2';
}

# Igenomes
#$IGENOMES_ROOT = ($ENV{'IGENOMES_ROOT'}) ? $ENV{'IGENOMES_ROOT'} : $IGENOMES_ROOT;
print "The following igenomes folder is used			: $IGENOMES_ROOT\n";

# Cores
print "Number of cores to use for analysis			: $cores\n";

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

# Define ENSEMBL SQLite DB
#my $db_ensembl = ($ens_db) ? $ens_db : $work_dir."/SQLite/"."ENS_".$spec_short."_".$version.".db";
my $db_ensembl = $ens_db;
print "The following ENSEMBL db is used			: $db_ensembl\n";

# Create output directory and files
my $out_dir = $work_dir."/output/";
if(!-d "$out_dir"){
    system ("mkdir -p ".$out_dir);
}

if(!defined($output_file)){$output_file = $out_dir.$run_name."_ens".$version."_".$mapper."_".$treated;}
if(!defined($out_table)){$out_table = $output_file."_genedistribution.txt";}
if(!defined($out_pdf1)){$out_pdf1 = $output_file."_rankedgenes.pdf";}
if(!defined($out_pdf2)){$out_pdf2 = $output_file."_cumulative.pdf";}
if(!defined($out_pdf3)){$out_pdf3 = $output_file."_density.pdf";}

# Get chromosomes and correct coord_system_id
print "Get chromosomes and coord_system_id...\n";
my $chromosome_sizes; my $coord_system_id; my @ch;

$chromosome_sizes = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Annotation/Genes/ChromInfo.txt";
@ch = @{get_chr($species)};
#@ch = (17,18,19);
$coord_system_id = get_coord_system_id($db_ensembl,$assembly);

################
## GENE DISTRIBUTION OF RIBO-SEQ DATA
################

# Get StartTime
my $startRun = time();

# Open files
open OUT,"+>>".$out_table or die $!;
print OUT "GeneID\tread_count\n";

# Analysis
#gene_distribution($db_ribo,$db_ensembl,$table_ribo,\@ch,$coord_system_id,$user,$pw);
gene_distribution($db_ribo,$db_ensembl,$table_ribo,\@ch,$coord_system_id,$user,$pw,$cores);

# Make Plots
make_plots($out_table,$out_pdf1,$out_pdf2,$out_pdf3,$tooldir);

timer($startRun);	# Get Run time

print "DONE!\n";

################
# 	THE SUBS   #
################

sub get_ARG_vars{
	# Catch
	my $db_ribo = $_[0];
	my $user = $_[1];
	my $pw = $_[2];

    my ($query,$sth);

    # Connect to db
    my $dbh_results = DBI->connect('DBI:SQLite:'.$db_ribo,$user,$pw,
    							{ RaiseError => 1},) || die "Database connection not made: $DBI::errstr";

    # Get ENS variables
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

    $query = "select value from arguments where variable = \'unique\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $unique = $sth->fetch()->[0];
    $sth->finish();

	$query = "select value from arguments where variable = \'run_name\'";
    $sth = $dbh_results->prepare($query);
	$sth->execute();
	my $run_name = $sth->fetch()->[0];
	$sth->finish();

	$query = "select value from arguments where variable = \'mapper\'";
    $sth = $dbh_results->prepare($query);
	$sth->execute();
	my $mapper = $sth->fetch()->[0];
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

	$query = "select value from arguments where variable = \'ens_db\'";
    $sth = $dbh_results->prepare($query);
	$sth->execute();
	my $ens_db = $sth->fetch()->[0];
	$sth->finish();

	$dbh_results -> disconnect();

	# Return ENS variables
	return($species,$version,$run_name,$mapper,$IGENOMES_ROOT,$cores,$ens_db,$unique);
} # Close sub

sub get_chr{
    # Catch
    my $species = $_[0];

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
    my @chr = keys %chr_sizes;
    foreach my $ch (@chr){
    	if ($species eq "fruitfly"){
    		if($ch eq "M"){
    			$ch = "dmel_mitochondrion_genome";
        	}else{
            	$ch = $ch;
        	}
     	}else {
     		$ch = $ch;
     	}
    }

    return(\@chr);
} # Close sub

sub get_coord_system_id{
	# Catch
	my $db_ensembl = $_[0];
	my $assembly = $_[1];

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

sub gene_distributionOLD{
	# Catch
	my $db_ribo = $_[0];
	my $db_ensembl = $_[1];
	my $table_ribo = $_[2];
	my @ch = @{$_[3]};
	my $coord_system_id = $_[4];
	my $user = $_[5];
	my $pw = $_[6];

	# Loop through chromosomes
	foreach my $chr(@ch){
		##############
		## ENSEMBL -> GENES
		##############

		# Connect to Ensembl db
		my $dbh = DBI->connect('DBI:SQLite:'.$db_ensembl,$user,$pw,
							{RaiseError => 1},) || die "Database connection not made: $DBI::errstr";

		# Get seq_region_id for specific chromosome from table 'seq_region'
		my $query = "SELECT seq_region_id FROM seq_region WHERE coord_system_id = '$coord_system_id' AND name = '$chr'";
		my $execute = $dbh->prepare($query);
		$execute->execute();

		my $seq_region;
		while(my @result = $execute->fetchrow_array()){
			#$result[0]: seq_region_id

			$seq_region = $result[0];
		}
		$execute->finish();

		# Get all genes with start and stop position for specific chromosome from table gene
		my $query1 = "SELECT stable_id,seq_region_start,seq_region_end,seq_region_strand FROM gene WHERE seq_region_id = '$seq_region'";
		my $execute1 = $dbh->prepare($query1);
		$execute1->execute();

		my %genes;
		while(my @result1 = $execute1->fetchrow_array()){
			#$result1[0]: gene stable_id
			#$result1[1]: gene seq_region_start
			#$result1[2]: gene seq_region_end
			#$result1[3]: gene seq_region_strand

			for(my $i=$result1[1];$i<=$result1[2];$i++){
				$genes{$chr.':'.$i}{$result1[0]} = $result1[3];
			}
		}
		$execute1->finish();

		# Disconnect Ensembl db
		$dbh->disconnect();

		##############
		## RIBO-SEQ -> READs (~A-site position): determine gene distribution
		##############
		# Connect to ribo-seq db
		my $dbh2 = DBI->connect('DBI:SQLite:'.$db_ribo,$user,$pw,
							{RaiseError => 1},) || die "Database connection not made: $DBI::errstr";

		# Get all reads for a specific chromosome, look if they fall in a gene and save the total read count per gene
		my $query2 = "SELECT chr,strand,start,count FROM $table_ribo WHERE chr = '$chr'";
		my $execute2 = $dbh2->prepare($query2);
		$execute2->execute();

		my %gene_count;
		my $intergenic_count;
		while(my @result2 = $execute2->fetchrow_array()){
			#$result2[0]: chr
			#$result2[1]: strand
			#$result2[2]: start (~A-site position)
			#$result2[3]: read_count

			my $pos = $result2[0].':'.$result2[2];
			if(defined($genes{$pos})){
				#print $pos."\n";
				foreach my $id(keys %{$genes{$pos}}){
					if($genes{$pos}{$id} == $result2[1]){
						#print $result2[3]."\n";
						$gene_count{$id} = $gene_count{$id} + $result2[3];
					}
				}
			}else{
				$intergenic_count = $intergenic_count + $result2[3];
			}
		}
		$execute2->finish();

		# Disconnect ribo-seq db
		$dbh2->disconnect();

		##############
		## RESULTS: make table
		##############

		foreach my $id(keys %gene_count){
			print OUT $id."\t".$gene_count{$id}."\n";
		}
		print $intergenic_count."\n";
	}
}

sub gene_distribution{
	print "Get Gene Distribution...\n";

	# Catch
	my $db_ribo = $_[0];
	my $db_ensembl = $_[1];
	my $table_ribo = $_[2];
	my @ch = @{$_[3]};
	my $coord_system_id = $_[4];
	my $user = $_[5];
	my $pw = $_[6];
	my $cores = $_[7];

	# Init multi core
	my $processes = $cores; # Nr of processes
	my $pm = new Parallel::ForkManager($processes); # Open fork

	# Loop through chromosomes
	foreach my $chr(@ch){
		print "\tStarting analysis for chr ".$chr."...\n";
		# Start parallel process
		$pm->start and next;

		##############
		## ENSEMBL -> GENES
		##############
		print "\t\tRead in Ensembl data of chr ".$chr."\n";

		# Connect to Ensembl db
		my $dbh = DBI->connect('DBI:SQLite:'.$db_ensembl,$user,$pw,
							{RaiseError => 1},) || die "Database connection not made: $DBI::errstr";

		# Get seq_region_id for specific chromosome from table 'seq_region'
		my $query = "SELECT seq_region_id FROM seq_region WHERE coord_system_id = '$coord_system_id' AND name = '$chr'";
		my $execute = $dbh->prepare($query);
		$execute->execute();

		my $seq_region;
		while(my @result = $execute->fetchrow_array()){
			#$result[0]: seq_region_id

			$seq_region = $result[0];
		}
		$execute->finish();

		# Get all genes with start and stop position for specific chromosome from table gene
		my $query1 = "SELECT stable_id,seq_region_start,seq_region_end,seq_region_strand FROM gene WHERE seq_region_id = '$seq_region'";
		my $execute1 = $dbh->prepare($query1);
		$execute1->execute();

		my %genes;
		while(my @result1 = $execute1->fetchrow_array()){
			#$result1[0]: gene stable_id
			#$result1[1]: gene seq_region_start
			#$result1[2]: gene seq_region_end
			#$result1[3]: gene seq_region_strand

			$genes{$chr.":".$result1[3]}{$result1[0]} = [$result1[1],$result1[2]];
		}
		$execute1->finish();

		# Disconnect Ensembl db
		$dbh->disconnect();

		##############
		## RIBO-SEQ -> READs (~A-site position): determine gene distribution
		##############
		print "\t\tAnnotating ribo-seq reads of chr ".$chr."\n";

		# Connect to ribo-seq db
		my $dbh2 = DBI->connect('DBI:SQLite:'.$db_ribo,$user,$pw,
							{RaiseError => 1},) || die "Database connection not made: $DBI::errstr";

		# Get all reads for a specific chromosome, look if they fall in a gene and save the total read count per gene
		my $query2 = "SELECT chr,strand,start,count FROM $table_ribo WHERE chr = '$chr'";
		my $execute2 = $dbh2->prepare($query2);
		$execute2->execute();

		my %gene_count;
		my $intergenic_count;
		while(my @result2 = $execute2->fetchrow_array()){
			#$result2[0]: chr
			#$result2[1]: strand
			#$result2[2]: start (~A-site position)
			#$result2[3]: read_count

			my $pos = $result2[2];
			my $count = $result2[3];
			my $def = 0;

			#print "chr".$result2[0]."\t".$result2[1]."\t".$result2[2]."\n";

			foreach my $id(keys %{$genes{$result2[0].":".$result2[1]}}){
				#print $id."\n";
				my $gene_start = $genes{$result2[0].":".$result2[1]}{$id}[0];
				my $gene_stop = $genes{$result2[0].":".$result2[1]}{$id}[1];

				if($pos >= $gene_start && $pos <= $gene_stop){
					#print $id."\n";
					$gene_count{$id} += $count;
					$def++; # defined in genes
				}
			}

			if($def == 0){
				$intergenic_count +=  $count;
			}
		}
		$execute2->finish();

		# Disconnect ribo-seq db
		$dbh2->disconnect();

		##############
		## RESULTS: make table
		##############

		foreach my $id(keys %gene_count){
			print OUT $id."\t".$gene_count{$id}."\n";
		}
		#print $intergenic_count."\n";

		# Close process
		$pm->finish;
	}

	# Finish forking
	$pm->wait_all_children;
}

sub make_plots{
	print "Make Gene Distribution Plots...\n";

	# Catch
	my $out_table = $_[0];
	my $out_pdf1 = $_[1];
	my $out_pdf2 = $_[2];
	my $out_pdf3 = $_[3];
	my $tooldir = $_[4];

	# Execute Rscript
    #system("Rscript ".$galaxydir."/tools/proteoformer/quality_plots.R ".$out_table." ".$out_pdf1." ".$out_pdf2." ".$out_pdf3);
    system("Rscript ".$tooldir."/quality_plots.R ".$out_table." ".$out_pdf1." ".$out_pdf2." ".$out_pdf3);
    #system("Rscript quality_plots.R ".$out_table." ".$out_pdf1." ".$out_pdf2." ".$out_pdf3);
}

sub timer {
    my $startRun = shift;

    my $endRun 	= time();
    my $runTime = $endRun - $startRun;

    printf("\nTotal running time: %02d:%02d:%02d\n\n", int($runTime / 3600), int(($runTime  % 3600) / 60), int($runTime % 60));
}
