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


# This script determines the translation-level of transcripts (exon_level) in ribo-seq experiments
################

################
# COMMAND LINE:
# ./ribo_translation.pl (--work_dir getcwd --in_sqlite SQLite/results.db --ens_db SQLite/ens.db --tmp $TMP --out_sqlite SQLite/results.db)
# ** example: ./ribo_translation.pl (--work_dir /path/to/workdir/ --in_sqlite SQLite/results.db --ens_db SQLite/ens.db --tmp path/to/tmp/ --out_sqlite SQLite/results.db)
# GALAXY
# ribo_translation.pl --in_sqlite "${sqlite_db}" --ens_db "${ensembl_db.fields.path}"  --out_sqlite "${out_sqlite_db}"
################

use strict;
use warnings;
use DBI;
use Data::Dumper;
use Parallel::ForkManager;
use Getopt::Long;
use Cwd;

################
## PARAMETERS & PREPARATION
################

# Parameters, get the command line arguments
my ($work_dir,$in_sqlite,$ens_db,$out_sqlite,$tmpfolder);
my $help;
GetOptions(
"work_dir:s" =>\$work_dir,              # Working directory											optional  argument (default = $CWD env setting)
"in_sqlite:s" =>\$in_sqlite,            # Results db (relative of CWD)                              optional  argument (default = SQLite/results.db)
"ens_db:s" =>\$ens_db,                  # ENSEMBL db (relative of CWD)                              optional  argument (default = SQLite/name_build_from_arguments.db)
"tmp:s" =>\$tmpfolder,					# Folder where temporary files are stored,                  optional  argument (default = $TMP env setting)
"out_sqlite:s" =>\$out_sqlite,          # Directory of out results db (relative of CWD)             optional  argument (default = SQLite/results.db)
"help" => \$help                        # Help text option
);

if ($help) {
    print_help_text();
    exit;
}

###########################################################################
#Check all input variable and/or get default values and set extra variables
###########################################################################
# Input variables
my $CWD = getcwd;
if($work_dir){
    print "Working directory					: $work_dir\n";
}else{
    $work_dir = $CWD."/";
    print "Working directory					: $CWD\n";
}
if($in_sqlite){
    print "Name of SQLite input db					: $in_sqlite\n";
}else{
    $in_sqlite = $work_dir."/SQLite/results.db";
    print "Name of SQLite input db					: $in_sqlite\n";
}
if($out_sqlite){
	print "Name of SQLite output db				: $out_sqlite\n";
}else{
	$out_sqlite = $in_sqlite;
	print "Name of SQLite output db				: $out_sqlite\n";
}
my $TMP             = ($ENV{'TMP'}) ? $ENV{'TMP'} : ($tmpfolder) ? $tmpfolder : "$CWD/tmp/" ; # First select the TMP environment variable, or secondly select the $tmpfolder variable, or finally select current_working_dir/tmp
if($tmpfolder){
	print "The following tmpfolder is used				: $tmpfolder\n";
}else{
	$tmpfolder = $TMP;
	print "The following tmpfolder is used				: $TMP\n";
}

#Check if tmpfolder exists, if not create it...
if(!-d "$TMP"){
    system ("mkdir ". $TMP);
}

# Create/define SQLite DBs/tables for the ribo-run
my $db_ribo= $in_sqlite;
my $table_ribo = "count_fastq1";
my $table_ribo_trans = "tr_translation";

my $user = "";
my $pw = "";
my $dbh_create = DBI->connect('DBI:SQLite:dbname='.$db_ribo,$user,$pw,
										{ RaiseError => 1},) || die "Database connection not made: ".$DBI::errstr;

my $transcript_drop = "DROP TABLE IF EXISTS $table_ribo_trans";
$dbh_create->do($transcript_drop);

my $transcript_create = "CREATE TABLE IF NOT EXISTS $table_ribo_trans (
		transcript_id VARCHAR(100) NOT NULL,
		stable_id VARCHAR(100) NOT NULL,
		chr char(50) NOT NULL default '',
		seq_region_id VARCHAR(10) NOT NULL,
		seq_region_strand VARCHAR(2) NOT NULL,
		seq_region_start FLOAT NOT NULL,
		seq_region_end FLOAT NOT NULL,
		read_counts FLOAT NOT NULL,
		normalized_counts FLOAT NOT NULL,
		biotype VARCHAR(100) NOT NULL,
		exon_coverage VARCHAR(5) NOT NULL,
		canonical VARCHAR(5) NOT NULL,
		ccds VARCHAR(20) NOT NULL,
		gene_stable_id VARCHAR(100) NOT NULL,
        FPKM FLOAT NOT NULL,
        coverage FLOAT NOT NULL,
		PRIMARY KEY (stable_id,gene_stable_id)
		)";
$dbh_create->do($transcript_create);

my $name_idx1 = $table_ribo_trans."_fastq1_exon_coverage";
my $transcript_idx1 = "CREATE INDEX IF NOT EXISTS $name_idx1 ON $table_ribo_trans (exon_coverage)";
$dbh_create->do($transcript_idx1);

my $name_idx2 = $table_ribo_trans."_fastq1_seq_region_id";
my $transcript_idx2 = "CREATE INDEX IF NOT EXISTS $name_idx2 ON $table_ribo_trans (seq_region_id)";
$dbh_create->do($transcript_idx2);

$dbh_create->disconnect();

# Get arguments vars
my ($species,$version,$IGENOMES_ROOT,$cores,$run_name,$mapper,$uniq) = get_ARG_vars($db_ribo,$user,$pw);

# Igenomes
print "The following igenomes folder is used			: $IGENOMES_ROOT\n";

# Cores
print "Number of cores to use for analysis			: $cores\n";

#Conversion for species terminology
my $spec = ($species eq "mouse") ? "Mus_musculus" 
: ($species eq "human") ? "Homo_sapiens" 
: (uc($species) eq "ARCTIC_SQUIRREL") ? "Urocitellus_parryii" 
: ($species eq "arabidopsis") ? "Arabidopsis_thaliana" 
: ($species eq "fruitfly") ? "Drosophila_melanogaster" 
: (uc($species) eq "CNECNA3") ? "Cryptococcus_neoformans_var_grubii_h99_gca_000149245" 
: ($species eq "SL1344") ? "SL1344" : "";
my $spec_short = ($species eq "mouse") ? "mmu" 
: ($species eq "human") ? "hsa" 
: (uc($species) eq "ARCTIC_SQUIRREL") ? "upa" 
: (uc($species) eq "CNECNA3") ? "cnecna3" 
: ($species eq "arabidopsis") ? "ath" : ($species eq "fruitfly") ? "dme" 
: ($species eq "SL1344") ? "SL1344" : "";
#Old mouse assembly = NCBIM37, new one is GRCm38. Old human assembly = GRCh37, the new one is GRCh38
my $assembly = (uc($species) eq "MOUSE" && $version >= 70 ) ? "GRCm38"
: (uc($species) eq "MOUSE" && $version < 70 ) ? "NCBIM37"
: (uc($species) eq "HUMAN" && $version >= 76) ? "GRCh38"
: (uc($species) eq "HUMAN" && $version < 76) ? "GRCh37"
: (uc($species) eq "ARCTIC_SQUIRREL" && $version > 95) ? "ASM342692v1"
: (uc($species) eq "ARABIDOPSIS") ? "TAIR10"
: (uc($species) eq "SL1344") ? "ASM21085v2"
: (uc($species) eq "CNECNA3") ? "CNA3"
: (uc($species) eq "FRUITFLY" && $version < 79) ? "BDGP5"
: (uc($species) eq "FRUITFLY" && $version >= 79) ? "BDGP6" : "";

# Define ENSEMBL SQLite DB
my $db_ensembl = ($ens_db) ? $ens_db : $work_dir."/SQLite/"."ENS_".$spec_short."_".$version.".db";
print "The following ENSEMBL db is used			: $db_ensembl\n\n";

# Get chromosomes and correct coord_system_id
print "Get chromosomes and coord_system_id...\n";
my $chromosome_sizes; my $coord_system_id; my @ch;

$chromosome_sizes = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Annotation/Genes/ChromInfo.txt";
@ch = @{get_chr($chromosome_sizes,$species)};
print Dumper (@ch);
#@ch = ("Y","MT");
$coord_system_id = get_coord_system_id($db_ensembl,$assembly,'chromosome');
my $coord_system_id_plasmid = '';
# Get coord system id for the SL1344 plasmids, if necessary
if(uc($species) eq "SL1344"){
    $coord_system_id_plasmid = get_coord_system_id($ens_db,$assembly,'plasmid');
}

# Store used path of ENSEMBL db in arguments table
store_ENS_var($db_ribo,$user,$pw,$db_ensembl);

# Construct sample name
my $sample_name = $run_name."_".$mapper."_".$uniq."_".$version.".fastq1";

################
## BASED ON RIBO-SEQ DATA: GET TRANSLATION-LEVEL FOR EACH GENE AND FOR EACH TRANSCRIPT
################

# Translation-info per chr in a file
translation_per_chr($sample_name);

# Combine output for all chr in SQLite dumpfile
make_sqlite_dumpfile();

# import SQLite dump
import_sqlite_dumpfile($tmpfolder."fastq1_transcript_all.csv");

# Add transcript calling method
update_arguments_table($in_sqlite);

# Move to galaxy history
system ("mv ".$in_sqlite." ".$out_sqlite);

################
# 	THE SUBS   #
################

sub update_arguments_table{
    #Catch
    my $db = $_[0];
    my $user = "";
    my $pw = "";
    
    # Connect to sqlite database
    my $dbh  = DBI->connect('DBI:SQLite:'.$db,$user,$pw,
    { RaiseError => 1},) || die "Database connection not made: $DBI::errstr";
    
    #Update arguments table
    my $query = "DELETE FROM arguments WHERE value='tr_calling';";
    my $execute = $dbh->prepare($query);
    $execute->execute();
    $execute->finish();
    $query = "INSERT INTO arguments (variable, value) VALUES ('tr_calling', 'rule_based');";
    $execute = $dbh->prepare($query);
    $execute->execute();
    $execute->finish();
    
    return;
}

sub get_chr{
    # Catch
    my $chromosome_sizes = $_[0];
    my $species = $_[1];
    
    # Work
    my @chr;
    my $count = 0;
    open (Q,"<".$chromosome_sizes) || die "Cannot open chr sizes input\n";
    while (<Q>){
        my @a = split(/\s+/,$_);
        $chr[$count] = $a[0];
        $count++;
    }

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
    my $name = $_[2];

	# Connect to ensembl sqlite database
	my $dbh  = DBI->connect('DBI:SQLite:'.$db_ensembl,$user,$pw,
						{ RaiseError => 1},) || die "Database connection not made: $DBI::errstr";

	# Get correct coord_system_id
	my $query = "SELECT coord_system_id FROM coord_system WHERE name = '$name' AND version = '$assembly'";
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

sub translation_per_chr{
    #Catch
    my $sample_name = $_[0];
    
	print "Starting translation analysis of transcripts per chromosome using ".$cores." cores...\n";
	# Init multi core
	my $processes = $cores; # Nr of processes
	my $pm = new Parallel::ForkManager($processes); # Open fork

	# Loop through chromosomes
	foreach my $chr(@ch){
		# Start parallel process
		$pm->start and next;

		###########
		## ENSEMBL_DB
		###########
		# Connect to ensembl sqlite database
		my $dbh  = DBI->connect('DBI:SQLite:'.$db_ensembl,$user,$pw,
								{ RaiseError => 1},) || die "Database connection not made: $DBI::errstr";

		# Get seq_region_id for specific chromosome from table 'seq_region'
        my $query1 = "";
        if ($chr =~ m/.+_SL1344$/){
            $query1 = "SELECT seq_region_id FROM seq_region WHERE coord_system_id = '$coord_system_id_plasmid' AND name = '$chr'";
        } else {
            $query1 = "SELECT seq_region_id FROM seq_region WHERE coord_system_id = '$coord_system_id' AND name = '$chr'";
        }
		my $execute1 = $dbh->prepare($query1);
		$execute1->execute();

		my $seq_region;
		while(my @result1 = $execute1->fetchrow_array()){

			$seq_region = $result1[0];
		}
		$execute1->finish();
        
		####
		# Get all genes and corresponding transcripts for specific chromosome from joined tables 'gene', 'transcript' and 'translation'
		####
		## First for PROTEIN-CODING transcripts -> UTR-info
		my $query2 = "SELECT b.gene_id, b.stable_id, a.transcript_id, a.seq_region_start, a.seq_region_end, a.seq_region_strand, a.stable_id, c.start_exon_id, c.end_exon_id, c.seq_start, c.seq_end, a.biotype, b.canonical_transcript_id
		FROM transcript a
		INNER JOIN gene b ON a.gene_id = b.gene_id
		INNER JOIN translation c ON a.transcript_id = c.transcript_id
		WHERE b.seq_region_id = '$seq_region'";

		my $execute2 = $dbh->prepare($query2);
		$execute2->execute();

		my %gene_count;
		my %transcript_count;
        my %transcript_covered_positions;
		my %transcript_id;
		my %transcript_coding;
		my %gene_id;
		my %transcript_attributes;
		while(my @result2 = $execute2->fetchrow_array()){

			# Gene info
			$gene_id{$result2[1]} = [$result2[0],$result2[5]];
			$gene_count{$result2[1]} = 0; # Initialize count hash

			# Transcript info
			$transcript_id{$result2[2]} = [$result2[6],$result2[11],$result2[3],$result2[4]];
			$transcript_coding{$result2[2]} = [$result2[1],$result2[5],$result2[7],$result2[8],$result2[9],$result2[10]];
			$transcript_count{$result2[1]}{$result2[2]} = 0; # Initialize count hash
            $transcript_covered_positions{$result2[1]}{$result2[2]} = 0; # Initialize cov positions hash
            # Attr1 = canonical; Attr2 = ccds
            my $canonical = ($result2[2] eq $result2[12]) ? "Yes" : "No";
            my $ccds = "No";
			$transcript_attributes{$result2[2]} = [$canonical,$ccds];
		}
		$execute2->finish();

		#### transcripts with CCDS info
		my $queryCCDS = "SELECT a.transcript_id, e.dbprimary_acc
		FROM transcript a
		INNER JOIN gene b ON a.gene_id = b.gene_id
		INNER JOIN translation c ON a.transcript_id = c.transcript_id
		INNER JOIN object_xref d ON a.transcript_id = d.ensembl_id
		INNER JOIN xref e ON d.xref_id = e.xref_id WHERE b.seq_region_id = '$seq_region' AND e.external_db_id = '3800'";

		my $executeCCDS = $dbh->prepare($queryCCDS);
		$executeCCDS->execute();

		while(my @resultCCDS = $executeCCDS->fetchrow_array()){

			my $ccds = $resultCCDS[1];
			$transcript_attributes{$resultCCDS[0]}[1] = $ccds;
		}
		$executeCCDS->finish();

		## Next for the NON-PROTEIN-CODING transcripts -> no UTR/translation-info
		my $query3 = "SELECT b.gene_id, b.stable_id, a.transcript_id, a.seq_region_start, a.seq_region_end, a.seq_region_strand, a.stable_id, a.biotype, b.canonical_transcript_id
		FROM transcript a
		INNER JOIN gene b ON a.gene_id = b.gene_id
		LEFT OUTER JOIN translation c ON a.transcript_id = c.transcript_id
		WHERE c.transcript_id IS NULL AND b.seq_region_id = '$seq_region'";

		my $execute3 = $dbh->prepare($query3);
		$execute3->execute();

		my %transcript_noncoding;
		while(my @result3 = $execute3->fetchrow_array()){
			# Gene info
			$gene_id{$result3[1]} = [$result3[0],$result3[5]];
			$gene_count{$result3[1]} = 0; # Initialize count hash

			# Transcript info
			$transcript_id{$result3[2]} = [$result3[6],$result3[7],$result3[3],$result3[4]];
			$transcript_noncoding{$result3[2]} = [$result3[1],$result3[5]];
			$transcript_count{$result3[1]}{$result3[2]} = 0; # Initialize count hash
            $transcript_covered_positions{$result3[1]}{$result3[2]} = 0; # Initialize cov positions hash
			# Attr1 = canonical; Attr2 = ccds
            my $canonical = ($result3[2] eq $result3[8]) ? "Yes" : "No";
            my $ccds  = "No";
			$transcript_attributes{$result3[2]} = [$canonical,$ccds];
		}
		$execute3->finish();

		#### transcripts with CCDS info
		my $queryCCDS2 = "SELECT a.transcript_id, a.seq_region_start, a.seq_region_end, a.seq_region_strand, a.stable_id, a.biotype, b.canonical_transcript_id, e.dbprimary_acc
		FROM transcript a
		INNER JOIN gene b ON a.gene_id = b.gene_id
		INNER JOIN object_xref d ON a.transcript_id = d.ensembl_id
		INNER JOIN xref e ON d.xref_id = e.xref_id
		LEFT OUTER JOIN translation c ON a.transcript_id = c.transcript_id
		WHERE c.transcript_id IS NULL AND b.seq_region_id = '$seq_region' AND e.external_db_id = '3800'";

		my $executeCCDS2 = $dbh->prepare($queryCCDS2);
		$executeCCDS2->execute();

		while(my @resultCCDS2 = $executeCCDS2->fetchrow_array()){

			my $ccds = $resultCCDS2[1];
			$transcript_attributes{$resultCCDS2[0]}[1] = $ccds;
		}
		$executeCCDS2->finish();

		####
		# Get all exons(positions) per transcript
		####
		## First for NON-PROTEIN-CODING transcripts from joined tables 'exon' and 'exon_transcript'
		my %exon_gene_transcript;
		my %transcript_length;
		my %exon_length;
		my %exon_id;
		my %exon_count;
		my %transcript_exon;
		foreach my $transcript(keys %transcript_noncoding){
			my $query4 = "SELECT b.exon_id,a.seq_region_start,a.seq_region_end,a.seq_region_strand,a.stable_id,b.rank FROM exon a JOIN exon_transcript b ON a.exon_id = b.exon_id WHERE b.transcript_id = '$transcript'";
			my $execute4 = $dbh->prepare($query4);
			$execute4->execute();

			while(my @result4 = $execute4->fetchrow_array()){

				# Exonic position info
				for(my $i=$result4[1];$i<=$result4[2];$i++){
					$exon_gene_transcript{$i}{$transcript_noncoding{$transcript}[0]}{$transcript}{$result4[0]} = $result4[3];
				}

				# Exon info
				$exon_id{$result4[0]}{$transcript} = $result4[4];
				$exon_length{$transcript}{$result4[0]} = $result4[2] - $result4[1] + 1;
				$exon_count{$transcript}{$result4[0]} = 0; # Initialize count hash

				# Transcript info
				$transcript_length{$transcript} += $exon_length{$transcript}{$result4[0]};

				# Transcript-exon info
				$transcript_exon{$transcript}{$result4[5]} = [$result4[0],$result4[1],$result4[2],$result4[3]];
			}
			$execute4->finish();
		}

		## Next for the PROTEIN-CODING transcripts
		foreach my $transcript(keys %transcript_coding){

			my $id_first = $transcript_coding{$transcript}[2];
			my $id_last = $transcript_coding{$transcript}[3];
			my $seq_start = $transcript_coding{$transcript}[4];
			my $seq_end = $transcript_coding{$transcript}[5];

			# Determine start codon and thus the protein-coding region of the first exon
			my $start_codon;
			my $query5 = "SELECT a.exon_id,a.seq_region_start,a.seq_region_end,a.seq_region_strand,a.stable_id,b.rank FROM exon a JOIN exon_transcript b ON a.exon_id = b.exon_id WHERE a.exon_id = '$id_first'";
			my $execute5 = $dbh->prepare($query5);
			$execute5->execute();

			while(my @result5 = $execute5->fetchrow_array()){

				# Exon info
				$exon_id{$result5[0]}{$transcript} = $result5[4];

				# Transcript-exon info
				if($result5[3] == 1){ # Forward strand
					$start_codon = $result5[1] + $seq_start - 1;
					$transcript_exon{$transcript}{$result5[5]} = [$result5[0],$start_codon,$result5[2],$result5[3]];
				}elsif($result5[3] == -1){ # Reverse strand
					$start_codon = $result5[2] - $seq_start + 1;
					$transcript_exon{$transcript}{$result5[5]} = [$result5[0],$result5[1],$start_codon,$result5[3]];
				}
			}
			$execute5->finish();

			# Determine stop codon and thus the protein-coding region of the last exon
			my $stop_codon;
			my $query6 = "SELECT a.exon_id,a.seq_region_start,a.seq_region_end,a.seq_region_strand,a.stable_id,b.rank FROM exon a JOIN exon_transcript b ON a.exon_id = b.exon_id WHERE a.exon_id = '$id_last'";
			my $execute6 = $dbh->prepare($query6);
			$execute6->execute();

			while(my @result6 = $execute6->fetchrow_array()){

				# Exon info
				$exon_id{$result6[0]}{$transcript} = $result6[4];

				# Transcript-exon info
				if($result6[3] == 1){ # Forward strand
					$stop_codon  = $result6[1] + $seq_end - 1;
					$transcript_exon{$transcript}{$result6[5]} = [$result6[0],$result6[1],$stop_codon,$result6[3]];
				}elsif($result6[3] == -1){ # Reverse strand
					$stop_codon = $result6[2] - $seq_end + 1;
					$transcript_exon{$transcript}{$result6[5]} = [$result6[0],$stop_codon,$result6[2],$result6[3]];
				}
			}
			$execute6->finish();

			# Determine other protein-coding exons
			my $query7 = "SELECT b.exon_id,a.seq_region_start,a.seq_region_end,a.seq_region_strand,a.stable_id,b.rank FROM exon a JOIN exon_transcript b ON a.exon_id = b.exon_id WHERE b.transcript_id = '$transcript' AND a.exon_id != '$id_first' AND a.exon_id != '$id_last'";
			my $execute7 = $dbh->prepare($query7);
			$execute7->execute();

			while(my @result7 = $execute7->fetchrow_array()){

				# Determine if exon is a protein_coding exon in this specific transcript
				if(($result7[1]>=$start_codon && $result7[2]<=$stop_codon) || ($result7[2]<=$start_codon && $result7[1]>=$stop_codon)){
					# Exon info
					$exon_id{$result7[0]}{$transcript} = $result7[4];

					# Transcript-exon info
					$transcript_exon{$transcript}{$result7[5]} = [$result7[0],$result7[1],$result7[2],$result7[3]];
				}
			}
			$execute7->finish();

			# ASSEMBLE TRANSCRIPT
			#####################
			my $rank1 = (sort {$a <=> $b} keys %{$transcript_exon{$transcript}})[0];
			my $rank2 = (sort {$b <=> $a} keys %{$transcript_exon{$transcript}})[0];

			# Clip piece of START & STOP codon
			my $clip_start = 15;
			my $clip_stop = 15;
			my ($new_start,$new_rank1); my ($new_stop,$new_rank2);

			# Forward strand
			if($transcript_exon{$transcript}{$rank1}[3] == 1){
				# Clip from START CODON
				my $bases = 0;
				foreach my $rank(sort {$a <=> $b} keys %{$transcript_exon{$transcript}}){
					for(my $i=$transcript_exon{$transcript}{$rank}[1];$i<=$transcript_exon{$transcript}{$rank}[2];$i++){
						$bases++;
						$new_start = $i;
						$new_rank1 = $rank;
						last if $bases == $clip_start + 1;
					}
					last if $bases == $clip_start + 1;
				}
				# Clip from STOP CODON
				$bases = 0;
				foreach my $rank(sort {$b <=> $a} keys %{$transcript_exon{$transcript}}){
					for(my $i=$transcript_exon{$transcript}{$rank}[2];$i>=$transcript_exon{$transcript}{$rank}[1];$i--){
						$bases++;
						$new_stop = $i;
						$new_rank2 = $rank;
						last if $bases == $clip_stop + 1;
					}
					last if $bases == $clip_stop + 1;
				}
				foreach my $rank(sort {$a <=> $b} keys %{$transcript_exon{$transcript}}){
					if($rank == $new_rank1){ # First exon
						# Exonic position info
						for(my $i=$new_start;$i<=$transcript_exon{$transcript}{$rank}[2];$i++){
							$exon_gene_transcript{$i}{$transcript_coding{$transcript}[0]}{$transcript}{$transcript_exon{$transcript}{$rank}[0]} = $transcript_exon{$transcript}{$rank}[3];
						}

						# Exon info
						$exon_length{$transcript}{$transcript_exon{$transcript}{$rank}[0]} = $transcript_exon{$transcript}{$rank}[2] - $new_start + 1;
						$exon_count{$transcript}{$transcript_exon{$transcript}{$rank}[0]} = 0; # Initialize count hash

						# Transcript info
						$transcript_length{$transcript} += $exon_length{$transcript}{$transcript_exon{$transcript}{$rank}[0]};
					}elsif($rank == $new_rank2){ # Last exon
						# Exonic position info
						for(my $i=$transcript_exon{$transcript}{$rank}[1];$i<=$new_stop;$i++){
							$exon_gene_transcript{$i}{$transcript_coding{$transcript}[0]}{$transcript}{$transcript_exon{$transcript}{$rank}[0]} = $transcript_exon{$transcript}{$rank}[3];
						}

						# Exon info
						$exon_length{$transcript}{$transcript_exon{$transcript}{$rank}[0]} = $new_stop - $transcript_exon{$transcript}{$rank}[1] + 1;
						$exon_count{$transcript}{$transcript_exon{$transcript}{$rank}[0]} = 0; # Initialize count hash

						# Transcript info
						$transcript_length{$transcript} += $exon_length{$transcript}{$transcript_exon{$transcript}{$rank}[0]};
					}elsif($rank > $new_rank1 && $rank < $new_rank2){ # Intermediate exons
						# Exonic position info
						for(my $i=$transcript_exon{$transcript}{$rank}[1];$i<=$transcript_exon{$transcript}{$rank}[2];$i++){
							$exon_gene_transcript{$i}{$transcript_coding{$transcript}[0]}{$transcript}{$transcript_exon{$transcript}{$rank}[0]} = $transcript_exon{$transcript}{$rank}[3];
						}

						# Exon info
						$exon_length{$transcript}{$transcript_exon{$transcript}{$rank}[0]} = $transcript_exon{$transcript}{$rank}[2] - $transcript_exon{$transcript}{$rank}[1] + 1;
						$exon_count{$transcript}{$transcript_exon{$transcript}{$rank}[0]} = 0; # Initialize count hash

						# Transcript info
						$transcript_length{$transcript} += $exon_length{$transcript}{$transcript_exon{$transcript}{$rank}[0]};
					}
				}
			}
			# Reverse strand
			elsif($transcript_exon{$transcript}{$rank1}[3] == -1){
				# Clip from START CODON
				my $bases = 0;
				foreach my $rank(sort {$a <=> $b} keys %{$transcript_exon{$transcript}}){
					for(my $i=$transcript_exon{$transcript}{$rank}[2];$i>=$transcript_exon{$transcript}{$rank}[1];$i--){
						$bases++;
						$new_start = $i;
						$new_rank1 = $rank;
						last if $bases == $clip_start + 1;
					}
					last if $bases == $clip_start + 1;
				}
				# Clip from STOP CODON
				$bases = 0;
				foreach my $rank(sort {$b <=> $a} keys %{$transcript_exon{$transcript}}){
					for(my $i=$transcript_exon{$transcript}{$rank}[1];$i<=$transcript_exon{$transcript}{$rank}[2];$i++){
						$bases++;
						$new_stop = $i;
						$new_rank2 = $rank;
						last if $bases == $clip_stop + 1;
					}
					last if $bases == $clip_stop + 1;
				}
				foreach my $rank(sort {$a <=> $b} keys %{$transcript_exon{$transcript}}){
					if($rank == $new_rank1){ # First exon
						# Exonic position info
						for(my $i=$transcript_exon{$transcript}{$rank}[1];$i<=$new_start;$i++){
							$exon_gene_transcript{$i}{$transcript_coding{$transcript}[0]}{$transcript}{$transcript_exon{$transcript}{$rank}[0]} = $transcript_exon{$transcript}{$rank}[3];
						}

						# Exon info
						$exon_length{$transcript}{$transcript_exon{$transcript}{$rank}[0]} = $new_start - $transcript_exon{$transcript}{$rank}[1] + 1;
						$exon_count{$transcript}{$transcript_exon{$transcript}{$rank}[0]} = 0; # Initialize count hash

						# Transcript info
						$transcript_length{$transcript} += $exon_length{$transcript}{$transcript_exon{$transcript}{$rank}[0]};
					}elsif($rank == $new_rank2){ # Last exon
						# Exonic position info
						for(my $i=$new_stop;$i<=$transcript_exon{$transcript}{$rank}[2];$i++){
							$exon_gene_transcript{$i}{$transcript_coding{$transcript}[0]}{$transcript}{$transcript_exon{$transcript}{$rank}[0]} = $transcript_exon{$transcript}{$rank}[3];
						}

						# Exon info
						$exon_length{$transcript}{$transcript_exon{$transcript}{$rank}[0]} = $transcript_exon{$transcript}{$rank}[2] - $new_stop + 1;
						$exon_count{$transcript}{$transcript_exon{$transcript}{$rank}[0]} = 0; # Initialize count hash

						# Transcript info
						$transcript_length{$transcript} += $exon_length{$transcript}{$transcript_exon{$transcript}{$rank}[0]};
					}elsif($rank > $new_rank1 && $rank < $new_rank2){ # Intermediate exons
						# Exonic position info
						for(my $i=$transcript_exon{$transcript}{$rank}[1];$i<=$transcript_exon{$transcript}{$rank}[2];$i++){
							$exon_gene_transcript{$i}{$transcript_coding{$transcript}[0]}{$transcript}{$transcript_exon{$transcript}{$rank}[0]} = $transcript_exon{$transcript}{$rank}[3];
						}

						# Exon info
						$exon_length{$transcript}{$transcript_exon{$transcript}{$rank}[0]} = $transcript_exon{$transcript}{$rank}[2] - $transcript_exon{$transcript}{$rank}[1] + 1;
						$exon_count{$transcript}{$transcript_exon{$transcript}{$rank}[0]} = 0; # Initialize count hash

						# Transcript info
						$transcript_length{$transcript} += $exon_length{$transcript}{$transcript_exon{$transcript}{$rank}[0]};
					}
				}
			}
		} # End transcript search

		# Close ensembl db_connection
		$dbh->disconnect();

		###########
		## RIBO_SEQ_DB
		###########
		# Connect to ribo_seq sqlite db
		my $dbh2 = DBI->connect('DBI:SQLite:'.$db_ribo,$user,$pw,
									{ RaiseError => 1},) || die "Database connection not made: $DBI::errstr";

		####
		# For each ribo_read, look if it falls in an exon en make a read_count for each gene as well as for each transcript
		####
        my $query8 = "SELECT chr,strand,start,count FROM $table_ribo WHERE chr = '$chr'";
		my $execute8 = $dbh2->prepare($query8);
		$execute8->execute();

		while(my @result8 = $execute8->fetchrow_array()){
			my $ribo = $result8[2];
			if(defined($exon_gene_transcript{$ribo})){
				foreach my $gene(keys %{$exon_gene_transcript{$ribo}}){
					my $strand = $gene_id{$gene}[1];

					foreach my $transcript(keys %{$exon_gene_transcript{$ribo}{$gene}}){
						if($strand == 1 && $result8[1] == 1){
							$transcript_count{$gene}{$transcript} += $result8[3];
                            $transcript_covered_positions{$gene}{$transcript} += 1;
						}elsif($strand == -1 && $result8[1] == -1){
							$transcript_count{$gene}{$transcript} += $result8[3];
                            $transcript_covered_positions{$gene}{$transcript} += 1;
						}

						foreach my $exon(keys %{$exon_gene_transcript{$ribo}{$gene}{$transcript}}){
							if($strand == 1 && $result8[1] == 1){
								$exon_count{$transcript}{$exon} += $result8[3];
							}elsif($strand == -1 && $result8[1] == -1){
								$exon_count{$transcript}{$exon} += $result8[3];
							}
						}
					}
				}
			}
		}
		$execute8->finish();

        #Get the sequencing depth for calculating FPKM
        my $query9 = "SELECT mapped_T FROM statistics WHERE sample='".$sample_name."' AND type='genomic';";
        my $execute9 = $dbh2->prepare($query9);
        $execute9->execute();
        my $seq_depth = $execute9->fetch()->[0];
        $execute9->finish();
        my $seq_depthM = $seq_depth/1000000; #Sequence depth in millions
        
		####
		# NORMALIZE the transcript-counts and write the information for each chr in a separate file
		####
		## Write all translation info to temporary csv-files
		### Init temporary csv-files
		my $temp_trans_csv = $tmpfolder."fastq1_transcript_".$chr."_tmp.csv";
		open(TMPtrans,'>>'.$temp_trans_csv);

		foreach my $gene(keys %transcript_count){
			foreach my $transcript(keys %{$transcript_count{$gene}}){
				# Transcript_info
				my $transcount = $transcript_count{$gene}{$transcript};
                my $trans_cov_pos = $transcript_covered_positions{$gene}{$transcript};
				my $length = $transcript_length{$transcript};
                my $lengthK = $length/1000; #Length in kilobases
                my $transFPM = $transcount/$seq_depthM; #Fragments per million reads sequence depth
                my $transFPKM = $transFPM/$lengthK; #Fragments per kilo read length and millions sequence depth
                my $trans_cov = (1.0 * $trans_cov_pos)/$length;

				# Now normalize for each transcript separately -> transcript_count/transcript_length
				my $norm_transcount = $transcount/$length;

				# See if every exon of the transcript has a more or less similar RIBO PROFILE COVERAGE
				my $exon_sum = 0;
				my $exon_number = keys %{$exon_count{$transcript}};
				my $exon_premise = "No";
				if($exon_number > 0){
					# Calculate mean normalized exon coverage
					foreach my $exon(keys %{$exon_count{$transcript}}){
						my $exoncount = $exon_count{$transcript}{$exon};
						my $exonlength = $exon_length{$transcript}{$exon};

						# Normalize exon coverage
						my $norm_exoncount = $exoncount/$exonlength;

						# Total coverage overall exons
						$exon_sum += $norm_exoncount;
					}
					my $exon_mean = $exon_sum/$exon_number;

					# Define how many exons fit the threshold
					my $noise = 0;
					foreach my $exon(keys %{$exon_count{$transcript}}){
						my $exoncount = $exon_count{$transcript}{$exon};
						my $exonlength = $exon_length{$transcript}{$exon};
						my $norm_exoncount = $exoncount/$exonlength;
						if($norm_exoncount < $exon_mean/5){
							$noise++;
						}
					}

					# Call transcript only if the amount of 'non-fitted exons' is lower than 15%
					if($noise/$exon_number < 0.15){
						$exon_premise = "Yes";
					}
				}
				# Write to csv-file
				if($transcount > 0){
					print TMPtrans $transcript.",".$transcript_id{$transcript}[0].",".$chr.",".$seq_region.",".$gene_id{$gene}[1].",".$transcript_id{$transcript}[2].",".$transcript_id{$transcript}[3].",".$transcount.",".$norm_transcount.",".$transcript_id{$transcript}[1].",".$exon_premise.",".$transcript_attributes{$transcript}[0].",".$transcript_attributes{$transcript}[1].",".$gene.",".$transFPKM.",".$trans_cov."\n";
				}
			}
		}

		# Close temporary csv-files
		close(TMPtrans);

		# Close ribo_seq db_connection
		$dbh2->disconnect();

		# Close process
		print "\tFinished analyzing chromosome ".$chr."\n";
		$pm->finish;

	} # Close chromosome loop

	# Finish forking
	$pm->wait_all_children;

} # Close sub

sub make_sqlite_dumpfile{
	# Gather all temporary csv-files and create SQLite dumpfile
	print "Gather all temporary csv-files and create SQLite dumpfile...\n";
	## Transcripts
	my $temp_trans_csv_all = $tmpfolder."fastq1_transcript_all.csv";
	system("touch ".$temp_trans_csv_all);

	## Put information in dumpfile and remove temporary chr_csv-files
	foreach my $chr(@ch){
		# Transcripts
		my $temp_trans_csv = $tmpfolder."fastq1_transcript_".$chr."_tmp.csv";
		system("cat ".$temp_trans_csv." >> ".$temp_trans_csv_all);
        system("rm -rf ".$temp_trans_csv);
	}

} # Close sub

sub import_sqlite_dumpfile{
	print "Import SQLite dumpfile and remove tmp files...\n";

	my $temp_trans_csv_all = $_[0];

	# Dump transcript import in SQLite transcript_translation_table
	system("/usr/bin/sqlite3 -separator , ".$db_ribo." \".import ".$temp_trans_csv_all." ".$table_ribo_trans." \"");

	# Remove temporary SQLite dump-files

    system("rm -rf ".$temp_trans_csv_all);

	print "DONE!\n";
} # Close sub

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

	$query = "select value from arguments where variable = \'nr_of_cores\'";
    $sth = $dbh_results->prepare($query);
	$sth->execute();
	my $cores = $sth->fetch()->[0];
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
    
    $query = "select value from arguments where variable = \'unique\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $uniq = $sth->fetch()->[0];
    $sth->finish();

	$dbh_results -> disconnect();

	# Return ARG variables
	return($species,$version,$IGENOMES_ROOT,$cores,$run_name,$mapper,$uniq);
} # Close sub

sub store_ENS_var{
	# Catch
	my $db_ribo = $_[0];
	my $user = $_[1];
	my $pw = $_[2];
	my $db_ensembl = $_[3];

    # Connect to db
    my $dbh_results = DBI->connect('DBI:SQLite:'.$db_ribo,$user,$pw,
    							{ RaiseError => 1},) || die "Database connection not made: $DBI::errstr";
    
    #Delete possible previous ensembl db in arguments table
    my $del_query = "DELETE FROM arguments WHERE variable='ens_db';";
    $dbh_results->do($del_query);

	# Store path of ENSEMBL db in arguments table
	my $query = "INSERT INTO arguments (variable,value) VALUES (\'ens_db\',\'".$db_ensembl."\')";
    $dbh_results->do($query);

	$dbh_results -> disconnect();
}

#Print help text
sub print_help_text{
    
    my $help_string = "\n\nTranslated transcript calling Proteoformer
    
This tool calls the translated transcripts out of the aligned RIBO-seq data based on an exon-coverage rule. Transcript without RIBO-seq counts are ignored from the start. Then, for each exon of the transcript, the counts of ribosome reads are calculated and normalized by exon length. If the normalized count of an exon is 5 times lower than the mean normalized exon count of the given transcript, the exon is classified as a noise exon. Transcripts with less than 15% noise exons, are called as actively-translated transcripts (exon_coverage = Yes in the output table).
    
Example:
    perl ribo_translation.pl --in_sqlite SQLite/results.db --out_sqlite SQLite/results.db --ens_db ENS_hsa_92.db

Input parameters:
    --work_dir                              The working directory (default: CWD env setting)
    --tmp                                   The temporary files folder (default: work_dir/tmp)
    --in_sqlite                             The SQLite results database from previous steps (default: SQLite/results.db)
    --out_sqlite                            The SQLite results database to store all results in (default: the same as in_sqlite argument)
    --ens_db                                The Ensembl database with annotation info (mandatory)
    --help                                  Generates this help message
    
";
    
    print $help_string;
    
    return
}
