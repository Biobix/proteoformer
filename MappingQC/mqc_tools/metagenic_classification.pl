#!/usr/bin/perl -w

# This script determines the amount of ribo_seq reads that lay in exons, introns and UTRs
# Output: tabular tables with counts per region and accompanying pie charts
###########

################
# COMMAND LINE:
# ./metagenic_classification.pl --cores 3 --treated untreated (or 'treated') (--work_dir getcwd --in_sqlite SQLite/results.db --output_file run_name_annotation)
# **example: ./metagenic_classification.pl --cores 3 --treated untreated (--work_dir /data/RIBO_runs/RIBO_Ingolia_GerbenM/ --in_sqlite SQLite/results.db --output_file mESC_GA_ens72_STAR_untreated_annotation)
# GALAXY
# metagenic_classification.pl --tooldir "${__tool_directory__}" --cores "${cores}" --treated "${treatment}" --in_sqlite "${sqlite_db}" --out_table1 "${out_table1}" --out_table2 "${out_table2}" --out_pdf1 "${out_pdf1}" --out_pdf2 "${out_pdf2}"
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
my ($tooldir,$cores,$work_dir,$in_sqlite,$output_folder,$treated, $ens_db);
GetOptions(
"tooldir:s"=>\$tooldir,                 # The directory the tool currently resides in (included for Galaxy) optional argument (default = $CWD env setting: for script version)
"cores=i"=>\$cores,                     # Number of cores to use 										mandatory argument
"work_dir:s" =>\$work_dir,              # Working directory												optional  argument (default = $CWD env setting)
"in_sqlite:s" =>\$in_sqlite,            # Results db (relative of CWD)                              	optional  argument (default = SQLite/results.db)
"ens_db:s" =>\$ens_db,                  # ENSEMBL db (relative of CWD)                              	mandatory argument
"output_folder:s" =>\$output_folder,		#Output folder                                              optional  argument (default = CWD/output)
#"igenomes_root:s" =>\$IGENOMES_ROOT,    # Igenomes root folder                                      	mandatory argument
"treated:s" => \$treated 				# Define data('no':no treatment/CHX or 'yes':treated(LTM/Harr))	mandatory argument
);

###########################################################################
#Check all input variable and/or get default values and set extra variables
###########################################################################
# Input variables
if($cores){
    print "Number of cores						: $cores\n";
}else{
    die "\nDon't forget to pass number of cores to use for mapping using the --cores argument!\n\n";
}
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
if ($ens_db){
    print "The Ensembl database     : $ens_db/n";
} else {
    print "Error: do not forget the Ensembl db!";
    die;
}


# Get ARG vars
my $user = "";
my $pw = "";
my $db_ribo= $in_sqlite;
my ($species,$version,$run_name,$mapper,$IGENOMES_ROOT,$unique) = get_ARG_vars($db_ribo,$user,$pw);


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
if ($output_folder){
    print "The output folder is set to : ".$output_folder;
} else {
    $output_folder = $work_dir."/mappingQC/";
    print "The output folder is set to : ".$output_folder;
}
my $out_dir = $output_folder;
if(!-d "$out_dir"){
    system ("mkdir -p ".$out_dir);
}

my $output_file = $run_name."_ens".$version."_".$mapper."_".$treated;
my $out_table1 = $out_dir."/".$output_file."_annotation_coding.txt";
my $out_table2 = $out_dir."/".$output_file."_annotation_noncoding.txt";
my $out_png1 = $out_dir."/".$output_file."_annotation_coding.png";
my $out_png2 = $out_dir."/".$output_file."_annotation_noncoding.png";
print "Output Table coding transcripts			    	: $out_table1\n";
print "Pie Chart for coding transcripts			: $out_png1\n";
print "Output Table non-coding transcripts	    		: $out_table2\n";
print "Pie Chart for non-coding transcripts			: $out_png2\n";

# Get chromosomes and correct coord_system_id
print "Get chromosomes and coord_system_id...\n";
#my $chromosome_sizes;
my $coord_system_id; my @ch;

#$chromosome_sizes = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Annotation/Genes/ChromInfo.txt";
@ch = @{get_chr($species)};
#@ch = (1,19);
$coord_system_id = get_coord_system_id($db_ensembl,$assembly);

# Get non protein-coding Ensembl protein biotypes
my %biotypes = %{get_nPCbiotypes($db_ensembl,$user,$pw)};

################
## METAGENIC CLASSIFICATION OF RIBO-SEQ DATA
################
# Open files
open OUT1,"+>>".$out_table1 or die $!;
open OUT2,"+>>".$out_table2 or die $!;

print OUT1 "chr\tribo\texon\t5utr\t3utr\tintron\tnon_protein_coding\tintergenic\n";
print OUT2 "chr\tnon_protein_coding";
foreach my $biotype(sort keys %biotypes){
	print OUT2 "\t".$biotype;
}
print OUT2 "\n";

# Analysis
metagenic_analysis($db_ribo,$db_ensembl,$table_ribo,\@ch,$cores,$coord_system_id,$user,$pw);

# Make Pie Charts
piecharts($out_table1,$out_table2,$out_png1,$out_png2,$tooldir);

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

	$query = "select value from arguments where variable = \'run_name\'";
    $sth = $dbh_results->prepare($query);
	$sth->execute();
	my $run_name = $sth->fetch()->[0];
	$sth->finish();

    $query = "select value from arguments where variable = \'unique\'";
    $sth = $dbh_results->prepare($query);
    $sth->execute();
    my $unique = $sth->fetch()->[0];
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

	$dbh_results -> disconnect();

	# Return ENS variables
	return($species,$version,$run_name,$mapper,$IGENOMES_ROOT,$unique);

} # Close sub

sub get_chr{
    # Catch
    my $species = $_[0];

    # Catch
    my %chr_sizes;
    my $filename = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Annotation/Genes/ChromInfo.txt";
    print "$filename\n";
    if (-e $filename) {
        # Work
        open (Q,"<".$filename) || die "Cannot open chr sizes input $filename\n";
        while (<Q>){
            my @a = split(/\s+/,$_);
            $chr_sizes{$a[0]} = $a[1];
        }
    }
    else
    {
        $filename = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Sequence/WholeGenomeFasta/GenomeSize.xml";
        open (Q,"<".$filename) || die "Cannot open chr sizes input $filename\n";
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

sub get_nPCbiotypes{

	# Catch
	my $db_ensembl = $_[0];
	my $user = $_[1];
	my $pw = $_[2];

	# Connect to Ensembl db
	my $dbh = DBI->connect('DBI:SQLite:'.$db_ensembl,$user,$pw,
					{RaiseError => 1},) || die "Database connection not made: $DBI::errstr";

	# Query db
	my $query = "SELECT biotype FROM transcript WHERE biotype NOT LIKE '%protein_coding%' GROUP BY biotype";
	my $execute = $dbh->prepare($query);
	$execute->execute();

	my %biotypes;
	while(my @results = $execute->fetchrow_array()){
		$biotypes{$results[0]} = 0;
	}

	$execute->finish();

	# Disconnect
	$dbh->disconnect();

	# Return
	return(\%biotypes)
}

sub metagenic_analysis{
	print "Metagenic Classification...\n";

	# Catch
	my $db_ribo = $_[0];
	my $db_ensembl = $_[1];
	my $table_ribo = $_[2];
	my @ch = @{$_[3]};
	my $cores = $_[4];
	my $coord_system_id = $_[5];
	my $user = $_[6];
	my $pw = $_[7];

	# Init multi core
	my $processes = $cores; # Nr of processes
	my $pm = new Parallel::ForkManager($processes); # Open fork

	# Loop through chromosomes
	foreach my $chr(@ch){
		# Start parallel process
		$pm->start and next;

		print "\tStarting analysis for chr ".$chr."...\n";

		###########
		## ENSEMBL DB
		###########
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

		###########
		## TRANSCRIPTs
		###########
		# Init
		my $trs_c = {};
		my $trs_nc = {};

		# Get all protein_coding transcripts
		my $query1 = "SELECT transcript_id,gene_id,seq_region_start,seq_region_end,seq_region_strand,biotype,stable_id FROM transcript WHERE seq_region_id = '$seq_region' AND biotype = 'protein_coding'";
		my $execute1 = $dbh->prepare($query1);
		$execute1->execute();

		$trs_c = $execute1->fetchall_hashref('transcript_id');
		$execute1->finish();


		# Get all other transcripts
		my $query2 = "SELECT transcript_id,gene_id,seq_region_start,seq_region_end,seq_region_strand,biotype,stable_id FROM transcript WHERE seq_region_id = '$seq_region' AND biotype NOT LIKE '%protein_coding%'";
		my $execute2 = $dbh->prepare($query2);
		$execute2->execute();

		$trs_nc = $execute2->fetchall_hashref('transcript_id');
		$execute2->finish();

		# Split transcripts in forward and reverse arrays
		my ($tr_c_for,$tr_c_rev);
		foreach my $tr_id (sort { $trs_c->{$a}{'seq_region_start'} <=> $trs_c->{$b}{'seq_region_start'} } keys %{$trs_c}){
			if($trs_c->{$tr_id}{'seq_region_strand'} eq '1'){	# Forward
				push(@$tr_c_for,$tr_id);
			}else{	# Reverse
				push(@$tr_c_rev,$tr_id);
			}
		}

		my ($tr_nc_for,$tr_nc_rev);
		foreach my $tr_id (sort { $trs_nc->{$a}{'seq_region_start'} <=> $trs_nc->{$b}{'seq_region_start'} } keys %{$trs_nc}){
			if($trs_nc->{$tr_id}{'seq_region_strand'} eq '1'){	# Forward strand
				push(@$tr_nc_for,$tr_id);
			}else{	# Reverse strand
				push(@$tr_nc_rev,$tr_id);
			}
		}

		###########
		## EXONs & UTRs
		###########
		## PROTEIN-CODING TRANSCRIPTS
		###########
		# Init
		my $exon = {};
		my (%cds_for,%cds_rev); # CDS = UTRs + EXONs
		foreach my $tr_id(keys %{$trs_c}){
			# Get all exons for obtained transcripts from tables exon and exon_transcript
			my $query3 = "SELECT a.exon_id,b.seq_region_start,b.seq_region_end,b.seq_region_strand FROM exon_transcript a JOIN exon b ON a.exon_id = b.exon_id WHERE a.transcript_id = '$tr_id'";
			my $execute3 = $dbh->prepare($query3);
			$execute3->execute();
			$exon = $execute3->fetchall_hashref('exon_id');

			$execute3->execute();
			while(my @result3 = $execute3->fetchrow_array()){
				#$result3[0]: exon_id
				#$result3[1]: seq_region_start
				#$result3[2]: seq_region_end
				#$result3[3]: seq_region_strand

				for(my $i1=$result3[1];$i1<=$result3[2];$i1++){
					if($result3[3] == 1){	# Forward strand
						$cds_for{$chr.':'.$i1}{'exon'}++;
					}else{	# Reverse strand
						$cds_rev{$chr.':'.$i1}{'exon'}++;
					}
				}
			}
			$execute3->finish();

			# Get first and last exon of each transcript from table translation
			my $query4 = "SELECT start_exon_id,end_exon_id,seq_start,seq_end FROM translation WHERE transcript_id = '$tr_id'";
			my $execute4 = $dbh->prepare($query4);
			$execute4->execute();

			while(my @result4 = $execute4->fetchrow_array()){
				#$result4[0]: start_exon_id
				#$result4[1]: end_exon_id
				#$result4[2]: seq_start (offset position translation start)
				#$result4[3]: seq_end (offset position translation stop)

				my ($start_id,$end_id,$seq_start,$seq_end) = ($result4[0],$result4[1],$result4[2],$result4[3]);
				my ($start_first_exon,$stop_first_exon,$start_last_exon,$stop_last_exon);
		 		my ($start_codon,$stop_codon);

		 		# Determine start codon and stop codon
		 		# Determine 5' and 3' UTRs

		 		if($trs_c->{$tr_id}{'seq_region_strand'} eq '1'){ # Forward strand
		 			$start_first_exon = $exon->{$start_id}{'seq_region_start'};
		 			$stop_first_exon = $exon->{$start_id}{'seq_region_end'};
		 			$start_last_exon = $exon->{$end_id}{'seq_region_start'};
		 			$stop_last_exon = $exon->{$end_id}{'seq_region_end'};
		 			$start_codon = $start_first_exon + $seq_start - 1;
		 			$stop_codon = $start_last_exon + $seq_end - 1;

		 			for(my $l1=$start_first_exon;$l1<$start_codon;$l1++){
		 				$cds_for{$chr.':'.$l1}{'5UTR'}++;
		 			}
		 			for(my $l2=($stop_codon+1);$l2<=$stop_last_exon;$l2++){
		 				$cds_for{$chr.':'.$l2}{'3UTR'}++;
		 			}
		 		}elsif($trs_c->{$tr_id}{'seq_region_strand'} eq '-1'){ # Reverse strand
		 			$start_first_exon = $exon->{$start_id}{'seq_region_end'};
		 			$stop_first_exon = $exon->{$start_id}{'seq_region_start'};
		 			$start_last_exon = $exon->{$end_id}{'seq_region_end'};
		 			$stop_last_exon = $exon->{$end_id}{'seq_region_start'};
		 			$start_codon = $start_first_exon - $seq_start + 1;
		 			$stop_codon = $start_last_exon - $seq_end + 1;

		 			for(my $l3=($start_codon+1);$l3<=$start_first_exon;$l3++){
		 				$cds_rev{$chr.':'.$l3}{'5UTR'}++;
		 			}for(my $l4=$stop_last_exon;$l4<$stop_codon;$l4++){
		 				$cds_rev{$chr.':'.$l4}{'3UTR'}++;
		 			}
		 		}
			}
			$execute4->finish();
		}

		###########
		## NON-CODING TRANSCRIPTS
		###########
		# Get biotypes
		my $query5 = "SELECT biotype FROM transcript WHERE biotype NOT LIKE '%protein_coding%' GROUP BY biotype";
		my $execute5 = $dbh->prepare($query5);
		$execute5->execute();

		my %biotypes_nc;
		while(my @results5 = $execute5->fetchrow_array()){
			$biotypes_nc{$results5[0]} = 0;
		}
		$execute5->finish();
		# Disconnect
		$dbh->disconnect();

		###########
		## RIBO-DATA
		###########
		###########
		## ANNOTATE RIBO-SEQ READS
		###########
		print "\t\tAnnotating ribo-seq reads of chr ".$chr."\n";

		# Get ribo-seq reads, split per strand
		my $ribo_for = get_reads($db_ribo,$user,$pw,$table_ribo,$chr,1);
		my $ribo_rev = get_reads($db_ribo,$user,$pw,$table_ribo,$chr,-1);

		# Init values
		my ($ribo_reads,$ribo_readsnc) = (0,0);
 		my ($ribo_exon,$ribo_intron,$ribo_5utr,$ribo_3utr,$ribo_intergenic) = (0,0,0,0,0);

		# Loop over forward ribo-seq reads
		my @window_c = (); # Init window with protein-coding transcripts
		my @window_nc = (); # Init window with non protein-coding transcripts
		foreach my $pos (sort {$a <=> $b} keys %{$ribo_for}){
			######
			## PROTEIN-CODING WINDOW
			######
			# Push all tr_ids to @window_c where tr_start < window_pos
			foreach my $tr_for_id (@$tr_c_for){
				if($trs_c->{$tr_for_id}{'seq_region_start'} <= $pos){
					push(@window_c,$tr_for_id);
				}else{last;} # Don't unnecessarily loop over all transcripts
			}

			# Get rid of tr_c_for elements already in @window_c
			@$tr_c_for = grep {$trs_c->{$_}{'seq_region_start'} > $pos} @$tr_c_for;

			# Get rid of tr_ids in @window_c where tr_end < window_pos
			@window_c = grep {$trs_c->{$_}{'seq_region_end'} >= $pos} @window_c;

			#####
			## NON-CODING WINDOW
			#####
			# Push all tr_ids to @window_nc where tr_start < window_pos
			foreach my $tr_for_id (@$tr_nc_for){
				if($trs_nc->{$tr_for_id}{'seq_region_start'} <= $pos){
					push(@window_nc,$tr_for_id);
				}else{last;} # Don't unnecessarily loop over all transcripts
			}

			# Get rid of tr_nc_for elements already in @window_nc
			@$tr_nc_for = grep {$trs_nc->{$_}{'seq_region_start'} > $pos} @$tr_nc_for;

			# Get rid of tr_ids in @window_nc where tr_end < window_pos
			@window_nc = grep {$trs_nc->{$_}{'seq_region_end'} >= $pos} @window_nc;

			#####
			## ANNOTATE
			#####
			# Loop over windows and check functional annotation of each ribo-seq read
			$ribo_reads++;
			if(@window_c){
				# Annotate reads in PROTEIN-CODING transcripts
				my $count = $ribo_for->{$pos}{'count'};

				# Check 5'UTR
				if(defined($cds_for{$chr.':'.$pos}{'5UTR'})){
					$count = 0;
					$ribo_5utr++;
				}

				# Check 3'UTR
				if($count > 0 && defined($cds_for{$chr.':'.$pos}{'3UTR'})){
					$count = 0;
					$ribo_3utr++;
				}

				# Check Exons
				if($count > 0 && defined($cds_for{$chr.':'.$pos}{'exon'})){
					$count = 0;
					$ribo_exon++;
				}

				# If still not defined -> intronic region
				if($count > 0){
					$ribo_intron++;
				}
			}elsif(@window_nc){
				# Annotate reads in NON PROTEIN-CODING transcripts
				$ribo_readsnc++;

				# Define biotype (if #biotypes>0, take a random/first one)
				my @biotypes; my $i = 0;
				foreach my $tr_for_id (@window_nc){
						$biotypes[$i] = $trs_nc->{$tr_for_id}{'biotype'};
						$i++;
				}
				my $random = int(rand(@biotypes));
				my $biotype = $biotypes[$random];
				$biotypes_nc{$biotype}++;
			}else{
				$ribo_intergenic++;
			}
		} # Close forward-loop

		# Loop over reverse ribo-seq reads
		@window_c = (); # Empty window_c
		@window_nc = (); # Empty window_nc
		foreach my $pos (sort {$a <=> $b} keys %{$ribo_rev}){
			######
			## PROTEIN-CODING WINDOW
			######
			# Push all tr_ids to @window_c where tr_start < window_pos
			foreach my $tr_rev_id (@$tr_c_rev){
				if($trs_c->{$tr_rev_id}{'seq_region_start'} <= $pos){
					push(@window_c,$tr_rev_id);
				}else{last;} # Don't unnecessarily loop over all transcripts
			}

			# Get rid of tr_c_rev elements already in @window_c
			@$tr_c_rev = grep {$trs_c->{$_}{'seq_region_start'} > $pos} @$tr_c_rev;

			# Get rid of tr_ids in @window_c where tr_end < window_pos
			@window_c = grep {$trs_c->{$_}{'seq_region_end'} >= $pos} @window_c;

			######
			## NON-CODING WINDOW
			######
			# Push all tr_ids to @window_nc where tr_start < window_pos
			foreach my $tr_rev_id (@$tr_nc_rev){
				if($trs_nc->{$tr_rev_id}{'seq_region_start'} <= $pos){
					push(@window_nc,$tr_rev_id);
				}else{last;} # Don't unnecessarily loop over all transcripts
			}

			# Get rid of tr_nc_for elements already in @window_nc
			@$tr_nc_rev = grep {$trs_nc->{$_}{'seq_region_start'} > $pos} @$tr_nc_rev;

			# Get rid of tr_ids in @window_nc where tr_end < window_pos
			@window_nc = grep {$trs_nc->{$_}{'seq_region_end'} >= $pos} @window_nc;

			######
			## ANNOTATE
			######
			# Loop over windows and check functional annotation of each ribo-seq read
			$ribo_reads++;
			if(@window_c){
				# Annotate reads in PROTEIN-CODING transcripts
				my $count = $ribo_rev->{$pos}{'count'};

				# Check 5'UTR
				if(defined($cds_rev{$chr.':'.$pos}{'5UTR'})){
					$count = 0;
					$ribo_5utr++;
				}

				# Check 3'UTR
				if($count > 0 && defined($cds_rev{$chr.':'.$pos}{'3UTR'})){
					$count = 0;
					$ribo_3utr++;
				}

				# Check Exons
				if($count > 0 && defined($cds_rev{$chr.':'.$pos}{'exon'})){
					$count = 0;
					$ribo_exon++;
				}

				# If still not defined -> intronic region
				if($count > 0){
					$ribo_intron++;
				}
			}elsif(@window_nc){
				# Annotate reads in NON PROTEIN-CODING transcripts
				$ribo_readsnc++;

				# Define biotype (if #biotypes>0, take a random/first one)
				my @biotypes; my $i = 0;
				foreach my $tr_rev_id (@window_nc){
						$biotypes[$i] = $trs_nc->{$tr_rev_id}{'biotype'};
						$i++;
				}
				my $random = int(rand(@biotypes));
				my $biotype = $biotypes[$random];
				$biotypes_nc{$biotype}++;

			}else{
				$ribo_intergenic++;
			}
		} # Close rev-loop

		# Print results
		print OUT1 $chr."\t".$ribo_reads."\t".$ribo_exon."\t".$ribo_5utr."\t".$ribo_3utr."\t".$ribo_intron."\t".$ribo_readsnc."\t".$ribo_intergenic."\n";
		print OUT2 $chr."\t".$ribo_readsnc;
		foreach my $biotype(sort keys %biotypes_nc){
			print OUT2 "\t".$biotypes_nc{$biotype};
		}
		print OUT2 "\n";

		# Close process
		$pm->finish;

	} # Close chr-loop

	# Finish forking
	$pm->wait_all_children;

} # Close sub

sub get_reads{

	# Catch
	my $db_ribo = $_[0];
	my $user = $_[1];
	my $pw = $_[2];
	my $table_ribo = $_[3];
	my $chr = $_[4];
	my $strand = $_[5];

	# Init
	my $ribo_reads = {};

	# Connect to ribo_seq database
	my $dbh2 = DBI->connect('DBI:SQLite:'.$db_ribo,$user,$pw,
						{RaiseError => 1},) || die "Database connection not made: $DBI::errstr";

	# Query
	my $query = "SELECT * FROM $table_ribo WHERE chr = '$chr' and strand = '$strand'";
	my $execute = $dbh2->prepare($query);
	$execute->execute();

	$ribo_reads = $execute->fetchall_hashref('start');

	$execute->finish();

	# Disconnect
	$dbh2->disconnect();

	# Return
	return($ribo_reads);
}

sub piecharts{
	print "Make Pie Charts...\n";

	# Catch
	my $out_table1 = $_[0];
	my $out_table2 = $_[1];
	my $out_png1 = $_[2];
	my $out_png2 = $_[3];
	my $tooldir = $_[4];

	# Execute Rscript
    system("Rscript ".$tooldir."/metagenic_piecharts.R ".$out_table1." ".$out_table2." ".$out_png1." ".$out_png2);
    #system("Rscript metagenic_piecharts.R ".$out_table1." ".$out_table2." ".$out_pdf1." ".$out_pdf2);
}
