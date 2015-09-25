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


use strict;
use warnings;
use POSIX ":sys_wait_h";
use Getopt::Long;
use Cwd;
use Cwd 'abs_path';
use Bio::SearchIO;
use Bio::SeqIO;
use DBI;
use LWP::UserAgent; 
use XML::Smart;
use Parallel::ForkManager;

# ---------------------------------------------------------------------
	##	GLOBAL VARIABLES

my $blastdb;		# Usearch/blastp formatted database
my $result_db;		# SQLite source database
my $tis_ids;		# TIS ID to generate table from
my $mapping;		# Ensembl database name to download biomart mapped transcripts 
my $mflag;			# Flag for biomart mappings, 1 = remote download, 2 = local file, 3 = sequence based mapping, 4 = no mapping.
my $total_tr = 0;
my $no_blast_tr = 0;
my $num_non_Red_trans = 0;

# ---------------------------------------------------------------------
	##	Other prameters
my $user = "";		# User name for sqlite database
my $password = "";	# Password for sqlite database
my $work_dir;		# Working directory
my $blast_pgm;		# The program to use for search, Default program = blastp [NB: Not valid when mflag = 4]
my $min_blast_length;	# minimum sequence length to perform blast
my $coverage;		# Minimum number of identical positions
my $mslength;		# Minimum transcript length to allow in the translation product database
my $evalue;			# blast e-value
my $identity;		# Minimum alignment score
my $gapopen;		# Gap opening penalty
my $gapextend;		# Gap extension penalty
my $matrix;			# Blast search matrix
my $word_size;		# word size
my $extra_info;		# Extra refererence
my $translation_db;	# FASTA file of non redundant derived translation products
my $tis_call;		# Allow annotated TIS that do not pass the TIS calling algorithm in the databases [Y or N]
my $db_config_version; 	# Ensembl databset confirguration version
my $external_ref;	# External reference in biomart to map transcripts to
my $tmp;				# temporary folder

# ---------------------------------------------------------------------
	##	GET command line arguments
GetOptions(
	'password=s'			=> \$password,
	'user=s'				=> \$user,
	'blast_db=s' 	 		=> \$blastdb,
	'result_db=s' 	 		=> \$result_db,
	'tis_ids=s' 			=> \$tis_ids,
	'work_dir=s' 	 		=> \$work_dir,
	'blast_pgm=s' 	 		=> \$blast_pgm,
	'evalue=f' 	 			=> \$evalue,
	'min_blast_length=i'	=> \$min_blast_length,
	'mslenght=i'	 		=> \$mslength,
	'mflag=i'	 			=> \$mflag,
	'identity=f'			=> \$identity,
	'coverage=i'			=> \$coverage,
	'mapping=s'			=> \$mapping,
	'db_config_version=f'		=> \$db_config_version,
	'external_ref=s'		=> \$external_ref,
	'extra_info=s'			=> \$extra_info,
	'gapopen=f'	 			=> \$gapopen,
	'gapextend=f'	 		=> \$gapextend,
	'word_size=i'	 		=> \$word_size,
	'matrix=s'	 			=> \$matrix,
	'translation_db=s'		=> \$translation_db,
	'tis_call=s'			=> \$tis_call,
	'tmp=s'					=> \$tmp
);

# ---------------------------------------------------------------------
	## EXECUTION

## Command line
# perl generate_translation_db.pl -blast_db /path/to/blast_db -result_db results.db -tis_ids 1 -blast_pgm ublast -mflag 0 -external_ref uniprot_swissprot_accession -mapping_db mmusculus_gene_ensembl -num_threads 3 -work_dir working_dir -tis_call Y -tmp temporary_dir
# mmusculus_gene_ensembl
# mart_export_mm.txt

# perl generate_translation_db.pl -blast_db /path/to/blast_db -result_db results.db -tis_ids 1 -blast_pgm ublast -mflag 1 -mapping_db mart_export_mm.txt -num_threads 3 -work_dir working_dir -tis_call Y -tmp temporary_dir -translation_db file_name
# mmusculus_gene_ensembl
# mart_export_mm.txt
#####


my $CWD = getcwd();
if (!($work_dir)) {
	$work_dir = $CWD;
} elsif (!-d "$work_dir") { 
	system ("mkdir ". $work_dir);
} 

print "\n";
if ($work_dir) {
	$work_dir = abs_path($work_dir);
	print "The following working directory : $work_dir\n";
}

my $TMP  = ($ENV{'TMP'}) ? $ENV{'TMP'} : ($tmp) ? $tmp : "$work_dir/tmp"; # (1) get the TMP environment variable, (2) get the $tmpfolder variable, (3) get current_working_dir/tmp
if (!-d "$TMP") { 
	system ("mkdir ".$TMP);
} 
print "The following tmpfolder is used: $TMP\n";

if ($result_db) {
	print "SQLite database containing transcripts : $result_db\n";
} 
if ($tis_ids) {
	print "Database is generated for TIS ids : $tis_ids\n";
}

if ($tis_call) {
	print "Allow annotated transcripts that do not pass the TIS calling algorithm: $tis_call\n";
} else {
	$tis_call = "Y";
	print "Allow annotated transcripts that do not pass the TIS calling algorithm: $tis_call\n";
}

if ($mflag) {

	if ($mflag == 1) {
		if ($mapping) {
			print "Ensembl database for remote mapping : $mapping\n";
			if ($db_config_version) {
				print "Ensembl dataset configuration version : $db_config_version\n";
			} else {
				print "Dataset configuration version is require.\n";
				exit;
			}
			
			if ($external_ref) {
			
				print "External reference to map transcript with : $external_ref\n";
				if ($extra_info) {
					print "External reference to map transcript with : $extra_info\n";
				}
				
			} else {
				print "Ensembl external attribute for biomart mapping is required.\n";
				exit;
			}

		} else {
			print "Mapping database require!\n";
			print "If mflag = 1, the name of the biomart database to download mapping to external Id is required.\n";
			exit;
		}
		
	} elsif ($mflag == 2) {
		if ($mapping) {
			if (-e $mapping) {
				print "Locale file containing ensembl to Swissprot mapping : $mapping\n";
			} else {
				print "No such file $mapping.\tEnsure the file exist and you have the required permission.\n";
				exit;
			}
		} else {
			print "If mflag = 2, a comma seperated file for the Transcript id to external reference file is require.\n";
			print "See the Readme file more information.\n";

		}
	} elsif ($mflag == 3) {
	
		print "Sequence based mapping of transcripts to canonical database by blast search.\n";
		
		if ($blast_pgm) {
			print "Blast program used for mapping : $blast_pgm\n"; 
		} else {
			$blast_pgm = "ublast";
			print "Blast program used for mapping : $blast_pgm\n"; 
		}
		
	} elsif ($mflag == 4) {
		print "Derived translation product  database will not be mapped to any canonical information.\n"
	}

} else {
	print "Derived translation product database will not be mapped to any canonical information.\n";
}

if ($blastdb) {
	print "The blast database is : $blastdb\n";
} else {
	print "No blast database supplied. Ensure you have choose the no blast search option.\n"
}

if ($evalue) {
	print "Blast e-value : $evalue\n";
} else {
	$evalue=1e-10;
	print "Blast e-value : $evalue\n";
}

if ($min_blast_length) {
	print "Minimum sequence length allowed for Blast search	: $min_blast_length\n";
} else {
	$min_blast_length = 32;
	print "Minimum sequence length allowed for Blast search : $min_blast_length\n";	
}

if ($identity) {
	print "Blast Identity value : $identity%\n";
} else {
	$identity=75;
	print "Blast Identity value : $identity%\n";
}

if ($coverage) {
	print "Minimum percentage of identical positions : $coverage\n";
} else {
	$coverage=30;
	print "Minimum percentage of identical positions : $coverage\n";
}

if ($word_size) {
	print "Minimum Word size: $word_size\n";
} else {
	$word_size = 3;
	print "Minimum Word size: $word_size\n";
}

if ($gapopen) {
	print "Cost of gap open	: $gapopen\n";
} else {
	$gapopen = 11;
	print "Cost of gap open	: $gapopen\n";
}

if ($mslength) {
	print "Minimum Word size: $mslength\n";
} else {
	$mslength=6;
	print "Minimum sequence length	: $mslength\n";
}

if ($gapextend) {
	print "Gap extension penalty	: $gapextend\n";
} else {
	$gapextend = 1;
	print "Gap extension penalty	: $gapextend\n";
}

if ($matrix) {
	print "Minimum Word size: $matrix\n";
} else {
	$matrix="BLOSUM62";
	print "Matrix blast search matrix : $matrix\n";
}

	# get arguments from SQLite DB
my $dsn_results = "DBI:SQLite:dbname=$result_db";
my $dbh_results = dbh($dsn_results,$user,$password);
my ($run_name,$species,$mapper,$nr_of_cores)=get_input_vars($dbh_results);
print "Number of cores used: $nr_of_cores\n";
print "\n";

my %annotation = ( "aTIS"=>1,"5UTR"=>2,"CDS"=>3,"ntr"=>4,"3UTR"=>5);
my %top_anno = ("aTIS","5UTR");

print "Annotations and their rankings (with 1 the most important) \n";
foreach my $key (sort {$annotation{$a} <=> $annotation{$b}} keys %annotation) {
    print "$key\t$annotation{$key}\n";
}

	# If Mapping to canonical database is allowed 
my $external_mapping = {};
if ($mflag) {

	if ($mflag == 1) {	# Id based mapping
	
		my $mapping_file = remote_biomart_mapping($mapping,$db_config_version,$external_ref,$extra_info);
		print "File containing maping information :	$mapping_file\n";
		$external_mapping = biomart_mapping($mapping_file);

	} elsif ($mflag == 2)  {
		$external_mapping = biomart_mapping($mapping);
	} 
}

my @idsref = get_analysis_ids($dbh_results,$tis_ids);  #$tis_ids is input variable

foreach (@idsref) {	# generate translation db for selected tis_ids 

	generate_trans_db($_);
}



#############################################
	#
	# SUBS
	#
#############################################

sub generate_trans_db {

	my $tis_id = shift;

	my $startRun = time();
	my $table = "tis_".$tis_id."_transcripts";
	my @tis = split '_', $tis_id;
	
	my ($transcript,$gene_transcript) = get_transcripts_from_resultdb($dbh_results,$table,$tis[0]);
	$total_tr = scalar(keys %$transcript);
	print "total number of TIS identified ".$total_tr."\n";
	print "Total number of genes with one or more identified TIS ".scalar(keys %$gene_transcript)."\n";

	$transcript = remove_redundancy($transcript,$gene_transcript);
	my $total_non_red_tr = scalar(keys %$transcript);
	print "Total number of TIS after removing redundancy ".$total_non_red_tr."\n";
	
	# 1 = remote download, 2 = local file, 3 = sequence based mapping, 4 = no mapping.
	if ($mflag == 1 or $mflag == 2) {		#  ID based mapping
		$transcript = id_based_mapping($transcript,$external_mapping);
		
	} elsif ($mflag == 3) {		# Sequence based mapping
		
		my $file_to_blast = transcripts_to_blast($transcript, $tis_id);
		if ($blast_pgm eq "blastp") {
			$transcript = blastp($transcript, $file_to_blast);
		} elsif ($blast_pgm eq "ublast")  {
			$transcript = ublast($transcript, $file_to_blast);
		}
		
	}
	
	my $output_dir = $work_dir."/derived_db";
	unless (-d "$output_dir") { system ("mkdir ".$output_dir)}
	$translation_db =  path($species."_".$table.".fasta",$output_dir);
	#unless ($translation_db) { $translation_db =  path($species."_".$table.".fasta",$output_dir)}
	
	write_output($transcript,$translation_db,$output_dir);
	
	timer($startRun);	# Get Run time
	print "\n";
}


##------ WRITE OUTPUT -------##
sub write_output {

	my $transcript = $_[0];
	my $output = $_[1];
	my $output_dir = $_[2];

	my $count_mapped = 0;		
	my %annotations_mapped = ();
	my %annotations = ();
	
	my $seq_out  = Bio::SeqIO->new(-file => ">$output", -format => "fasta");
	foreach my $tr (sort keys %$transcript) {

		my $id = "generic|".$tr."|".$transcript->{$tr}->{'gene'};
		my $desc = $transcript->{$tr}->{'codon'}." ".$transcript->{$tr}->{'aTIS_call'}." ".$transcript->{$tr}->{'peak_shift'};	
		
		if ($mflag == 1 or $mflag == 2) {
			if ($transcript->{$tr}->{'em'}) {$desc = $transcript->{$tr}->{'em'}." ".$desc}		# Add the percentage of the mapping to the description field
			$annotations_mapped{$transcript->{$tr}->{'anno'}}++;		# count mapping annotations

		} elsif ($mflag == 3) {
					
			if ($transcript->{$tr}->{'mp'}) {$desc = $desc." ".$transcript->{$tr}->{'mp'}."%"}		# Add the percentage of the mapping to the description field
			if ($transcript->{$tr}->{'or'}) {$desc = $transcript->{$tr}->{'or'}." ".$desc}
			if ($transcript->{$tr}->{'em'}) {
				$desc = $transcript->{$tr}->{'em'}." ".$desc;				# Add the percentage of the mapping to the description field
				$annotations_mapped{$transcript->{$tr}->{'anno'}}++;		# count mapping per annotations
				$count_mapped++;											# count mapped transcripts			
			}		

		}
		
		my $seq_obj = Bio::Seq->new(-display_id => $id, -desc => $desc, -seq => $transcript->{$tr}->{'seq'});
		$seq_out->write_seq($seq_obj);

		$annotations{$transcript->{$tr}->{'anno'}}++;
  	}

		# write Summary of external mapping to file
	my $percent_overall = sprintf '%.1f', 100*($count_mapped/$total_tr);
	my $nonRedundantTrans = scalar(keys %$transcript);
	my $percent_mapped = sprintf '%.1f', 100*($count_mapped/$nonRedundantTrans);
	my $percent = sprintf '%.1f', 100*($nonRedundantTrans/$total_tr);

	print "Totat transcripts in database:  ", $total_tr, "\n";
	print "Number of transcripts after cleaning out redundancy ",$nonRedundantTrans," ($percent%)\n";
	print "Number of transcripts mapped to external refernece: $count_mapped ($percent_mapped%)\n\n";
	
	foreach my $key (sort {$annotations{$b} <=> $annotations{$a} }keys %annotations) {
		if ($annotations_mapped{$key}) {
			print "Number of ", $key," annotations ", $annotations{$key}." ($annotations_mapped{$key})","\n";
		} else {
			print "Number of ", $key," annotations ", $annotations{$key}." (0)\n";
		}
	}

	print "Output files written to directory: $output_dir\n";
}

##------ SEQUENCE MAPPING -------##


	# Blast Search
sub blastp {

	my $transcript = $_[0];
	my $trans_to_blast = $_[1];
	
	print "Running Local blastp search...\n";
	my ($blast_rpt) = $trans_to_blast =~ /^(.*)\./;
	$blast_rpt = $blast_rpt.".bls";
			
	my $command = "blastp -query ".$trans_to_blast." -db ".$blastdb." -out ".$blast_rpt." -num_alignments 1 -num_descriptions 1 -evalue ".$evalue." -num_threads ".$nr_of_cores." -gapopen ".$gapopen." -gapextend ".$gapextend." -matrix ".$matrix." > /dev/null 2>&1";
	system($command);

	return blastp_parser($transcript, $blast_rpt);
}

sub ublast {

	my $transcript = $_[0];
	my $file_to_blast = $_[1];
	

	print "Running Local ublast search...\n";
	if ($identity > 1) {$identity = $identity/100;}

	my $file_size = (stat($file_to_blast))[7];
	$file_size = $file_size/1048576;	
	my $max_file_size = 1;
	
	if ($file_size <= $max_file_size) {		# if file size is <= 1mb don't split run the file directly
		my ($blast_rpt) = $file_to_blast =~ /^(.*)\./;
		$blast_rpt = $blast_rpt.".bls";
		system( "usearch -ublast ".$file_to_blast." -db ".$blastdb." -evalue ".$evalue." -id ".$identity." -alnout ".$blast_rpt." -maxhits 1 -lopen ".$gapopen." -lext ".$gapextend." -quiet");
		$transcript = ublast_parser($transcript, $blast_rpt);
		
	} else {	# if file size if > 1mb split and process files
		my $fasta_files = split_fasta_file($file_to_blast, $file_size, $max_file_size);
		
		foreach my $file (keys %$fasta_files) {		# loop thru all the files to perform ublast
		
			my ($blast_rpt) = $file =~ /^(.*)\./;
			$blast_rpt = $blast_rpt.".bls";
			system( "usearch -ublast ".$file." -db ".$blastdb." -evalue ".$evalue." -id ".$identity." -alnout ".$blast_rpt." -maxhits 1 -lopen ".$gapopen." -lext ".$gapextend." -quiet");
			$transcript = ublast_parser($transcript, $blast_rpt);
		}
	}

	return $transcript;
}

	# Blastp PARSER
sub blastp_parser {

	my ($transcript, $blast_file) = @_;
	my $in = new Bio::SearchIO( -format => 'blast', -file => $blast_file);

	while(my $result = $in->next_result) {
		while(my $hit = $result->next_hit) {
			while(my $hsp = $hit->next_hsp) {

				next if ($hsp->length('query') < $coverage);			# skip if length of identical positions < $coverage
				my $percent = sprintf '%.2f', $hsp->percent_identity;
				next if ($percent < $identity);
				my @en_id = split '\|', $result->query_name;

				my @target = split('\|',(substr($hit->name, 3)));

				$transcript->{$en_id[1]}->{'em'} = $target[0];

				$transcript->{$en_id[1]}->{'mp'} = $percent;		# Store sequence similarity percentage
				$transcript->{$en_id[1]}->{'or'} = $target[1];			# get other reference if available in canonical db
			}
		}
	}

	return $transcript;
}

	## ublast m8 PARSER
sub ublast_parser {

	# this function parsers the usearh ublast m8 output format

	my $transcript = $_[0];
	my $blast_rpt = $_[1];
		
	my %hash;					# Hash to store the parsed results results
	my $qry_id;					# store ID of the query sequences
	my $qry_seq;				# Query Sequence
	my $qry_range;				# Alignment range of the Query sequences
	my $tgt_id;					# ID of the matched/target sequence
	my $tgt_seq;				# Target Sequence from ublast
	my $tgt_range;				# lignment range of the Target sequences 
	my $length_of_alignment;	# length of the alignment 
	my $identical;				# Percentage of Identical match

	open (IN, $blast_rpt) or die "error reading file: ", $blast_rpt,"\n";
	while (<IN>) {

		next if (/usearch/);				
		next if (/^\s*$/);			

		if (/Query\s>/) {	# check if we are at the beginning of a new record.
					# if at a new record populate the hash with current values.
			if ($tgt_id) {
			
				my @target = split('\|',$tgt_id);
				$transcript->{$qry_id}->{'em'} = $target[1];		# mapped canonical transccript id
				$transcript->{$qry_id}->{'or'} = $target[2];		# add gene info to other refreence field
			}

			if ($identical) {$transcript->{$qry_id}->{'mp'} = $identical;}

			$qry_seq 		= "";			# clear out old values
			$tgt_seq 		= "";
			$qry_id			= "";
			$tgt_id 		= "";
			$tgt_range 		= "";
			$qry_range 		= "";
			$identical		= "";
			$length_of_alignment 	= "";

			my @id = split '\|', $_;
			$qry_id = $id[1];			# get the query id
			
			if ($qry_id) {$qry_id =~ s/\s+$//;}		# remove all white spaces at the end

		} elsif (/^\s+\d+.*/) {
			my @value 	= split(/\s+/, $_);
			$tgt_id 	= $value[scalar(@value) - 1]; 	
			($tgt_range 	= $value[scalar(@value) - 2]) =~ s/(\(\d+\))|\s//;	# Get target range
			($qry_range 	= $value[scalar(@value) - 3]) =~ s/(\(\d+\))|\s//;	# Get query range
	
		} elsif (/Qry\s/) {
			s/^Qry|\s+|\d+//g;
			$qry_seq .= $_;		# combine the query sequences (if in multiple lines)

	   	} elsif (/Tgt\s/) {
			s/^Tgt|\s+|\d+//g;
			$tgt_seq .= $_;		# combine the target sequences (if in multiple lines)

		} elsif (/^\d+\s+[cols].*/) {
		
			($length_of_alignment) = $_ =~ /\s*(\d+)\s+cols/;
			my @ident = split ',', $_;
			my @id = split '\(', $ident[1]; 	# get the blast percentage
			$id[1] =~ s/\)|\%//g;
			
			$identical = sprintf '%.2f', $id[1];
		}
	 }    
	 close IN;

	if ($qry_id) { 		# handle last result
		unless ($transcript->{$qry_id}) {
			my @t_range = split('-',$tgt_range);
			my @q_range = split('-',$qry_range);

			my @target = split('\|',(substr($tgt_id, 3)));
			$transcript->{$qry_id}->{'em'} = $target[0];
			$transcript->{$qry_id}->{'mp'} = $identical;
		}
	}
	 
	return $transcript;
}

sub transcripts_to_blast {

	my $transcript = $_[0];
	my $tis = $_[1];
	
	my $fname = "TB_".$species."_".$tis.".fa";
	my $blast_file = path($fname, $TMP); # genearte file to store output
	my $out = new Bio::SeqIO(-file => ">".$blast_file, -format=>'fasta');

	foreach my $tr (keys %$transcript) {
		next if (length($transcript->{$tr}->{'seq'}) < $min_blast_length);	# skip all transripts with sequences < $min_blast_length
		next unless ($top_anno{$transcript->{$tr}->{'anno'}});	# skip all transripts which are not aTIS or 5UTR
		unless ($transcript->{$tr}->{'em'}) {
			my $seq_obj = Bio::Seq-> new (-display_id => "generic\|".$tr, -seq => $transcript->{$tr}->{'seq'});
			$out->write_seq($seq_obj);
			$no_blast_tr++;
		}
	}
	
	return $blast_file;
}

sub split_fasta_file {

	# this function splits the transcript files for use with usearch (32bit version)

	my ($in_file,$tmp_fld, $f_size, $max_f_size) = @_;

	my $num_files = int(($f_size/$max_f_size) + 0.5);
	print "number of files:$num_files \n";
	my $seqs_per_file = int(($no_blast_tr/$num_files) + 0.5);	# number of sequences in each file

	my $fcount = 1;			# Count number of files
	my $count_seq = 1;		# Count number of sequences in each file

	my $in  = new Bio::SeqIO(-file  => $in_file, -format=>'fasta');	
	my $out_file = path("split_",$TMP);
	my $out = new Bio::SeqIO(-file => ">".$out_file.$fcount.".fasta", -format=>'fasta');
	
	my $files = {};
	while (my $seq = $in->next_seq) {
	
		if ($count_seq % $seqs_per_file == 0) {		# if number of sequences = required number of sequences 
		    $fcount++;
			my $file = $out_file.$fcount.".fasta";
		    $out = new Bio::SeqIO(-file => ">".$file, -format=>'fasta');
			$files->{$file} = 1;
			$count_seq = 1;
		}
		
		$out->write_seq($seq);
		$count_seq++;	
	}

	return $files;
}

##------ ID MAPPING -------##
sub id_based_mapping {

	my $transcript = $_[0];
	my $mapping = $_[1];
	
	my $count = 0;
	foreach my $tr (keys %$transcript) {
		if ($mapping->{$transcript->{$tr}->{'tr'}}) {
			#if ($top_anno{$transcript->{$tr}->{'anno'}}) {
			
			$transcript->{$tr}->{'em'} = $mapping->{$transcript->{$tr}->{'tr'}}->{'em'};
			$transcript->{$tr}->{'or'} = $mapping->{$transcript->{$tr}->{'tr'}}->{'or'};
			$count++;
			#}
		}
	}
	
	print "Total number of transcripts mapped to a canonical ID.\n";
	return $transcript;
}

##------ REMOVE REDUNDANCY -------##
sub remove_redundancy {

	my $transcript = $_[0];
	my $gene_transcript = $_[1];
	
	print "Removing redundant sequences, this might take a while. Please wait...\n";
	foreach my $tr1 (keys %$transcript) {
		
		next if ($transcript->{$tr1}->{'red'} eq "Y");	# skip if transcript is already a subset of another
		my $gene = $transcript->{$tr1}->{'gene'};
		my $seq1 = $transcript->{$tr1}->{'seq'};
		my $anno_rank1 = $annotation{$transcript->{$tr1}->{'anno'}};
		
		foreach my $tr2 (keys %$transcript) {
		#foreach my $tr2 (@{$gene_transcript->{$gene}}) {
			next if ($transcript->{$tr2}->{'red'} eq "Y");	# skip if transcript is already a subset of another
			my $seq2 = $transcript->{$tr2}->{'seq'};
			my $anno_rank2 = $annotation{$transcript->{$tr2}->{'anno'}};

			if ($seq1 eq $seq2) {
				# check annotation
				if ($anno_rank1 < $anno_rank2) {		# if tr1 annotation is ranked higer that tr2
					$transcript->{$tr2}->{'red'} = "Y";
				} elsif ($anno_rank1 > $anno_rank2) {
					$transcript->{$tr1}->{'red'} = "Y";
				}
			} elsif (index($seq1,$seq2) > 0) {
				$transcript->{$tr2}->{'red'} = "Y";
			} elsif (index($seq2,$seq1) > 0) {
				$transcript->{$tr1}->{'red'} = "Y";
			}
		}
	}

	my $non_red_trans = {};
	foreach my $tr (keys $transcript) {
		next if ($transcript->{$tr}->{'red'} eq "Y");
		$non_red_trans->{$tr} = $transcript->{$tr};
	}
	
	return $non_red_trans;
}


##------ GET TRANSCRIPTS -------##
sub get_transcripts_from_resultdb {

	my ($dbh,$tbl,$tis) = @_;

	my $transcript = {}; 			# collect sequences from SQLite DB
	my $gene_transcript = {};
	
	print "Extracting transcripts form SQLite database. Please wait ....\n";
	my $query = "SELECT DISTINCT a.tr_stable_id, a.chr, a.start, a.start_codon, a.dist_to_aTIS, a.aTIS_call, a.annotation, a.peak_shift, a.SNP, a.aa_seq, b.biotype, c.gene_stable_id
				 FROM ".$tbl." a JOIN TIS_".$tis." b JOIN tr_translation c
				 WHERE a.tr_stable_id = b.stable_id AND b.transcript_id = c.transcript_id";

 	my $sth = $dbh->prepare($query);
	$sth->execute();
	
	while ( my ($tr_stable_id, $chr, $start, $start_codon, $dist_to_aTIS, $aTIS_call, $annotation, $peak_shift, $snp, $aa_seq, $biotype,$gene_stable_id) = $sth->fetchrow()) {
		
		if (uc($tis_call) eq "N") {
			next if (uc($aTIS_call) eq 'NO_DATA' or uc($aTIS_call) eq 'FALSE');
		}
			# to be removed
		if ($annotation eq "no_translation") {$annotation = "ntr"}
		
		$aa_seq =~ s/\*//g;
		next if (length($aa_seq) < $mslength);		# skip if sequence is less than minimum allowed amino acid length
		my $tr = $tr_stable_id."_".$chr."_".$start."_".$annotation;

		$transcript->{$tr}->{'chr'} = $chr;
		$transcript->{$tr}->{'codon'} = $start_codon;
		$transcript->{$tr}->{'aTIS_call'} = $aTIS_call;
		$transcript->{$tr}->{'anno'} = $annotation;
		$transcript->{$tr}->{'peak_shift'} = $peak_shift;
		$transcript->{$tr}->{'seq'} = $aa_seq;
		$transcript->{$tr}->{'biotype'} = $biotype;
		$transcript->{$tr}->{'gene'} = $gene_stable_id;
		$transcript->{$tr}->{'tr'} = $tr_stable_id;
		$transcript->{$tr}->{'red'} = "N";			# Set all transcript to default redundant value of N
		
		push @{$gene_transcript->{$gene_stable_id}}, $tr;
	
	}
	$sth->finish();
	
	print "Done.\n";
	return $transcript,$gene_transcript;
}


##------ MAPPING TO CANONICAL -------##
sub remote_biomart_mapping {
	
	# This function remotely searches Biomart website for ensemble transcript to swissprot ID mapppings
	# input is the ensemble database name
	
	my $mapping_db = $_[0];
	my $version = $_[1];
	my $ref_attribute = $_[2];
	my $other_attribute = $_[3];
	
	print "Retriving the list of Biomart Ensembl transcript to Swissprot ids mappings.\n";
	print "Please wait.....\n";
	
		# Construct the XML query
	my $xml = XML::Smart->new();
	$xml->{Query}{virtualSchemaName} = "default";
	$xml->{Query}{formatter} = "CSV";
	$xml->{Query}{header} = "0";
	$xml->{Query}{uniqueRows} = "0";
	$xml->{Query}{count} = "";
	$xml->{Query}{datasetConfigVersion} = $version;
	$xml->{Query}{Dataset}{name} = $mapping_db;
	$xml->{Query}{Dataset}{interface} = "default";
	$xml->{Query}{Dataset}{Attribute}[0]{name} = "ensembl_transcript_id";
	$xml->{Query}{Dataset}{Attribute}[1]{name} = $ref_attribute;
	if ($extra_info) {
		$xml->{Query}{Dataset}{Attribute}[2]{name} = $other_attribute;
	}
	my $xmlQuery = $xml->data;
	$xmlQuery =~ s/\n//g;
	$xmlQuery =~ s/>\s+</></g;

	my $path="http://www.biomart.org/biomart/martservice?";
	my $request = HTTP::Request->new("POST",$path,HTTP::Headers->new(),'query='.$xmlQuery."\n"); # post the query to the biomart website
	my $ua = LWP::UserAgent->new;
	my $response;
	my %trans;

	my $mpfile = path($species."_biomart_mapped.txt", $work_dir);
	open FH, ">".$mpfile or die file_err($mpfile, 1);
	$ua->request($request, 
		sub {   
			my($data, $response) = @_;
			if ($response->is_success) {
				print FH $_[0];
			}
			else {
				warn ("Problems with the web server: ".$response->status_line);
			}
		});
	close(FH);
	
	return $mpfile;
}

sub biomart_mapping {

	my $map_file = shift;
	
	my $mapping_hash = ();

	open (MYIDS, $map_file) or die file_err($map_file,0);
	while (<MYIDS>) {

		chomp;
		my @line = split ',', $_ ;
		my $tran_id = $line[0];
		my $em = $line[1];
		my $other = $line[2]; 

		if ($em) {	
			$mapping_hash->{$tran_id}->{'em'} = $em;
			if (defined $other) {$mapping_hash->{$tran_id}->{'or'} = $other;}
		}
	}
	close (MYIDS);

	#print "Number of transcripts mapped to a Canonical Protein sequence ".scalar(keys %$mapping_hash)."\n";
	return $mapping_hash;
	
}


##------ DBH -------##
sub dbh {
    # Catch

    my $db  = $_[0];
    my $us	= $_[1];
    my $pw	= $_[2];

    # Init DB
    my $dbh = DBI->connect($db,$us,$pw,{ RaiseError => 1 },) || die "Cannot connect: " . $DBI::errstr;
    return($dbh);
}

sub get_analysis_ids {
    
	# Catch
	my $dbh    = $_[0];
	my $ids_in = $_[1]; #Either comma separated list of identifiers or "all"
	my @idsref;

	if ($ids_in eq "all") {

		my $query = "select ID, SNP from TIS_overview";
	    my $sth = $dbh->prepare($query);
		$sth->execute();

		while ( my ($id, $snp) = $sth->fetchrow()) {
			if (uc($snp) eq "NO") {
				push @idsref, $id;
			} else {
				push @idsref, $id."_".$snp;
			}
		}
		$sth->finish();
    } else {
        	
		my @sel_ids  = split(',',$ids_in);
		foreach (@sel_ids) {

			my $query = "select ID, SNP from TIS_overview where ID = $_";
		    my $sth = $dbh->prepare($query);
			$sth->execute();
			my ($id, $snp) = $sth->fetchrow();
			
			if (uc($snp) eq "NO") {
				push @idsref, $id;
			} else {
				push @idsref, $id."_".$snp;
			}
			$sth->finish();
		}
    }
	
    return @idsref;
}


##------ GET_INPUT_VARS -------##
sub get_input_vars {

    # Catch

    my $dbh_results = $_[0];

    my ($query,$sth);

    # Get input variables

    $query = "select value from arguments where variable = \'run_name\'";
    $sth = $dbh_results->prepare($query);
	$sth->execute();
	my $run_name = $sth->fetch()->[0];
	
    $query = "select value from arguments where variable = \'species\'";
    $sth = $dbh_results->prepare($query);
	$sth->execute();
	my $species = $sth->fetch()->[0];

    $query = "select value from arguments where variable = \'mapper\'";
    $sth = $dbh_results->prepare($query);
	$sth->execute();
	my $mapper = $sth->fetch()->[0];

    $query = "select value from arguments where variable = \'nr_of_cores\'";
    $sth = $dbh_results->prepare($query);
	$sth->execute();
	my $nr_of_cores = $sth->fetch()->[0];

    # Return input variables

    return($run_name,$species,$mapper,$nr_of_cores);

}

### CHECK REQUIRED PARAMETERS ###
sub uninitialized_param {
	
	my ($v) = @_;
	not ( defined($v) and length $v );
}

### HANDLE PATH TO FILES ###
sub path {

	my ($file, $dir) = @_;
	return (File::Spec->catfile( $dir, $file));
}

### FILE ERROR HANDLING ###
sub file_err {

	my ($file, $flag) = @_;		# flag determines if the operation is open (1) or read (0)
	if ($flag == 1) {
		print "Error creating file ", $file, ". Ensure you have the required permission.\n";
	} else {
		print "Error reading file ", $file, ". Ensure the file exist.\n";
	}
	
}

sub timer {
	my $startRun = shift;

	my $endRun 	= time();
	my $runTime = $endRun - $startRun;

	printf("\nTotal running time: %02d:%02d:%02d\n\n", int($runTime / 3600), int(($runTime  % 3600) / 60), int($runTime % 60));
}
