#!/usr/bin/perl -w

use strict;
use warnings;
use POSIX ":sys_wait_h";
use Getopt::Long;
use Cwd;
use Bio::SearchIO;
use Bio::SeqIO;
use DBI;
use LWP::UserAgent; 
use XML::Smart;
use Parallel::ForkManager;

# ---------------------------------------------------------------------
	##	GLOBAL VARIABLES

my $countTranscripts = 0;
my $num_trans_to_Blast = 0;
my $num_non_Red_trans = 0;

# ---------------------------------------------------------------------
	##	REQUIRED PARAMETERS

my $blastdb;		# Usearch/blastp formatted database
my $sqlitedb;		# SQLite source database
my $tis_ids;		# TIS ID to generate table from
my $mapping_db;		# Ensembl database name to download biomart mapped transcripts 
my $mflag;		# Flag for biomart mappings, 0 = remote download, 1 = local file, 2 = no biomart mapping.

# ---------------------------------------------------------------------
	##	DEFAULT prameters

my $password		= '';		# DB password
my $user		= '';		# user name to access SQL DB
my $pgm 		= 'blastp';	# The program to use for search, Default program = blastp
my $minBlastLength 	= 32;		# minimum sequence length to perform blast
my $coverage 		= 30;		# Minimum number of identical positions
my $minSeqLength 	= 6;		# Minimum transcript length to allow in the translation product database
my $evalue 		= 1e-10;	# blast e-value
my $identity 		= 75;		# Minimum alignment score
my $gapopen		= 11;		# Gap opening penalty
my $gapextend		= 1;		# Gap extension penalty
my $matrix		= 'BLOSUM62';	# Blast search matrix
my $word_size		= 3;		# word size
my $numthreads		= 3;		# Number of processes to use
my $summary;				# Summary of the mapping to external reference and annotations
my $translation_db;			# FASTA file of non redundant derived translation products
#my $blast_report;			# blast report file
my $work_dir;				# Working directory

my $external_ref	= "uniprot_swissprot_accession";
my $extra_info 		= "uniprot_swissprot";

# ---------------------------------------------------------------------
	##	GET command line arguments
GetOptions(
	'password=s'		=> \$password,
	'user=s'		=> \$user,
	'blast_db=s' 	 	=> \$blastdb,
	'sqlite_db=s' 	 	=> \$sqlitedb,
	'tis_ids=s' 		=> \$tis_ids,
	'work_dir=s' 	 	=> \$work_dir,
	'blast_pgm=s' 	 	=> \$pgm,
	'evalue=f' 	 	=> \$evalue,
	'min_blast_length=i'	=> \$minBlastLength,
	#'num_threads=i'	 	=> \$numthreads,
	'mslenght=i'	 	=> \$minSeqLength,
	'mflag=i'	 	=> \$mflag,
	'identity=f'		=> \$identity,
	'coverage=i'		=> \$coverage,
	'mapping_db=s'		=> \$mapping_db,
	'external_ref=s'	=> \$external_ref,
	'extra_info=s'		=> \$extra_info,
	'gapopen=f'	 	=> \$gapopen,
	'gapextend=f'	 	=> \$gapextend,
	'word_size=i'	 	=> \$word_size,
	'matrix=s'	 	=> \$matrix,
	'summary=s'		=> \$summary,
	'translation_db=s'	=> \$translation_db,
);

# ---------------------------------------------------------------------
	## EXECUTION

## Command line
# perl generate_db.pl -blast_db /home/ndah/blastdb/Mus_Musculus/Mus_Musculus.udb -sqlite_db results_STARclip_26-35_1-14MM.db -tis_ids 1 -blast_pgm ublast -mflag 0 -external_ref uniprot_swissprot_accession -mapping_db mmusculus_gene_ensembl -num_threads 3 -work_dir work_directory
# mmusculus_gene_ensembl
# mart_export_mm.txt

# perl generate_db.pl -blast_db /home/ndah/blastdb/Mus_Musculus/Mus_Musculus.udb -sqlite_db results_STARclip_26-35_1-14MM.db -tis_ids 1 -blast_pgm ublast -mflag 1 -mapping_db mart_export_mm.txt -num_threads 3 -work_dir work_directory
# mmusculus_gene_ensembl
# mart_export_mm.txt
#####

###################
# For scripts
#
###################
#my %params = (
#	blast_db 	=> $blastdb,
#	sqlite_db 	=> $sqlitedb,
#	tis_ids		=> $tis_ids,
#	blast_pgm 	=> $pgm,
#	mflag		=> $mflag	
#);

	#check if para,eters are properly initialized
#my @invalid = grep uninitialized_param($params{$_}), keys %params;	
#die "Not properly initialized: @invalid\n" if @invalid;

my $cwd = getcwd();
if (!($work_dir)) {
	$work_dir = $cwd;
} elsif (!-d "$work_dir") { 
	system ("mkdir ". $work_dir);
} 

print "\n";
if ($work_dir) {
	print "Working directory .............................................. $work_dir\n";
} 
if ($sqlitedb) {
	print "SQLite database containing transcripts ......................... $sqlitedb\n";
} 
if ($tis_ids) {
	print "TIS id to generate translation database......................... $tis_ids\n";
}
if ($pgm) {
	print "Blast program .................................................. $pgm\n"; 
} 
if ($mflag) {

	if ($mflag == 0) {
		print "Ensemble database for remote transcript to Swissprot mapping ... $mapping_db\n";
		print "External reference to map transcript with ...................... $external_ref\n";
	} elsif ($mflag == 1) {
		print "Locale file containing ensembl to Swissprot mapping ............ $mapping_db\n";
	} else {
		print "All transcripts will be mapped by $pgm search\n"
	}

} else {
	print "Resulting derived translation product database will not be mapped to any known canonical database\n";
}
if ($blastdb) {
	print "The blast database is .......................................... $blastdb\n";
} 
if ($evalue) {
	print "Blast e-value .................................................. $evalue\n";
} 
if ($minBlastLength) {
	print "Minimum sequence length allowed for Blast search ............... $minBlastLength\n";
} 
if ($identity) {
	print "Blast Identity value ........................................... $identity%\n";
} 
if ($coverage) {
	print "Minimum percentage of identical positions ...................... $coverage\n";
}
print "\n";

my %annotationHash = ( "aTIS",1,"5UTR",2,"CDS",3,"ntr",4,"3UTR",5);
my @TopAnno = ("aTIS","5UTR");

print "Annotations and their rankings (with 1 the most important) \n";
foreach my $key (sort {$annotationHash{$a} <=> $annotationHash{$b}} keys %annotationHash) {
    print "$key\t$annotationHash{$key}\n";
}

my $external_mapping = {};
if ($mflag) {

	if ($mflag == 0) {
		my $mapping_file = remote_biomart_mapping();
		$external_mapping = biomart_mapping($mapping_file);

	} elsif ($mflag == 1)  {
		$external_mapping = biomart_mapping($mapping_db);

	} 
}

# GET TABLES TO GENERATE DB
my $dsn_results = "DBI:SQLite:dbname=$sqlitedb";
my $dbh_results = dbh($dsn_results,$user,$password);
my @idsref = get_analysis_ids($dbh_results,$tis_ids);  #$tis_ids is input variable

$numthreads = get_cores($dbh_results);
if ($numthreads) {
	print "Number of cores used ........................................... $numthreads\n";
}

foreach (@idsref) {	# generate translation db for selected tis_ids 
	generate_trans_db($_);
}


# ---------------------------------------------------------------------
	## SUBROUTINES

sub generate_trans_db {

	my $tis_id = shift;

	my $startRun = time();
	my $table = "tis_".$tis_id."_transcripts";

	print "Generating translation product database for table ", $table, "\n";

	my $tmp_folder = $work_dir."/temp_".$tis_id;	# Create temporary folder
	if (!-d "$tmp_folder") { system ("mkdir ". $tmp_folder)} 

	unless ($translation_db) {$translation_db = path(substr($table, rindex($table,".")+1).".fasta", $work_dir);}
	unless ($summary) {$summary = path($tis_id."_summary.txt", $work_dir);}

	my $blast_report = path($tis_id."_blast_report.bls", $tmp_folder);
	my $transcriptHOH = transcripts_from_sqldb($table, $external_mapping);
	my $transcriptsToBlast = transcripts_to_blast($transcriptHOH, $tmp_folder);

	my $final_transcriptsHOH;
	if ($pgm eq 'ublast') {

		$final_transcriptsHOH = ublast($transcriptHOH, $transcriptsToBlast, $tmp_folder, $blast_report);
		
	} elsif ($pgm eq 'blastp') {

		$final_transcriptsHOH = blastp($transcriptHOH, $transcriptsToBlast, $blast_report);
		
	} elsif ($pgm eq 'none') {

		$final_transcriptsHOH = $transcriptHOH;
	} 
	
	write_output($final_transcriptsHOH, $translation_db, $summary);

	timer($startRun);	# Get Run time
	print "\n";
}


####################################################################
	#
	# GET tis_ids FOR PROCESSING
	#
####################################################################

sub get_cores {

	my $dbh = shift;

	my $query = "select value from arguments where variable = \'nr_of_cores\'";
	my $sth = $dbh_results->prepare($query);
	$sth->execute();
	my $cores = $sth->fetch()->[0];
	$sth->finish();

	return $cores;
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


sub dbh {
    
    # Catch
    my $db  = $_[0];
    my $us	= $_[1];
    my $pw	= $_[2];
    
    # Init DB
    my $dbh = DBI->connect($db,$us,$pw,{ RaiseError => 1 },) || die "Cannot connect: " . $DBI::errstr;
    
    return($dbh);
}

####################################################################
	#
	# GET BIOMART MAPPING
	#
####################################################################

	  
sub remote_biomart_mapping {
	
	# This function remotely searches Biomart website for ensemble transcript to swissprot ID mapppings
	# input is the ensemble database name
	
	print "Retriving the list of Biomart Ensembl transcript to Swissprot ids mappings.\n";
	
		# Construct the XML query
	my $xml = XML::Smart->new();
	$xml->{Query}{virtualSchemaName} = "default";
	$xml->{Query}{formatter} = "CSV";
	$xml->{Query}{header} = "0";
	$xml->{Query}{uniqueRows} = "0";
	$xml->{Query}{count} = "";
	$xml->{Query}{datasetConfigVersion} = "0.6";
	$xml->{Query}{Dataset}{name} = $mapping_db;
	$xml->{Query}{Dataset}{interface} = "default";
	$xml->{Query}{Dataset}{Attribute}[0]{name} = "ensembl_transcript_id";
	$xml->{Query}{Dataset}{Attribute}[1]{name} = $external_ref;
	if ($extra_info) {
		$xml->{Query}{Dataset}{Attribute}[2]{name} = $extra_info;
	}
	my $xmlQuery = $xml->data;
	$xmlQuery =~ s/\n//g;
	$xmlQuery =~ s/>\s+</></g;

	my $path="http://www.biomart.org/biomart/martservice?";
	my $request = HTTP::Request->new("POST",$path,HTTP::Headers->new(),'query='.$xmlQuery."\n"); # post the query to the biomart website
	my $ua = LWP::UserAgent->new;
	my $response;
	my %trans;

	my $mpfile = path("biomart_mapped.txt", $work_dir);
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
	my %mapping_hash = ();

	if ($mflag == 0 or $mflag == 1) { 	# No pre-mapping from file
	
		open (MYIDS, $map_file) or die file_err($map_file,0);
		while (<MYIDS>) {

			chomp;
			my @line = split ',', $_ ;
			my $tran_id = shift @line;
			my $ext_ref = shift @line;
			my $other = join(' ', @line); 

			if (defined $tran_id && defined $tran_id) {
			
				$mapping_hash{$tran_id}{'ext_ref'} = $ext_ref;
				$mapping_hash{$tran_id}{'other_ref'} = $other;
			}
		}
		
		close (MYIDS);
		
	} else {
		print "No external pre-mapping perform.\n";
	}

	return(\%mapping_hash);
}


####################################################################
	#
	# PERFORM BLAST SEARCH
	#
####################################################################

	# Blast Search
sub blastp {

	my ($transHoH, $transFile, $blast_rpt) 	= @_;
	
	print "Running Local blastp search...\n";
	
	my $command = "blastp -query ".$transFile." -db ".$blastdb." -out ".$blast_rpt." -num_alignments 1 -num_descriptions 1 -evalue ".$evalue." -num_threads ".$numthreads." -gapopen ".$gapopen." -gapextend ".$gapextend." -matrix ".$matrix." > /dev/null 2>&1";
	system($command);

	return blastp_parser($transHoH, $blast_rpt);
}


sub ublast {

	my ($trans_HoH, $trans_file, $split_fdr, $blast_report) = @_;

	print "Running Local ublast search...\n";
	if ($identity > 1) {$identity = $identity/100;}

	my $fasta_dir = split_fasta_file($trans_file, $split_fdr);

	opendir(DIRECTORY, $fasta_dir) or die $!;
	while (my $file = readdir(DIRECTORY)) { 

		if($file=~/^\./){next;}
		my $file_in = path($file, $fasta_dir);

		system( "usearch -ublast ".$file_in." -db ".$blastdb." -evalue ".$evalue." -id ".$identity." -alnout ".$blast_report." -maxhits 1 -lopen ".$gapopen." -lext ".$gapextend." -quiet");
		
		$trans_HoH = ublast_parser($trans_HoH, $blast_report);

	}
	close DIRECTORY;

	return $trans_HoH;
}


####################################################################
	#
	# PARSE BLAST REPORTS
	#
####################################################################

	# Blastp PARSER
sub blastp_parser {

	my ($transcriptsHoH, $blast_file) = @_;
	my $in = new Bio::SearchIO( -format => 'blast', -file => $blast_file);

	while(my $result = $in->next_result) {
		while(my $hit = $result->next_hit) {
			while(my $hsp = $hit->next_hsp) {

				next if ($hsp->percent_identity < $identity);			# skip if hsp_identity < $identity
				next if ($hsp->length('query') < $coverage);			# skip if length of identical positions < $coverage
				my @en_id = split '\|', $result->query_name;

				my @target = split('\|',(substr($hit->name, 3)));

				$$transcriptsHoH{$en_id[1]}{'ext_ref'} = $target[0];
				$$transcriptsHoH{$en_id[1]}{'des'} = $target[1]." ".$$transcriptsHoH{$en_id[1]}{'des'};
				my $percent = sprintf '%.2f', $hsp->percent_identity;
				$$transcriptsHoH{$en_id[1]}{'des'}	= $$transcriptsHoH{$en_id[1]}{'des'}." ".$percent."%";
			}
		}
	}

	return $transcriptsHoH;
}


	## ublast m8 PARSER
sub ublast_parser {

	# this function parsers the usearh ublast m8 output format

	my ($transcriptsHoH, $fasta_file) = @_;
	my %hash;				# Hash to store the parsed results results
	my $qry_id;				# store ID of the query sequences
	my $qry_seq;				# Query Sequence
	my $qry_range;				# Alignment range of the Query sequences
	my $tgt_id;				# ID of the matched/target sequence
	my $tgt_seq;				# Target Sequence from ublast
	my $tgt_range;				# lignment range of the Target sequences 
	my $length_of_alignment;		# length of the alignment 
	my $identical;				# Percentage of Identical match

	open (IN, $fasta_file) or die "error reading file: ", $fasta_file,"\n";
	while (<IN>) {

		next if (/usearch/);				
		next if (/^\s*$/);			

		if (/Query\s>/) {	# check if we are at the beginning of a new record.
					# if at a new record populate the hash with current values.
			if ($tgt_id) {
			
				my @target = split('\|',(substr($tgt_id, 3)));
				$$transcriptsHoH{$qry_id}{'ext_ref'} = $target[0];
				$$transcriptsHoH{$qry_id}{'des'} = $target[1]." ".$$transcriptsHoH{$qry_id}{'des'};
			}

			if ($identical) {$$transcriptsHoH{$qry_id}{'des'} = $$transcriptsHoH{$qry_id}{'des'}." ".$identical."%";}

			$qry_seq 		= "";			# clear out old values
			$tgt_seq 		= "";
			$qry_id			= "";
			$tgt_id 		= "";
			$tgt_range 		= "";
			$qry_range 		= "";
			$identical		= "";
			$length_of_alignment 	= "";

			my @id = split '\|', $_;
			($qry_id) = $id[1];
			if ($qry_id) {$qry_id =~ s/\s+$//;}

		} elsif (/^\s+\d+.*/) {
			my @value 	= split(/\s+/, $_);
			$tgt_id 	= $value[scalar(@value) - 1]; 				# Target swissprot ID
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
			my @id = split '\(', $ident[1]; $id[1] =~ s/\)|\%//g;
			$identical = sprintf '%.2f', $id[1];
			
		}
	 }    
	 close IN;

	if ($qry_id) { 		# handle last result
		my @t_range = split('-',$tgt_range);
		my @q_range = split('-',$qry_range);

		my @target = split('\|',(substr($tgt_id, 3)));
		$$transcriptsHoH{$qry_id}{'ext_ref'} = $target[0];
		$$transcriptsHoH{$qry_id}{'des'} = $target[1]." ".$$transcriptsHoH{$qry_id}{'des'};
		$$transcriptsHoH{$qry_id}{'des'} = $$transcriptsHoH{$qry_id}{'des'}." ".$identical."%";
		
	}
	 
	return $transcriptsHoH;
}


####################################################################
	#
	# PREPARE BLAST INPUT FILES
	#
####################################################################

sub transcripts_to_blast {

	# This subroutine write transcripts with no biomart external mapping to file

	my ($transcriptHOH, $temp_fdr) = @_;
	my $trans_2_blast = path("trans2blast.fa", $temp_fdr);

	my $out = new Bio::SeqIO(-file => ">".$trans_2_blast, -format=>'fasta');

	foreach my $key (keys %$transcriptHOH) {
		next if (length($$transcriptHOH{$key}{'seq'}) < $minBlastLength);

		unless ($$transcriptHOH{$key}{'ext_ref'}) {
			my $seq_obj = Bio::Seq-> new (-display_id => "generic\|".$key, -seq => $$transcriptHOH{$key}{'seq'});
			$out->write_seq($seq_obj);
			$num_trans_to_Blast++;
		}
	}
	
	return $trans_2_blast;
}

sub split_fasta_file {

	# this function splits the transcript files foor use with usearch (32bit version)

	my ($inFile,$tmp_fld) = @_;
	my $numFiles = 4;
	my $seqsPerFile = int(($num_trans_to_Blast/$numFiles) + 0.5);	# number of sequences in each file
	my $fcount = 1;								# Count number of files
	my $countSeq = 1;							# Exit count

	my $fasta_dir = $tmp_fld."/trans2blast/";
	if (!-d "$fasta_dir") { system ("mkdir ". $fasta_dir) } 

	my $in  = new Bio::SeqIO(-file  => $inFile, -format=>'fasta');
	
	my $outFile = path("split_",$fasta_dir);
	my $out = new Bio::SeqIO(-file => ">".$outFile.$fcount."fasta", -format=>'fasta');
	
	while (my $seq = $in->next_seq) {
	
		if ($countSeq % $seqsPerFile == 0) {
		
		    	$fcount++;
		    	$out = new Bio::SeqIO(-file => ">".$outFile.$fcount."fasta", -format=>'fasta');
			$countSeq = 1;
		}
		
		$out->write_seq($seq);
		$countSeq++;	
	}

	return $fasta_dir;
}


####################################################################
	#
	# WRITE OUTPUT FILES
	#
####################################################################

	
	# get mutation information
sub get_variations {

	my ($qry_seq, $qry_start, $qry_end, $tgt_seq, $tgt_start, $tgt_end, $length_of_alignment) = @_;

	my $qry_length 	= $qry_end - $qry_start + 1;		# get the length of the query sequence				
	if ($qry_length < $length_of_alignment) {$qry_end = $qry_end + ($length_of_alignment - $qry_length)}
	my @query 	= ($qry_start...$qry_end);

	my $tgt_length 	= $tgt_end - $tgt_start + 1;		# get length of targe sequence			
	if ($tgt_length < $length_of_alignment) {$tgt_end = $tgt_end + ($length_of_alignment - $tgt_length)}
	my @tgts 	= ($tgt_start...$tgt_end);

	my $diff 	= $qry_seq ^ $tgt_seq;			# Get positions of mismatched
	$diff 		=~ tr{\x00-\xff}{.*};

	my @unmatched;						# Array to store mismatched positions
	for (my $i = 0; $i < length($diff); $i++) {
		if (substr($diff,$i,1) eq "*") {push(@unmatched, $i)} 
	}

	my $variation = '';
	if (scalar(@unmatched) > 0) {
		for  (my $i = 0 ; $i < @unmatched; $i++) {
			my $tgt_char 	= substr($tgt_seq,$unmatched[$i],1);
			my $qry_char 	= substr($qry_seq,$unmatched[$i],1);
			my $tgt_pos 	= $tgts[$unmatched[$i]];
			my $qry_pos 	= $query[$unmatched[$i]];
			$variation .= $tgt_pos."\\".$qry_pos."_".$tgt_char."\\".$qry_char.":";			
		}
	}

	return $variation;
}

sub write_output {

	# Write output FASTA files
	my ($transcript, $output, $summary) = @_;

	my $count_mapped = 0;		
	my %annotations_mapped = ();
	my %annotations = ();

	
	my $seq_out  = Bio::SeqIO->new(-file => ">$output", -format => "fasta");

	foreach my $enID (sort keys %$transcript) {

		my $id = "generic|".$enID;

		if (defined $$transcript{$enID}{'ext_ref'}) {

			$id =$id."|".$$transcript{$enID}{'ext_ref'};
			$annotations_mapped{$$transcript{$enID}{'anno'}}++;

			$count_mapped++;

		} else { $id =$id."|" }

		my $seq_obj = Bio::Seq->new(-display_id => $id, -desc => $$transcript{$enID}{'des'}, -seq => $$transcript{$enID }{'seq'});
		$seq_out->write_seq($seq_obj);

		$annotations{$$transcript{$enID}{'anno'}}++;
  	}

		# write Summary of external mapping to file
	open SUM, ">$summary" or die file_err($summary,1);
	my $percent_overall = sprintf '%.2f', 100*($count_mapped/$countTranscripts);
	
	my $nonRedundantTrans = scalar(keys %$transcript);
	my $percent_mapped = sprintf '%.2f', 100*($count_mapped/$nonRedundantTrans);
	my $percent = sprintf '%.2f', 100*($nonRedundantTrans/$countTranscripts);

	print SUM
		"Totat transcripts in database:  ", $countTranscripts, "\n". 
		"Number of transcripts after cleaning out redundancy ",$nonRedundantTrans," ($percent%)\n".
		"Number of transcripts mapped to external refernece: $count_mapped ($percent_mapped%)\n\n";
	
		foreach my $key (sort {$annotations{$b} <=> $annotations{$a} }keys %annotations) {
			
			if ($annotations_mapped{$key}) {
				print SUM "Number of ", $key," annotations ", $annotations{$key}." ($annotations_mapped{$key})","\n";
			} else {
				print SUM "Number of ", $key," annotations ", $annotations{$key}." (0)\n";
			}
		}

	close(SUM);

	print "Output files written to directory ............. $work_dir\n";
}


####################################################################
	#
	# GET TRANSCRIPTS FROM SQLITE DATABASE
	#
####################################################################

sub transcripts_from_sqldb {

	my ($tbl, $external_mapping_hash) = @_;

	my %transcriptHoH = (); 			# collect sequences from SQLite DB

	my $dsn = "dbi:SQLite:dbname=$sqlitedb";
	my %attr = ( RaiseError => 1 );
	my $dbh = DBI->connect($dsn, $user, $password, \%attr) or die "Can't connect to database: $DBI::errstr";

	my $query = "	SELECT DISTINCT a.tr_stable_id, a.chr, a.start, a.start_codon, a.dist_to_aTIS, a.aTIS_call, a.annotation, a.peak_shift, a.SNP, a.aa_seq, b.biotype, c.exon_coverage
 			FROM $tbl a JOIN TIS_1 b JOIN tr_translation c
 			WHERE a.tr_stable_id = b.stable_id AND b.transcript_id = c.transcript_id";

 	my $sth = $dbh->prepare($query);
	$sth->execute();
	print "Extracting transcripts form SQLite database.\n";

	while ( my ($tr_stable_id, $chr, $start, $start_codon, $dist_to_aTIS, $aTIS_call, $annotation, $peak_shift, $snp, $aa_seq, $biotype, $exon_coverage) = $sth->fetchrow()) {
		
		$countTranscripts++;	
		$aa_seq =~ s/\*//g;
		next if (length($aa_seq) < $minSeqLength);
		next if ($aTIS_call eq "no_data");

		my $enID = $tr_stable_id."_".$chr."_".$start."_".$annotation;
		$start_codon = $aTIS_call." ".$start_codon." ".$dist_to_aTIS." ".$biotype." ".$exon_coverage;
		if ($snp) {$start_codon.=" ".$snp}

		unless ($transcriptHoH{$aa_seq}) {
			if (grep $_ eq $annotation, @TopAnno) {
			
				if ($$external_mapping_hash{$tr_stable_id}{'ext_ref'}) {
					$transcriptHoH{$aa_seq}{'ext_ref'} = $$external_mapping_hash{$tr_stable_id}{'ext_ref'};
					$start_codon .= " BM";
				}
				if (defined $$external_mapping_hash{$tr_stable_id}{'other_ref'}) {$start_codon = $$external_mapping_hash{$tr_stable_id}{'other_ref'}." ".$start_codon}
			}
			
			$transcriptHoH{$aa_seq}{'enID'} = $enID;
			$transcriptHoH{$aa_seq}{'des'} 	= $start_codon;
			$transcriptHoH{$aa_seq}{'anno'} = $annotation;
			if ($snp) { $transcriptHoH{$aa_seq}{'snp'} = $snp;}

		} else {
			if ($$external_mapping_hash{$tr_stable_id}{'ext_ref'}) { 	# check if current sequence have a swissprot ID
				if (grep $_ eq $annotation, @TopAnno) {	# check if current sequence is a aTIS or 5UTR
					
					if ($annotationHash{$annotation} == $annotationHash{$transcriptHoH{$aa_seq}{'anno'}}) { # check if SNP

						if ($transcriptHoH{$aa_seq}{'snp'}) {
							next;
						} elsif ($snp) {

							$transcriptHoH{$aa_seq}{'ext_ref'} = $$external_mapping_hash{$tr_stable_id}{'ext_ref'};
							$transcriptHoH{$aa_seq}{'enID'} = $enID;

							if (defined $$external_mapping_hash{$tr_stable_id}{'other_ref'}) {$start_codon = $$external_mapping_hash{$tr_stable_id}{'other_ref'}." ".$start_codon}
							$transcriptHoH{$aa_seq}{'des'} 	= $start_codon." BM";
							$transcriptHoH{$aa_seq}{'anno'} = $annotation;
							$transcriptHoH{$aa_seq}{'snp'} = $snp;
						}
						
					} elsif ($annotationHash{$annotation} < $annotationHash{$transcriptHoH{$aa_seq}{'anno'}}) {

						if ($snp) {$transcriptHoH{$aa_seq}{'snp'} = $snp}
						$transcriptHoH{$aa_seq}{'ext_ref'} = $$external_mapping_hash{$tr_stable_id}{'ext_ref'};
						$transcriptHoH{$aa_seq}{'enID'} = $enID;

						if (defined $$external_mapping_hash{$tr_stable_id}{'other_ref'}) {$start_codon = $$external_mapping_hash{$tr_stable_id}{'other_ref'}." ".$start_codon}
						$transcriptHoH{$aa_seq}{'des'} 	= $start_codon." BM";
						$transcriptHoH{$aa_seq}{'anno'} = $annotation;
					}
					
				} else {
					if ($transcriptHoH{$aa_seq}{'ext_ref'}) { # Skip if sequence in hash has a swissprot ID
						next;
 					} elsif ($annotationHash{$annotation} < $annotationHash{$transcriptHoH{$aa_seq}{'anno'}}) {  

						delete $transcriptHoH{$aa_seq};
						
						if (defined $$external_mapping_hash{$tr_stable_id}{'other_ref'}) {$start_codon = $$external_mapping_hash{$tr_stable_id}{'other_ref'}." ".$start_codon}
						$transcriptHoH{$aa_seq}{'ext_ref'} = $$external_mapping_hash{$tr_stable_id}{'ext_ref'};
						$transcriptHoH{$aa_seq}{'enID'} = $enID;
						$transcriptHoH{$aa_seq}{'des'} 	= $start_codon;
						$transcriptHoH{$aa_seq}{'anno'} = $annotation;
						if ($snp) {$transcriptHoH{$aa_seq}{'snp'} = $snp;}
						
					} elsif ($annotationHash{$annotation} == $annotationHash{$transcriptHoH{$aa_seq}{'anno'}}) { # check if SNP
						if ($transcriptHoH{$aa_seq}{'snp'}) {
							next;
						} elsif ($snp) {

							$transcriptHoH{$aa_seq}{'ext_ref'} = $$external_mapping_hash{$tr_stable_id}{'ext_ref'};
							$transcriptHoH{$aa_seq}{'enID'} = $enID;
							if (defined $$external_mapping_hash{$tr_stable_id}{'other_ref'}) {$start_codon = $$external_mapping_hash{$tr_stable_id}{'other_ref'}." ".$start_codon}
							$transcriptHoH{$aa_seq}{'des'} 	= $start_codon." BM";
							$transcriptHoH{$aa_seq}{'anno'} = $annotation;
							$transcriptHoH{$aa_seq}{'snp'} = $snp;
						}
					}
				}
				
			} else { 	# if current sequence is not aTIS or 5UTR keep sequence with best annotation rank
				if ($transcriptHoH{$aa_seq}{'ext_ref'}) {
					next;
				} elsif ($annotationHash{$annotation} < $annotationHash{$transcriptHoH{$aa_seq}{'anno'}}) {
					delete $transcriptHoH{$aa_seq};
					
					$transcriptHoH{$aa_seq}{'enID'} = $enID;
					$transcriptHoH{$aa_seq}{'des'} 	= $start_codon;
					$transcriptHoH{$aa_seq}{'anno'} = $annotation;
					if ($snp) {$transcriptHoH{$aa_seq}{'snp'} = $snp;}
					
				} elsif ($annotationHash{$annotation} == $annotationHash{$transcriptHoH{$aa_seq}{'anno'}}) {

					if ($transcriptHoH{$aa_seq}{'snp'}) {
						next;
					} elsif ($snp) {
						$transcriptHoH{$aa_seq}{'enID'} = $enID;
						$transcriptHoH{$aa_seq}{'des'} 	= $start_codon;
						$transcriptHoH{$aa_seq}{'anno'} = $annotation;
						$transcriptHoH{$aa_seq}{'snp'} = $snp;
					}	
				}
			}
		}
	}

	$sth->finish();
	$dbh->disconnect();

	unless (%transcriptHoH) {print "\nNo transcript found in table $tbl in database $sqlitedb"; exit}
	print "Total non duplicate sequences ",scalar(keys %transcriptHoH),".\n";

	return remove_redundancy(\%transcriptHoH);
}


sub remove_redundancy {

	my $transcript = shift;

	my %non_redundant;
	my $counter = my $percentage = 0;
	my $rows = keys (%$transcript);
	
	print "Removing redundant sequences. Please wait...... \n";
	foreach my $key1 (sort keys %$transcript) {

		my $flag = 0; 	# flag when current sequence is added to %non_redundant
		my $seq1 = substr $key1, 1;	# remove the first M of the sequence
		if (%non_redundant) {
			foreach my $key2 (keys %non_redundant) {
				my $seq2 = substr $non_redundant{$key2}{'seq'}, 1;	# remove the first M of the sequence
				
				if (index($seq1, $seq2) > -1) {

					delete $non_redundant{$key2};
					$non_redundant{$$transcript{$key1}{'enID'}}{'seq'} = $key1;
					
					if (defined $$transcript{$key1}{'ext_ref'}) {
						$non_redundant{$$transcript{$key1}{'enID'}}{'ext_ref'} = $$transcript{$key1}{'ext_ref'};
					}
					
					$non_redundant{$$transcript{$key1}{'enID'}}{'anno'} = $$transcript{$key1}{'anno'};
					$non_redundant{$$transcript{$key1}{'enID'}}{'des'} = $$transcript{$key1}{'des'};
					$flag = 1;
					
				} elsif (index($seq2, $seq1) > -1) {	# If $seq1 is a subsequence of a sequence in %non_rudendant drop it
					$flag = 1;
					last;
				}	
			}

			if ($flag == 0) {	# if not a subsequence of any sequence in %non_redundant add it into the hash

				$non_redundant{$$transcript{$key1}{'enID'}}{'seq'} = $key1;
				if (defined $$transcript{$key1}{'ext_ref'}) {
					$non_redundant{$$transcript{$key1}{'enID'}}{'ext_ref'} = $$transcript{$key1}{'ext_ref'};
				}
				$non_redundant{$$transcript{$key1}{'enID'}}{'anno'} = $$transcript{$key1}{'anno'};
				$non_redundant{$$transcript{$key1}{'enID'}}{'des'} = $$transcript{$key1}{'des'};
			}

		} else {
		
			if (defined $$transcript{$key1}{'ext_ref'}) {
				$non_redundant{$$transcript{$key1}{'enID'}}{'ext_ref'} = $$transcript{$key1}{'ext_ref'};
			}

			$non_redundant{$$transcript{$key1}{'enID'}}{'seq'} = $key1;
			$non_redundant{$$transcript{$key1}{'enID'}}{'anno'} = $$transcript{$key1}{'anno'};
			$non_redundant{$$transcript{$key1}{'enID'}}{'des'} = $$transcript{$key1}{'des'};
		}

		#if (int($counter++ / $rows * 100) > $percentage) {
		#	printf STDERR "\r%3d%%", $percentage ;
		#	$percentage = int($counter++ / $rows * 100) ;
		#}
	}

	$num_non_Red_trans = keys %non_redundant;
	print "\nRedundancy removal completed. $num_non_Red_trans non redundant sequences left\n";

	return (\%non_redundant);
}

####################################################################
	#
	# ERROR HANDLING
	#
####################################################################

	# Check required parameters
sub uninitialized_param {
	
	my ($v) = @_;
	not ( defined($v) and length $v );
}

	# Handle path to files
sub path {

	my ($file, $dir) = @_;
	return (File::Spec->catfile( $dir, $file));
}

	# File Error handling
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
