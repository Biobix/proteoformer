$|=1;
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
use Data::Dumper;
use File::Which;

# ---------------------------------------------------------------------
	##	GLOBAL VARIABLES

my $user 			  = "";		# User name for sqlite database
my $password 		  = "";		# Password for sqlite database
my $total_tr 		  = 0;
my $total_gene		  = 0;
my $no_blast_tr 	  = 0;
my $num_non_Red_trans = 0;

my $blastdb;		# Usearch/blastp formatted database
my $result_db;		# SQLite source database
my $tis_ids;		# TIS ID to generate table from
my $mapping;		# Ensembl database name to download biomart mapped transcripts 
my $mflag;			# Flag for biomart mappings, 1 = remote download, 2 = local file, 3 = sequence based mapping, 4 = no mapping but remove redundancy, 5 = no mapping and no redundancy removal.
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
my $translation_db;	# FASTA file of non redundant derived translation products OR name for output PEFF file
my $var_file;		# File to store SNP and indel info for sequences in protein database
my $tis_call;		# Allow annotated TIS that do not pass the TIS calling algorithm in the databases [Y or N]
my $db_config_version; 	# Ensembl databset confirguration version
my $external_ref;	# External reference in biomart to map transcripts to
my $tmp;			# temporary folder
my $peff;           # Generate PEFF instead of FASTA format [Y or N] (default N)

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
	'mapping=s'				=> \$mapping,
	'db_config_version=f'	=> \$db_config_version,
	'external_ref=s'		=> \$external_ref,
	'extra_info=s'			=> \$extra_info,
	'gapopen=f'	 			=> \$gapopen,
	'gapextend=f'	 		=> \$gapextend,
	'word_size=i'	 		=> \$word_size,
	'matrix=s'	 			=> \$matrix,
	'translation_db=s'		=> \$translation_db,
	'var_file=s'			=> \$var_file,
	'tis_call=s'			=> \$tis_call,
    'peff=s'                => \$peff,
);

# ---------------------------------------------------------------------
	## EXECUTION

## Command line
# perl generate_translation_db.pl -blast_db /path/to/blast_db -result_db results.db -tis_ids 1 -blast_pgm ublast -mflag 0 -external_ref uniprot_swissprot_accession -mapping_db mmusculus_gene_ensembl -num_threads 3 -work_dir new -tis_call Y
# mmusculus_gene_ensembl
# mart_export_mm.txt

# perl generate_translation_db.pl -blast_db /path/to/blast_db -result_db results.db -tis_ids 1 -blast_pgm ublast -mflag 1 -mapping_db mart_export_mm.txt -num_threads 3 -work_dir working_dir -tis_call Y -tmp temporary_dir -translation_db file_name
# mmusculus_gene_ensembl
# mart_export_mm.txt

# perl generate_translation_db.pl -result_db SQLite/results.db -tis_ids 1 -tis_call Y
#####


my $CWD = getcwd();
if (!($work_dir)) {
	$work_dir = $CWD;
} elsif (!-d "$work_dir") { 
	system ("mkdir ". $work_dir);
} 

print STDOUT "\n";
if ($work_dir) {
	$work_dir = abs_path($work_dir);
	print STDOUT "The following working directory : $work_dir\n";
}

my $TMP  = ($ENV{'TMP'}) ? $ENV{'TMP'} :  "$work_dir/tmp"; # (1) get the TMP environment variable, (2) get the $tmpfolder variable, (3) get current_working_dir/tmp
if (!-d "$TMP") { 
	system ("mkdir ".$TMP);
} 
print STDOUT "The following tmpfolder is used: $TMP\n";

if ($result_db) {
	print STDOUT "SQLite database containing transcripts : $result_db\n";
} 
if ($tis_ids) {
	print STDOUT "Database is generated for TIS ids : $tis_ids\n";
}

if ($tis_call) {
	print STDOUT "Allow annotated transcripts that do not pass the TIS calling algorithm: $tis_call\n";
} else {
	$tis_call = "Y";
	print STDOUT "Allow annotated transcripts that do not pass the TIS calling algorithm: $tis_call\n";
}

if ($mslength) {
    print STDOUT "Minimum Word size: $mslength\n";
} else {
    $mslength=6;
    print STDOUT "Minimum sequence length	: $mslength\n";
}
if ($evalue) {
    print STDOUT "Blast e-value : $evalue\n";
} else {
    $evalue=1e-10;
    print STDOUT "Blast e-value : $evalue\n";
}

if ($min_blast_length) {
    print STDOUT "Minimum sequence length allowed for Blast search	: $min_blast_length\n";
} else {
    $min_blast_length = 32;
    print STDOUT "Minimum sequence length allowed for Blast search : $min_blast_length\n";
}

if ($identity) {
    print STDOUT "Blast Identity value : $identity%\n";
} else {
    $identity=75;
    print STDOUT "Blast Identity value : $identity%\n";
}

if ($coverage) {
    print STDOUT "Minimum percentage of identical positions : $coverage\n";
} else {
    $coverage=30;
    print STDOUT "Minimum percentage of identical positions : $coverage\n";
}

if ($word_size) {
    print STDOUT "Minimum Word size: $word_size\n";
} else {
    $word_size = 3;
    print STDOUT "Minimum Word size: $word_size\n";
}

if ($gapopen) {
    print STDOUT "Cost of gap open	: $gapopen\n";
} else {
    $gapopen = 11;
    print STDOUT "Cost of gap open	: $gapopen\n";
}

if ($gapextend) {
    print STDOUT "Gap extension penalty	: $gapextend\n";
} else {
    $gapextend = 1;
    print STDOUT "Gap extension penalty	: $gapextend\n";
}

if ($matrix) {
    print STDOUT "Minimum Word size: $matrix\n";
} else {
    $matrix="BLOSUM62";
    print STDOUT "Matrix blast search matrix : $matrix\n";
}

if ($mflag) {

	if ($mflag == 1) {
		if ($mapping) {
			print STDOUT "Ensembl database for remote mapping : $mapping\n";
			if ($db_config_version) {
				print STDOUT "Ensembl dataset configuration version : $db_config_version\n";
			} else {
				print STDOUT "Dataset configuration version is require.\n";
				exit;
			}
			
			if ($external_ref) {
			
				print STDOUT "External reference to map transcript with : $external_ref\n";
				if ($extra_info) {
					print STDOUT "External reference to map transcript with : $extra_info\n";
				}
				
			} else {
				print STDOUT "Ensembl external attribute for biomart mapping is required.\n";
				exit;
			}

		} else {
			print STDOUT "Mapping database require!\n";
			print STDOUT "If mflag = 1, the name of the biomart database to download mapping to external Id is required.\n";
			exit;
		}
		
	} elsif ($mflag == 2) {
		if ($mapping) {
			if (-e $mapping) {
				print STDOUT "Locale file containing ensembl to Swissprot mapping : $mapping\n";
			} else {
				print STDOUT "No such file $mapping.\tEnsure the file exist and you have the required permission.\n";
				exit;
			}
		} else {
			print STDOUT "If mflag = 2, a comma seperated file for the Transcript id to external reference file is require.\n";
			print STDOUT "See the Readme file more information.\n";

		}
	} elsif ($mflag == 3) {
	
		print STDOUT "Sequence based mapping of transcripts to canonical database by blast search.\n";
		
		if ($blast_pgm) {
			print STDOUT "Blast program used for mapping : $blast_pgm\n"; 
		} else {
			$blast_pgm = "ublast";
			print STDOUT "Blast program used for mapping : $blast_pgm\n"; 
		}
        
        if ($blastdb) {
            print STDOUT "The blast database is : $blastdb\n";
        } else {
            print STDOUT "No blast database supplied. Ensure you have choose the no blast search option.\n";
        }
		
	} elsif ($mflag == 4) {
		print STDOUT "Derived translation product  database will not be mapped to any canonical information. Redundancy will be removed. \n";
    } elsif ($mflag == 5) {
        print STDOUT "Derived translation product database will not be mapped to any canonical information. Redundancy will not be removed. \n";
    }

} else {
	print STDOUT "Derived translation product database will not be mapped to any canonical information.\n";
    $mflag=4;
}

if ($peff){
    if ($peff ne "Y" and $peff ne "N"){
        print STDOUT "Peff argument should be 'Y' or 'N'!\n";
        die;
    } else {
        if ($peff eq 'N'){
            print STDOUT "The program will generate : FASTA\n";
        } elsif ($peff eq 'Y'){
            print STDOUT "The program will generate : PEFF\n";
            $mflag = 5;
            print STDOUT "Mflag automatically turned to '5' as no mapping is yet available for PEFF file generation.\n";
        }
    }
} else {
    #Default peff: N
    $peff = 'N';
    print STDOUT "The program will generate : FASTA\n";
}


	# get arguments from SQLite DB
my $dsn_results = "DBI:SQLite:dbname=$result_db";
my $dbh_results = dbh($dsn_results,$user,$password);
my ($run_name,$species,$mapper,$nr_of_cores)=get_input_vars($dbh_results);
print STDOUT "Number of cores used: $nr_of_cores\n";
print STDOUT "\n";

my %annotation = ( "aTIS"=>1,"5UTR"=>3,"CDS"=>2,"ntr"=>4,"3UTR"=>5);
my %top_anno = ("aTIS","5UTR");

print STDOUT "Annotations and their rankings (with 1 the most important) \n";
foreach my $key (sort {$annotation{$a} <=> $annotation{$b}} keys %annotation) {
    print STDOUT "$key\t$annotation{$key}\n";
}

	# If Mapping to canonical database is allowed 
my $external_mapping = {};
if ($mflag) {

	if ($mflag == 1) {	# Id based mapping
	
		my $mapping_file = remote_biomart_mapping($mapping,$db_config_version,$external_ref,$extra_info);
		print STDOUT "File containing maping information :	$mapping_file\n";
		$external_mapping = biomart_mapping($mapping_file);

	} elsif ($mflag == 2)  {
		$external_mapping = biomart_mapping($mapping);
	} 
}

STDOUT->flush();

my @idsref = get_analysis_ids($dbh_results,$tis_ids);  #$tis_ids is input variable

foreach (@idsref) {	# generate translation db for selected tis_ids 
    if ($peff eq 'N'){
        generate_trans_db($_);
    } else {
        generate_peff($_, $dsn_results, $user, $password, $TMP);
    }
}


#############################################
	#
	# SUBS
	#
#############################################

sub generate_peff {
    
    my $tis_id = $_[0];
    my $dsn = $_[1];
    my $us = $_[2];
    my $pw = $_[3];
    my $TMP = $_[4];
    
    my $startRun = time();
    my $table = "TIS_".$tis_id."_transcripts";
    my @tis = split('_', $tis_id); #SPLIT TIS id from SNP and INDEL underscore info
    my $ref_tis = \@tis;
    
    #Check if program can find clustalOmega
    my $tool_name = "clustalo";
    my $tool_path = which($tool_name);
    if (!-e $tool_path){
        print "ERROR: Could not find Clustal Omega installation!\n";
        print "Please install Clustal Omega or contact admin\n";
        die;
    }
    
    #Create tmp folder if not exists
    my $tmp_folder = $TMP."/peff";
    if (!-d $tmp_folder){
        system("mkdir ".$tmp_folder);
    }
    
    print "Get underlying info for multiprocessing\n";
    my ($chr_sizes, $cores, $ens_db, $species) = get_multiprocess_info();
    
    #Convert species for PEFF
    my ($species_peff, $taxid) = convert_species($species, $tmp_folder);
    
    # Init multi core
    my $pm = new Parallel::ForkManager($cores);
    print "START multiprocessing\n";
    print "   Using ".$cores." core(s)\n   ---------------\n";
    
    for my $chr (keys %{$chr_sizes}){
        
        ### Start parallel process
        $pm->start and next;
        
        ### PEFF generation per chr
        generate_peff_per_chr($table, $ref_tis, $chr, $dsn, $us, $pw, $tmp_folder, $ens_db, $species_peff, $taxid);
        
        ### Finish
        print "* Finished chromosome ".$chr."\n";
        $pm->finish;
    }
    
    # Finish all subprocesses
    $pm->wait_all_children;
    print "\n";
    
    #Define output
    my $output_dir = $work_dir;
    unless ($translation_db) {
        unless (-d "$output_dir") { system ("mkdir ".$output_dir)}
        $translation_db =  path($species."_".$table.".peff",$output_dir);
    }
    system("rm -rf ".$translation_db);
    system("touch ".$translation_db);
    
    #Concat all chr tmp files
    for my $chr (keys %{$chr_sizes}){
        system("cat ".$translation_db." ".$tmp_folder."/out_".$chr.".peff >> ".$work_dir."/out.peff");
        system("mv ".$work_dir."/out.peff ".$translation_db);
        system("rm -rf ".$tmp_folder."/out_".$chr.".peff");
    }
    
    #Calculate the number of entries based on line count and the date
    my $number_of_entries = (readpipe("wc -l < ".$translation_db))/2;
    my $ymd = sub{sprintf '%04d-%02d-%02d',$_[5]+1900, $_[4]+1, $_[3]}->(localtime);
    
    #Print header into file
    open(my $FWheader, '>', $tmp_folder."/header.peff");
    my $header = "# PEFF 1.0\n# //\n# DbName=genericDB\n# DbSource=http://www.biobix.be/proteoform/\n# DbVersion=".$ymd."\n# HasAnnotationIdentifiers=true \n# ProteoformDb=true\n# Prefix=gen\n# NumberOfEntries=".$number_of_entries."\n# SequenceType=AA\n# GeneralComment=Proteogenomics application to generate protein sequences from RIBO-seq\n# //\n";
    print $FWheader $header;
    close($FWheader);
    
    #Cat header to rest of peff file
    system("cat ".$tmp_folder."/header.peff ".$translation_db." > ".$work_dir."/out.peff");
    system("mv ".$work_dir."/out.peff ".$translation_db);
    
    #Remove tmp folder
    #system("rm -rf ".$TMP."/peff");
    
    
    
    
    
    #### TEST for  (tr_stable_id='ENST00000000412' or tr_stable_id='ENST00000429644');
    # Delete in get transcripts_with_orfs (2x) and transcript_gene_id_per_chr
    
    
    
    timer($startRun);
    print " ---done---\n\n";
    
    return;
}

sub generate_peff_per_chr{
    
    my $table = $_[0];
    my $tis = $_[1];
    my $chr = $_[2];
    my $dsn = $_[3];
    my $us = $_[4];
    my $pw = $_[5];
    my $tmp_folder = $_[6];
    my $ens_db = $_[7];
    my $species_peff = $_[8];
    my $taxid = $_[9];
    
    #Make dbh
    my $dbh_results = dbh($dsn, $us, $pw);
    
    #Get all transcripts with orfs in out of transcripts table (unique argument in sqlite)
    my ($orfs_per_transcript, $transcript2geneid) = get_transcripts_with_orfs($dbh_results, $table, ${$tis}[0], $chr, $ens_db);
    
    #Group orfs in proteoform structures
    my $proteoform_structs = group_orfs($orfs_per_transcript);
    
    #Search base ORF per proteoform structure (i.e. 'canonical' ORF where other proteoforms will be compared to)
    my $base_orfs = search_base_orfs($proteoform_structs, $orfs_per_transcript);
    
    #Align ORFs and construct proteoform structure
    my ($simpleVars, $complexVars, $proteoforms, $proteoformGenerals) = construct_proteoforms($proteoform_structs, $base_orfs, $tmp_folder);
    
    #Write proteoforms to chr tmp peff file
    write_chr_tmp_output($proteoforms, $proteoformGenerals, $simpleVars, $complexVars, $chr, $tmp_folder, $species_peff, $transcript2geneid, $proteoform_structs, $taxid);
    
    #Remove all tax files again
    #system("rm -rf ".$tmp_folder."/readme.txt");
    #system("rm -rf ".$tmp_folder."/gc.prt");
    #system("rm -rf ".$tmp_folder."/*.dmp");
    #system("rm -rf ".$tmp_folder."/taxdump.tar");
    
    return;
}

sub write_chr_tmp_output {
    
    #Catch
    my $proteoforms = $_[0];
    my $proteoformGenerals = $_[1];
    my $simpleVars = $_[2];
    my $complexVars = $_[3];
    my $chr = $_[4];
    my $tmp_folder = $_[5];
    my $species_peff = $_[6];
    my $transcript2geneid = $_[7];
    my $proteoform_structs = $_[8];
    my $taxid = $_[9];
    
    #Open chromosomal tmp peff file
    my $chr_tmp_file = $tmp_folder."/out_".$chr.".peff";
    system("rm -rf ".$chr_tmp_file);
    open(my $FW, '>', $chr_tmp_file);
    
    #Go over all transcript id's
    for my $transcript_id (keys(%{$proteoforms})){
        
        #Get gene name
        my $gene_name = $transcript2geneid->{$transcript_id}->{'description'};
        
        #Go over all structure id's of transcript id (these form the entries for the PEFF
        for my $struct_id (sort keys(%{$proteoforms->{$transcript_id}})){
            
            #Convert gene name to protein name
            my $pname = get_pname($gene_name, $proteoform_structs, $transcript_id, $struct_id);
            
            #Construct variant strings
            my $simpleVarStr = "";
            my $complexVarStr = "";
            my $var_number = 1;
            my $var_number_hash = {};
            #Simple variations
            for my $var_description (sort keys(%{$simpleVars->{$transcript_id}->{$struct_id}})){
                $simpleVarStr = $simpleVarStr."(".$var_number.":".$var_description.")";
                $var_number_hash->{$var_description} = $var_number;
                $var_number++;
            }
            if ($simpleVarStr ne ""){
                $simpleVarStr = "\\VariantSimple=".$simpleVarStr." ";
            }
            #Complex variations
            for my $var_description (sort keys(%{$complexVars->{$transcript_id}->{$struct_id}})){
                $complexVarStr = $complexVarStr."(".$var_number.":".$var_description.")";
                $var_number_hash->{$var_description} = $var_number;
                $var_number++;
            }
            if ($complexVarStr ne ""){
                $complexVarStr = "\\VariantComplex=".$complexVarStr." ";
            }
            
            #Sort orf ids for printing sequence: get minimum var number for each orf
            my %min_var_number;
            for my $orf_id (keys %{$proteoforms->{$transcript_id}->{$struct_id}}){
                my $min=100000;
                for my $simpleVarDescr (sort keys %{$proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'simpleVars'}}){
                    if($var_number_hash->{$simpleVarDescr}<$min){
                        $min = $var_number_hash->{$simpleVarDescr};
                    }
                }
                for my $complexVarDescr (sort keys %{$proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'complexVars'}}){
                    if($var_number_hash->{$complexVarDescr}<$min){
                        $min = $var_number_hash->{$complexVarDescr};
                    }
                }
                $min_var_number{$orf_id}=$min;
            }
            #Do sorting of proteoforms based on (1) minimum amount of variations (2) minimum var number hash
            my @sorted_orf_ids = sort {((scalar(keys %{$proteoforms->{$transcript_id}->{$struct_id}->{$a}->{'simpleVars'}})+scalar(keys %{$proteoforms->{$transcript_id}->{$struct_id}->{$a}->{'complexVars'}})) <=> (scalar(keys %{$proteoforms->{$transcript_id}->{$struct_id}->{$b}->{'simpleVars'}})+scalar(keys %{$proteoforms->{$transcript_id}->{$struct_id}->{$b}->{'complexVars'}}))) || $min_var_number{$a} cmp $min_var_number{$b}} keys %{$proteoforms->{$transcript_id}->{$struct_id}};
            
            #Construct proteoforms string
            my $proteoformStr = "";
            for my $orf_id (@sorted_orf_ids){
                my $varPartsStr = "";
                #Add simple variation numbers
                for my $simpleVarDescr (sort keys %{$proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'simpleVars'}}){
                    $varPartsStr = $varPartsStr.$var_number_hash->{$simpleVarDescr}.",";
                }
                #Add complex variation numbers
                for my $complexVarDescr (sort keys %{$proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'complexVars'}}){
                    $varPartsStr = $varPartsStr.$var_number_hash->{$complexVarDescr}.",";
                }
                #Cut trailing comma
                $varPartsStr =~ s/,$//;
                #Add the canoical form to the beginning of the proteoform string
                if ($proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'} eq "canonical form"){
                    $proteoformStr = "(".$proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'accession'}."|".$proteoformGenerals->{$transcript_id}->{$struct_id}->{'AZstart'}."-".$proteoformGenerals->{$transcript_id}->{$struct_id}->{'AZstop'}."|".$varPartsStr."|".$proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'}.")".$proteoformStr;
                } else {
                    #Add proteoform to the end of proteoforms string
                    $proteoformStr = $proteoformStr."(".$proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'accession'}."|".$proteoformGenerals->{$transcript_id}->{$struct_id}->{'AZstart'}."-".$proteoformGenerals->{$transcript_id}->{$struct_id}->{'AZstop'}."|".$varPartsStr."|".$proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'}.")";
                }
            }
            
            #Build accession string
            my $accession_string = ">gen:".$transcript_id."-".$struct_id." \\PName=".$pname." \\GName=".$gene_name." \\NcbiTaxId=".$taxid." \\TaxName=".$species_peff." \\Length=".length($proteoformGenerals->{$transcript_id}->{$struct_id}->{'main_seq'})." ".$simpleVarStr.$complexVarStr." \\Proteoform=".$proteoformStr;
            
            #Write to output
            print $FW $accession_string."\n";
            print $FW $proteoformGenerals->{$transcript_id}->{$struct_id}->{'main_seq'}."\n";
        }
    }
    
    close($FW);
    
    return;
}

sub get_pname {
    
    #Catch
    my $gene_name = $_[0];
    my $proteoform_structs = $_[1];
    my $transcript_id = $_[2];
    my $struct_id = $_[3];
    
    #Init
    my $pname = "";
    my $aTIS_found='N';
    my $UTR_found='N';
    my $only_ntr='Y';
    my $only_cds='Y';
    
    #Search for aTIS in this structure id, then PName=GName.
    for my $orf (keys %{$proteoform_structs->{$transcript_id}->{$struct_id}}){
        if ($proteoform_structs->{$transcript_id}->{$struct_id}->{$orf}->{'annotation'} ne 'ntr'){
            $only_ntr='N';
        }
        if ($proteoform_structs->{$transcript_id}->{$struct_id}->{$orf}->{'annotation'} ne 'CDS'){
            $only_ntr='N';
        }
        if ($proteoform_structs->{$transcript_id}->{$struct_id}->{$orf}->{'annotation'} eq 'aTIS'){
            $pname = $gene_name;
            $aTIS_found='Y';
            return $pname;
        }
    }
    #If non aTIS found, search for 5UTR or 3UTR
    if ($aTIS_found eq 'N'){
        for my $orf (keys %{$proteoform_structs->{$transcript_id}->{$struct_id}}){
            if ($proteoform_structs->{$transcript_id}->{$struct_id}->{$orf}->{'annotation'} eq '5UTR'){
                $pname = "uORF of ".$gene_name;
                $UTR_found='Y';
                return $pname;
            }
            if ($proteoform_structs->{$transcript_id}->{$struct_id}->{$orf}->{'annotation'} eq '3UTR'){
                $pname = "3'UTR downstream ORF of ".$gene_name;
                $UTR_found='Y';
                return $pname;
            }
        }
    }
    #search in the other structure id's of the same transcript for the canonical orf
    if($only_ntr eq 'Y'){
        for my $other_struct_id (keys %{$proteoform_structs->{$transcript_id}}){
            if($other_struct_id==$struct_id){
                next;
            }
            for my $orf_id (keys %{$proteoform_structs->{$transcript_id}->{$other_struct_id}}){
                if($proteoform_structs->{$transcript_id}->{$other_struct_id}->{$orf_id}->{'annotation'} eq 'aTIS'){
                    $pname = "Product from non-coding part of ".$gene_name;
                    return $pname;
                }
            }
        }
    } elsif ($only_cds='Y') {
        for my $other_struct_id (keys %{$proteoform_structs->{$transcript_id}}){
            if($other_struct_id==$struct_id){
                next;
            }
            for my $orf_id (keys %{$proteoform_structs->{$transcript_id}->{$other_struct_id}}){
                if($proteoform_structs->{$transcript_id}->{$other_struct_id}->{$orf_id}->{'annotation'} eq 'aTIS'){
                    $pname = "Alternative reading frame product of ".$gene_name;
                    return $pname;
                }
            }
        }
    }
    
    $pname = "Undeterminable but in ".$gene_name;
    return $pname;
}


sub construct_proteoforms {
    
    #Catch
    my $proteoform_structs = $_[0];
    my $base_orfs = $_[1];
    my $tmp_folder = $_[2];
    
    #Init
    my $simpleVars = {};
    my $complexVars = {};
    my $proteoforms = {};
    my $proteoformGenerals = {};
    
    #For all transcripts
    for my $transcript_id (keys %{$proteoform_structs}){
        #For all proteoform structures
        for my $struct_id (keys %{$proteoform_structs->{$transcript_id}}){
            #Run over all orfs in proteoform struct
            my $var_count=0;
            my $proteoform_count=1; #Start at 1 as the base orf get number 0
            for my $orf_id (keys %{$proteoform_structs->{$transcript_id}->{$struct_id}}){
                #Check if this is the base ORF
                if ($orf_id == $base_orfs->{$transcript_id}->{$struct_id}){
                    #Add base orf to proteoforms
                    #Add proteoform number 1 to accession
                    $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'accession'} = $proteoform_structs->{$transcript_id}->{$struct_id}->{$orf_id}->{'tr_id'}."_pf".$struct_id.".0";
                    $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'simpleVars'} = {};
                    $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'complexVars'} = {};
                    $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'} = 'canonical form';
                    #All ORFs in the proteoform structure are defined as variations to the base ORF, so coordinates of start and stop of base ORF are needed for all ORFs. Also the sequence of the base ORF is saved as main sequence
                    $proteoformGenerals->{$transcript_id}->{$struct_id}->{'AZstart'} = 1; #As this is base ORF, these are the reference AZ coordinates
                    $proteoformGenerals->{$transcript_id}->{$struct_id}->{'AZstop'} = length($proteoform_structs->{$transcript_id}->{$struct_id}->{$orf_id}->{'aa_seq'});
                    $proteoformGenerals->{$transcript_id}->{$struct_id}->{'main_seq'} =$proteoform_structs->{$transcript_id}->{$struct_id}->{$orf_id}->{'aa_seq'};
                } else {
                    #Prepare input fasta for aligning base ORF with the other ORF
                    my $input_fasta = prepare_fasta($proteoform_structs, $transcript_id, $struct_id, $orf_id, $base_orfs, $tmp_folder);
                    #Execute clustalOmega
                    my $output_fasta = execute_clustalO($input_fasta, $transcript_id, $orf_id, $tmp_folder);
                    #Parse output
                    my $simpleVarsInOrf = {};
                    my $complexVarsInOrf = {};
                    ($simpleVars, $complexVars, $simpleVarsInOrf, $complexVarsInOrf, $var_count) = parse_clustal($output_fasta, $simpleVars, $complexVars, $transcript_id, $struct_id, $orf_id, $input_fasta, $var_count);
                    #Add proteoform to proteoform hash based on vars
                    ($proteoforms, $proteoform_count) = add_proteoform($simpleVars, $complexVars, $transcript_id, $struct_id, $orf_id, $proteoform_count, $simpleVarsInOrf, $complexVarsInOrf, $proteoforms, $proteoform_structs, $base_orfs);
                }
            }
        }
    }
    
    return ($simpleVars, $complexVars, $proteoforms, $proteoformGenerals);
}

sub add_proteoform {
    
    #Catch
    my $simpleVars = $_[0];
    my $complexVars = $_[1];
    my $transcript_id = $_[2];
    my $struct_id = $_[3];
    my $orf_id = $_[4];
    my $proteoform_count = $_[5];
    my $simpleVarsInOrf = $_[6];
    my $complexVarsInOrf = $_[7];
    my $proteoforms = $_[8];
    my $proteoform_structs = $_[9];
    my $base_orfs = $_[10];
    
    #Add features to proteoforms hash
    $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'accession'} = $proteoform_structs->{$transcript_id}->{$struct_id}->{$orf_id}->{'tr_id'}."_pf".$struct_id.".".$proteoform_count;
    $proteoform_count++;
    $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'simpleVars'} = $simpleVarsInOrf;
    $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'complexVars'} = $complexVarsInOrf;
    
    #Get end coordinate of base ORF AA sequence
    my $end_base_orf = length($proteoform_structs->{$transcript_id}->{$struct_id}->{$base_orfs->{$transcript_id}->{$struct_id}}->{'aa_seq'});
    
    #Define a name for the optional tag in \proteoform of PEFF
    for my $var_descr (keys(%{$complexVarsInOrf})){
        my ($start_var, $end_var, $alt_seq) = split(/\|/,$var_descr);
        #Check for N-terminal extension
        if ($start_var eq '1' && $end_var eq '1' && length($alt_seq)>1){
            if (exists $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'}){
                $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'} = $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'}.", N-terminal extension";
            } else {
                $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'} = 'N-terminal extension';
            }
        } elsif($start_var eq $end_base_orf && $end_var eq $end_base_orf && length($alt_seq)>1){#Check for C-terminal extension
            if (exists $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'}){
                $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'} = $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'}.", C-terminal extension";
            } else {
                $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'} = 'C-terminal extension';
            }
        } elsif($start_var eq '1' && $alt_seq eq 'M'){#Check for N-terminal truncation
            if (exists $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'}){
                $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'} = $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'}.", N-terminal truncation";
            } else {
                $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'} = 'N-terminal truncation';
            }
        } elsif($end_var eq $end_base_orf && $alt_seq eq ""){#Check for C-terminal truncation
            if (exists $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'}){
                $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'} = $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'}.", C-terminal truncation";
            } else {
                $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'} = 'C-terminal truncation';
            }
        }
    }
    if (%{$simpleVarsInOrf}){
        if (exists $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'}){
            $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'} = $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'}.", proteoform with SAV";
        } else {
            $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'} = 'proteoform with SAV';
        }
    }
    #If none of the higher mentioned proteoform classes were found, give optional tag 'other'
    unless (exists $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'}){
        $proteoforms->{$transcript_id}->{$struct_id}->{$orf_id}->{'name'} = 'other';
    }
    
    return($proteoforms, $proteoform_count);
}

sub parse_clustal {
    
    #Catch
    my $output_fasta = $_[0];
    my $simpleVars = $_[1];
    my $complexVars = $_[2];
    my $transcript_id = $_[3];
    my $struct_id = $_[4];
    my $orf_id = $_[5];
    my $input_fasta = $_[6];
    my $var_count = $_[7];
    
    #Init
    my $simpleVarsInOrf = {};
    my $complexVarsInOrf = {};
    
    #Read in aligned strings
    open(my $fr, '<', $output_fasta) or die "Could not open ".$output_fasta." ".$!;
    my $base_header = <$fr>;
    my $base_aligned_seq = <$fr>;
    my $other_header = <$fr>;
    my $other_aligned_seq = <$fr>;
    close($fr);
    
    #Cut string in arrays of separate amino acids
    chomp($base_aligned_seq);
    chomp($other_aligned_seq);
    
    my @base_AAs = split(//, $base_aligned_seq);
    my @other_AAs = split(//, $other_aligned_seq);
    
    #Init, stat values during run over alignment
    my $pos_in_ref_seq = 1; #To keep the base position in the canonical sequence
    my $in_truncation = 'N';
    my $base_pos_start_truncation=0;
    my $base_pos_end_truncation=0;
    my $alt_seq = "";
    my $org_seq = "";
    my $last_common_base = "";
    my $writableAsSimpleVars = 'N';
    #Go over alignment and construct extensions, truncations, variations
    for (my $pos_base=0; $pos_base<scalar(@base_AAs); $pos_base++){
        if ($base_AAs[$pos_base] ne $other_AAs[$pos_base]){
            if ($in_truncation eq 'N'){
                #Start new truncation
                $in_truncation = 'Y';
                #For the first altered base, if there is a '-' in the base ORF, start from one base earlier, unless we are still at the start of the alignment
                if($base_AAs[$pos_base] eq '-' && $pos_base!=0){
                    $base_pos_start_truncation = $pos_in_ref_seq - 1;
                    $alt_seq = $last_common_base;
                } else {
                    $base_pos_start_truncation = $pos_in_ref_seq;
                }
                #Prime the paramater that checks if the variation is writable as simpleVars
                if ($base_AAs[$pos_base] ne '-' && $other_AAs[$pos_base] ne '-'){
                    $writableAsSimpleVars='Y';
                }
            }
            #Keep track of the altered and original sequence
            #At the same time, as soon as there is an indel, variation is not writable as a simpleVar (i.e. a combination of SAVs)
            if ($other_AAs[$pos_base] ne '-'){
                $alt_seq = $alt_seq.$other_AAs[$pos_base]
            } else {
                $writableAsSimpleVars='N';
            }
            if ($base_AAs[$pos_base] ne '-'){
                $org_seq = $org_seq.$base_AAs[$pos_base];
            } else {
                $writableAsSimpleVars='N';
            }
        } else {
            #Keep this as last common base
            $last_common_base = $base_AAs[$pos_base];
            if ($in_truncation eq 'Y'){
                #End truncation
                #If no bases in base ORF over truncation yet, take the first common base as well in the variation seq
                if($org_seq eq ""){
                    $base_pos_end_truncation = $pos_in_ref_seq;
                    if($other_AAs[$pos_base] ne '-'){
                        $alt_seq = $alt_seq.$other_AAs[$pos_base];
                    } else {
                        $writableAsSimpleVars = 'N';
                    }
                    if($base_AAs[$pos_base] ne '-'){
                        $org_seq = $org_seq.$base_AAs[$pos_base];
                    } else {
                        $writableAsSimpleVars = 'N';
                    }
                } else {
                    $base_pos_end_truncation = $pos_in_ref_seq-1;#Previous base AA was the last AA with difference
                }
                #Check if simple or complex variation
                if ($base_pos_start_truncation==$base_pos_end_truncation && length($alt_seq)==1 && length($org_seq)==1){
                    #Construct description
                    my $description = $base_pos_start_truncation."|".$alt_seq;
                    #Check if new array needs to be defined, i.e. not seen variation
                    unless (exists $simpleVars->{$transcript_id}->{$struct_id}->{$description}){
                        $simpleVars->{$transcript_id}->{$struct_id}->{$description} = [];
                        $var_count++; #also in the count, define next variation number
                    }
                    #Add
                    push(@{$simpleVars->{$transcript_id}->{$struct_id}->{$description}}, $var_count);
                    $simpleVarsInOrf->{$description} = $var_count;
                } else {
                    #Check if a complex variation is anyhow writable as a combination of simpleVars
                    if($writableAsSimpleVars eq 'Y'){
                        for(my $i=0; $i<length($alt_seq); $i++){
                            #Construct description
                            my $description = ($base_pos_start_truncation+$i)."|".substr($alt_seq, $i, 1);
                            #Check if new array needs to be defined, i.e. not seen variation
                            unless (exists $simpleVars->{$transcript_id}->{$struct_id}->{$description}){
                                $simpleVars->{$transcript_id}->{$struct_id}->{$description} = [];
                                $var_count++; #also in the count, define next variation number
                            }
                            #Add
                            push(@{$simpleVars->{$transcript_id}->{$struct_id}->{$description}}, $var_count);
                            $simpleVarsInOrf->{$description} = $var_count;
                        }
                    } else {
                        #Construct description
                        my $description = $base_pos_start_truncation."|".$base_pos_end_truncation."|".$alt_seq;
                        #Check if new array needs to be defined
                        unless (exists $complexVars->{$transcript_id}->{$struct_id}->{$description}){
                            $complexVars->{$transcript_id}->{$struct_id}->{$description} = [];
                            $var_count++;
                        }
                        #Add
                        push(@{$complexVars->{$transcript_id}->{$struct_id}->{$description}}, $var_count);
                        $complexVarsInOrf->{$description} = $var_count;
                    }
                }
                #Run values back to default
                $in_truncation = 'N';
                $base_pos_start_truncation = 0;
                $base_pos_end_truncation = 0;
                $alt_seq = "";
                $org_seq = "";
            }
        }

        #End of sequence and still in truncation mode: then this should added as variation as well
        if($pos_base==(scalar(@base_AAs)-1) && $in_truncation eq 'Y'){
            #End truncation
            #The position of the end coordinate in the base sequence depends on if the last base AA is '-' or a real base
            if($base_AAs[$pos_base] eq '-'){
                $base_pos_end_truncation = $pos_in_ref_seq-1;
            } else {
                $base_pos_end_truncation = $pos_in_ref_seq;
            }
            #Check if simple or complex variation
            if ($base_pos_start_truncation==$base_pos_end_truncation && length($alt_seq)==1){
                #Construct description
                my $description = $base_pos_start_truncation."|".$alt_seq;
                #Check if new array needs to be defined
                unless (exists $simpleVars->{$transcript_id}->{$struct_id}->{$description}){
                    $simpleVars->{$transcript_id}->{$struct_id}->{$description} = [];
                    $var_count++;
                }
                #Add
                push(@{$simpleVars->{$transcript_id}->{$struct_id}->{$description}}, $var_count);
                $simpleVarsInOrf->{$description} = $var_count;
            } else {
                #Check if a complex variation is anyhow writable as a combination of simpleVars
                if($writableAsSimpleVars eq 'Y'){
                    for(my $i=0; $i<length($alt_seq); $i++){
                        #Construct description
                        my $description = ($base_pos_start_truncation+$i)."|".substr($alt_seq, $i, 1);
                        #Check if new array needs to be defined, i.e. not seen variation
                        unless (exists $simpleVars->{$transcript_id}->{$struct_id}->{$description}){
                            $simpleVars->{$transcript_id}->{$struct_id}->{$description} = [];
                            $var_count++; #also in the count, define next variation number
                        }
                        #Add
                        push(@{$simpleVars->{$transcript_id}->{$struct_id}->{$description}}, $var_count);
                        $simpleVarsInOrf->{$description} = $var_count;
                    }
                } else {
                    #Construct description
                    my $description = $base_pos_start_truncation."|".$base_pos_end_truncation."|".$alt_seq;
                    #Check if new array needs to be defined
                    unless (exists $complexVars->{$transcript_id}->{$struct_id}->{$description}){
                        $complexVars->{$transcript_id}->{$struct_id}->{$description} = [];
                        $var_count++;
                    }
                    #Add
                    push(@{$complexVars->{$transcript_id}->{$struct_id}->{$description}}, $var_count);
                    $complexVarsInOrf->{$description} = $var_count;
                }
            }
        }
        #Keep track of the AA number in the base ORF sequence (i.e. +1 only if not a '-' in the base ORF sequence)
        if ($base_AAs[$pos_base] ne '-'){
            $pos_in_ref_seq++;
        }
    }
    #Check if the process ran over the whole alignment file
    my $filtered_seq = $base_aligned_seq;
    $filtered_seq =~ s/\-//g;
    if ($pos_in_ref_seq != length($filtered_seq)+1){
        print "ERROR in running over alignment for transcript ID ".$transcript_id.". ORF ID ".$orf_id."\n";
        print "Last seen position is ".$pos_in_ref_seq.", while length of alignment is ".scalar(@base_AAs)."\n";
        die;
    }
    
    #Clean up
    system("rm -rf ".$input_fasta);
    system("rm -rf ".$output_fasta);
    
    return($simpleVars, $complexVars, $simpleVarsInOrf, $complexVarsInOrf, $var_count);
}

sub execute_clustalO{
    
    #Catch
    my $input_fasta = $_[0];
    my $transcript_id = $_[1];
    my $orf_id = $_[2];
    my $tmp_folder = $_[3];
    
    #Define output fasta path
    my $output_fasta = $tmp_folder."/output_".$transcript_id."_".$orf_id.".fa";
    
    #Execute
    my $command = "clustalo -i ".$input_fasta." -o ".$output_fasta." --outfmt fa --force --wrap=1000000000"; #Force to overwrite, output in fasta format
    #print $command."\n";
    system($command);
    
    return $output_fasta;
}

sub prepare_fasta{
    
    #Catch
    my $proteoform_structs = $_[0];
    my $transcript_id = $_[1];
    my $struct_id = $_[2];
    my $orf_id = $_[3];
    my $base_orfs = $_[4];
    my $tmp_folder = $_[5];
    
    #Define path to input fasta file
    my $input_fasta = $tmp_folder."/input_".$transcript_id."_".$orf_id.".fa";
    
    #Open fasta file and write aa seqs of base and other ORF
    open(my $fw, '>', $input_fasta) or die "Could not open ".$input_fasta." ".$!."\n";
    print $fw ">base\n";
    print $fw $proteoform_structs->{$transcript_id}->{$struct_id}->{$base_orfs->{$transcript_id}->{$struct_id}}->{'aa_seq'}."\n";
    print $fw ">other\n";
    print $fw $proteoform_structs->{$transcript_id}->{$struct_id}->{$orf_id}->{'aa_seq'}."\n";
    close $fw;
    
    return $input_fasta;
}


sub search_base_orfs {
    
    #Catch
    my $proteoform_structs = $_[0];
    my $orfs_per_transcript = $_[1];
    
    #Init
    my $base_orfs;
    
    #For all transcripts
    for my $transcript_id (keys %{$proteoform_structs}){
        #For all proteoform structures
        for my $struct_id (keys %{$proteoform_structs->{$transcript_id}}){
            #Run over all ORFs present in the proteoform struct
            for my $orf_id (keys %{$proteoform_structs->{$transcript_id}->{$struct_id}}){
                #Check if already a tmp base ORF present
                if (defined $base_orfs->{$transcript_id}->{$struct_id}){
                    #Check for aTIS annotation in tmp base ORF and considered ORF
                    if ($proteoform_structs->{$transcript_id}->{$struct_id}->{$base_orfs->{$transcript_id}->{$struct_id}}->{'annotation'} ne 'aTIS' && $proteoform_structs->{$transcript_id}->{$struct_id}->{$orf_id}->{'annotation'} eq 'aTIS'){
                        #tmp base ORF is not aTIS but considered is, change considered aTIS to base ORF
                        $base_orfs->{$transcript_id}->{$struct_id} = $orf_id;
                    } elsif ($proteoform_structs->{$transcript_id}->{$struct_id}->{$base_orfs->{$transcript_id}->{$struct_id}}->{'annotation'} eq 'aTIS' && $proteoform_structs->{$transcript_id}->{$struct_id}->{$orf_id}->{'annotation'} eq 'aTIS'){
                        #If both tmp base ORF and considered ORF are aTIS annotated, take the longest sequence as base ORF
                        if (length($proteoform_structs->{$transcript_id}->{$struct_id}->{$base_orfs->{$transcript_id}->{$struct_id}}->{'aa_seq'})<length($proteoform_structs->{$transcript_id}->{$struct_id}->{$orf_id}->{'aa_seq'})){
                            $base_orfs->{$transcript_id}->{$struct_id} = $orf_id;
                        } elsif (length($proteoform_structs->{$transcript_id}->{$struct_id}->{$base_orfs->{$transcript_id}->{$struct_id}}->{'aa_seq'})==length($proteoform_structs->{$transcript_id}->{$struct_id}->{$orf_id}->{'aa_seq'})){
                            #For equal lengths, take the sequence without SNP id
                            if ($orfs_per_transcript->{$transcript_id}->{$base_orfs->{$transcript_id}->{$struct_id}}->{'SNP'} ne "" && $orfs_per_transcript->{$transcript_id}->{$orf_id}->{'SNP'} eq ""){
                                $base_orfs->{$transcript_id}->{$struct_id} = $orf_id;
                            }
                        }
                    } elsif ($proteoform_structs->{$transcript_id}->{$struct_id}->{$base_orfs->{$transcript_id}->{$struct_id}}->{'annotation'} ne 'aTIS' && $proteoform_structs->{$transcript_id}->{$struct_id}->{$orf_id}->{'annotation'} ne 'aTIS') {
                        #Both are not aTIS, take the longest sequence as base ORF
                        if (length($proteoform_structs->{$transcript_id}->{$struct_id}->{$base_orfs->{$transcript_id}->{$struct_id}}->{'aa_seq'})<length($proteoform_structs->{$transcript_id}->{$struct_id}->{$orf_id}->{'aa_seq'})){
                            $base_orfs->{$transcript_id}->{$struct_id} = $orf_id;
                        } elsif (length($proteoform_structs->{$transcript_id}->{$struct_id}->{$base_orfs->{$transcript_id}->{$struct_id}}->{'aa_seq'})==length($proteoform_structs->{$transcript_id}->{$struct_id}->{$orf_id}->{'aa_seq'})){
                            #For equal lengths, take the sequence without SNP id
                            if ($orfs_per_transcript->{$transcript_id}->{$base_orfs->{$transcript_id}->{$struct_id}}->{'SNP'} ne "" && $orfs_per_transcript->{$transcript_id}->{$orf_id}->{'SNP'} eq ""){
                                $base_orfs->{$transcript_id}->{$struct_id} = $orf_id;
                            }
                        }
                    } else { #Base ORF is at annotated TIS and considered ORF not
                        next;
                    }
                } else {
                    #For the first ORF, take this as the tmp base ORF
                    $base_orfs->{$transcript_id}->{$struct_id} = $orf_id;
                }
            }
        }
    }

    return $base_orfs;
}


sub group_orfs {
    
    #Catch
    my $orfs_per_transcript = $_[0];
    
    #Init
    my $proteoform_structs = {};
    
    #Construct proteoforms per transcript ID
    for my $transcript_id (keys %{$orfs_per_transcript}){
        #Init
        my $struct_count=1;
        for my $orf_id (keys %{$orfs_per_transcript->{$transcript_id}}){
            #Check if there is already a proteoform_struct for that transcript ID
            if ($proteoform_structs->{$transcript_id}){
                #Init boolean to check if new ORF falls in existing proteoform structure
                my $in_existing_structure = "N";
                #Run over existing structures and ORFs already in these structures
                for my $struct (keys %{$proteoform_structs->{$transcript_id}}){
                    for my $orf_in_struct (keys %{$proteoform_structs->{$transcript_id}->{$struct}}){
                        #ORFs in same proteoform struct if they have start or stop coordinate in common (C- or N-terminal truncation or SAV inside)
                        if ($orfs_per_transcript->{$transcript_id}->{$orf_id}->{'start'}==$proteoform_structs->{$transcript_id}->{$struct}->{$orf_in_struct}->{'start'} || $orfs_per_transcript->{$transcript_id}->{$orf_id}->{'stop'}==$proteoform_structs->{$transcript_id}->{$struct}->{$orf_in_struct}->{'stop'}){
                            $proteoform_structs->{$transcript_id}->{$struct}->{$orf_id} = $orfs_per_transcript->{$transcript_id}->{$orf_id};
                            $in_existing_structure = "Y";
                        }
                    }
                }
                #If not in existing structure, define new proteoform structure
                if ($in_existing_structure eq "N"){
                    $proteoform_structs->{$transcript_id}->{$struct_count}->{$orf_id} = $orfs_per_transcript->{$transcript_id}->{$orf_id};
                    $struct_count++;
                }
            } else {
                #Define first proteoform struct if no struct exist yet
                $proteoform_structs->{$transcript_id}->{$struct_count}->{$orf_id} = $orfs_per_transcript->{$transcript_id}->{$orf_id};
                $struct_count++;
            }
        }
        
        #Make sure the aTIS proteoform structure gets structure id 1
        my $copy_proteoform_structs={};
        for my $struct_id (keys %{$proteoform_structs->{$transcript_id}}){
            my $aTIS_found = 'N';
            for my $orf_id (keys %{$proteoform_structs->{$transcript_id}->{$struct_id}}){
                if ($proteoform_structs->{$transcript_id}->{$struct_id}->{$orf_id}->{'annotation'} eq 'aTIS'){
                    $aTIS_found = 'Y';
                    last;
                }
            }
            if ($aTIS_found eq 'Y'){
                #Put aTIS proteoform struct as struct id 1
                $copy_proteoform_structs->{$transcript_id}->{1} = $proteoform_structs->{$transcript_id}->{1};
                $proteoform_structs->{$transcript_id}->{1} = $proteoform_structs->{$transcript_id}->{$struct_id};
                $proteoform_structs->{$transcript_id}->{$struct_id} = $copy_proteoform_structs->{$transcript_id}->{1};
            }
        }
    }

    return $proteoform_structs;
}

sub get_transcripts_with_orfs{
    
    #Catch
    my ($dbh, $tbl, $tis, $chr, $ens_db) = @_;
    
    #Init
    my $orfs_per_transcript = {};
    my $var_tracker = {};

    #Get gene ID info of all transcripts and get all annotated aTIS ORFs
    my ($transcript2geneid,$annotated_tr) = transcript_gene_id_per_chr($dbh,$tbl,$chr,$ens_db);
    
    my $query_transcripts = "SELECT DISTINCT tr_stable_id FROM ".$tbl." WHERE chr='".$chr."';";
    my $sth_transcripts = $dbh->prepare($query_transcripts);
    $sth_transcripts->execute();
    my $transcripts = $sth_transcripts->fetchall_arrayref();

    my $query_orfs = "SELECT tr_stable_id, chr, start, start_codon, dist_to_aTIS, aTIS_call, annotation, peak_shift, SNP, INDEL, aa_seq, stop, strand FROM ".$tbl." WHERE chr='".$chr."';";
    my $sth_orfs = $dbh->prepare($query_orfs);
    $sth_orfs->execute();
    my $orfs = $sth_orfs->fetchall_arrayref();
    
    for my $orf (@$orfs) {
        my $orf_id=scalar(keys(%{$orfs_per_transcript->{$$orf[0]}}));
        foreach (@$orf) {$_ = '' unless defined}; #Change empty SQLite values to empty strings
        
        # if instructed to not keep aTIS with not enough coverage to call TIS
        if (uc($tis_call) eq "N") {next if (uc($$orf[5]) eq 'NO_DATA' or uc($$orf[5]) eq 'FALSE')}
        $$orf[10] =~ s/\*//g;
        next if (length($$orf[10]) < $mslength);	# skip if sequence is less than minimum allowed amino acid length
        
        #Skip if exact same sequence already exists for that transcript id
        my $skip_orf = 'N';
        for my $other_orf (keys %{$orfs_per_transcript->{$$orf[0]}}){
            if ($$orf[10] eq $orfs_per_transcript->{$$orf[0]}->{$other_orf}->{'aa_seq'}){
                $skip_orf = 'Y';
                last;
            }
        }
        if ($skip_orf eq 'Y'){
            next;
        }

        # Skip all non aTIS without SNP or indel information that corresponds to an annotated TIS i.e redundant non annotated TIS
        my $red_tis = 0;
        if ($$orf[6] ne 'aTIS' and $$orf[8] eq "" and $$orf[9] eq "") {
            if ($transcript2geneid->{$$orf[0]}) {
                my $gene = $transcript2geneid->{$$orf[0]}->{'gene'};
                foreach my $start1 (keys %{$annotated_tr->{$gene}}) {
                    if ($start1 == $$orf[2]) {$red_tis = 1}
                }
            }
        }
        next if ($red_tis == 1);
        
        # create unique transcript ID
        my $tr = $$orf[0]."_".$$orf[1]."_".$$orf[2]."_".$$orf[6];
        if ($$orf[8] ne "" or $$orf[9] ne "") {		# if snp or indel info exist for current record
            if ($var_tracker->{$tr}) {
                $var_tracker->{$tr}++;
                $tr = $tr."_".$var_tracker->{$tr};
            } else {
                $var_tracker->{$tr} = 1;
                $tr = $tr."_".$var_tracker->{$tr};
            }
        }

        #Save in hash
        $orfs_per_transcript->{$$orf[0]}->{$orf_id}->{'tr_stable_id'} = $$orf[0];
        $orfs_per_transcript->{$$orf[0]}->{$orf_id}->{'chr'} = $$orf[1];
        $orfs_per_transcript->{$$orf[0]}->{$orf_id}->{'start'} = $$orf[2];
        $orfs_per_transcript->{$$orf[0]}->{$orf_id}->{'start_codon'} = $$orf[3];
        $orfs_per_transcript->{$$orf[0]}->{$orf_id}->{'dist_to_aTIS'} = $$orf[4];
        $orfs_per_transcript->{$$orf[0]}->{$orf_id}->{'aTIS_call'} = $$orf[5];
        $orfs_per_transcript->{$$orf[0]}->{$orf_id}->{'annotation'} = $$orf[6];
        $orfs_per_transcript->{$$orf[0]}->{$orf_id}->{'peak_shift'} = $$orf[7];
        $orfs_per_transcript->{$$orf[0]}->{$orf_id}->{'SNP'} = $$orf[8];
        $orfs_per_transcript->{$$orf[0]}->{$orf_id}->{'INDEL'} = $$orf[9];
        $orfs_per_transcript->{$$orf[0]}->{$orf_id}->{'aa_seq'} = $$orf[10];
        $orfs_per_transcript->{$$orf[0]}->{$orf_id}->{'stop'} = $$orf[11];
        $orfs_per_transcript->{$$orf[0]}->{$orf_id}->{'strand'} = $$orf[12];
        $orfs_per_transcript->{$$orf[0]}->{$orf_id}->{'tr_id'} = $tr;
    }
    
    
    return $orfs_per_transcript, $transcript2geneid;
}

sub get_multiprocess_info{
    
    #Get arguments from arguments table
    my ($ensemblversion,$ens_db,$species,$IGENOMES_ROOT,$cores) = get_arguments($dbh_results);
    
    #Conversion for species terminology
    my $spec = ($species eq "mouse") ? "Mus_musculus" : ($species eq "human") ? "Homo_sapiens" : ($species eq "arabidopsis") ? "Arabidopsis_thaliana" : ($species eq "fruitfly") ? "Drosophila_melanogaster" : "";
    my $spec_short = ($species eq "mouse") ? "mmu" : ($species eq "human") ? "hsa" : ($species eq "arabidopsis") ? "ath" : ($species eq "fruitfly") ? "dme" : "";
    
    #Old mouse assembly = NCBIM37, new one is GRCm38. Old human assembly = GRCh37, the new one is GRCh38
    my $assembly = (uc($species) eq "MOUSE" && $ensemblversion >= 70 ) ? "GRCm38"
    : (uc($species) eq "MOUSE" && $ensemblversion < 70 ) ? "NCBIM37"
    : (uc($species) eq "HUMAN" && $ensemblversion >= 76) ? "GRCh38"
    : (uc($species) eq "HUMAN" && $ensemblversion < 76) ? "GRCh37"
    : (uc($species) eq "ARABIDOPSIS") ? "TAIR10"
    : (uc($species) eq "FRUITFLY" && $ensemblversion < 79) ? "BDGP5"
    : (uc($species) eq "FRUITFLY" && $ensemblversion >= 79) ? "BDGP6" : "";
    
    # Get chromosomes and correct coord_system_id
    my $chromosome_sizes; my $coord_system_id; my @ch;
    
    $chromosome_sizes = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Annotation/Genes/ChromInfo.txt";
    ## Get chromosome sizes and cDNA identifiers #############
    print "Getting chromosome sizes and cDNA to chromosome mappings ...\n";
    my %chr_sizes = %{get_chr_sizes($chromosome_sizes)};
    
    my $ref_to_chr_sizes = \%chr_sizes;
    $coord_system_id = get_coord_system_id($ens_db,$assembly);
    
    return ($ref_to_chr_sizes, $cores, $ens_db, $species);
}


########################
# GENERATE FASTA SUBS  #
########################


sub generate_trans_db {

	my $tis_id = shift;

	my $startRun = time();
	my $table = "TIS_".$tis_id."_transcripts";
    my @tis = split '_', $tis_id; #SPLIT TIS id from SNP and INDEL underscore info
    
    #Get TIS calling method for TIS ID
    my $TIS_calling_method = get_TIS_calling_method($dbh_results, $tis[0]);
    my $ens_db;
    
    #Get ens_db for PRICE and SPECtre
    if($TIS_calling_method ne "Yes"){
        $ens_db = get_ens_db($dbh_results);
    }
	
    print "Get transcript out of results DB\n";
    my ($transcript,$gene_transcript) = get_transcripts_from_resultdb($dbh_results,$table,$tis[0], $TIS_calling_method, $ens_db); #Use TIS id itself (stored in $tis[0])
	$total_tr = scalar(keys %$transcript);
	$total_gene = scalar(keys %$gene_transcript);
    
    my $non_redundant_transcript;
    if($mflag!=5){
        print "Removin redundant sequences from TIS fasta database\n";
        $non_redundant_transcript = remove_redundancy($transcript);
        my $total_non_red_tr = scalar(keys %$non_redundant_transcript);
    }
	
	# 1 = remote download, 2 = local file, 3 = sequence based mapping, 4 = no mapping.
	if ($mflag == 1 or $mflag == 2) {		#  ID based mapping
		$non_redundant_transcript = id_based_mapping($non_redundant_transcript,$external_mapping);
		
	} elsif ($mflag == 3) {		# Sequence based mapping
		
		my $file_to_blast = transcripts_to_blast($non_redundant_transcript, $tis_id);
		if ($blast_pgm eq "blastp") {
			$non_redundant_transcript = blastp($non_redundant_transcript, $file_to_blast);
		} elsif ($blast_pgm eq "ublast")  {
			$non_redundant_transcript = ublast($non_redundant_transcript, $file_to_blast);
		}
	}	

	my $output_dir = $work_dir;
	unless ($translation_db) {
		unless (-d "$output_dir") { system ("mkdir ".$output_dir)}
		$translation_db =  path($species."_".$table.".fasta",$output_dir);
	} 

	unless ($var_file) {
		unless (-d "$output_dir") { system ("mkdir ".$output_dir)}
		$var_file =  path($species."_".$table."_VAR.txt",$output_dir);
	}
    
    print "Write output\n";
    if ($mflag==5){
        write_output($transcript, $translation_db, $var_file);
    } else {
        write_output($non_redundant_transcript,$translation_db,$var_file);
    }
	
	timer($startRun);	# Get Run time
	print STDOUT "-- Done --\n";
}

##Get ens db
sub get_ens_db{
    
    #Catch
    my $dbh = $_[0];
    
    my $query = "SELECT value FROM arguments WHERE variable='ens_db';";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    
    my @result = $sth->fetchrow_array();
    my $ens_db = $result[0];

    return $ens_db;
}

##Get transcript calling method from arguments table
sub get_TIS_calling_method{
    
    #Catch
    my $dbh = $_[0];
    my $tis_id = $_[1];
    
    my $query = "SELECT TIS_calling FROM TIS_overview WHERE ID=".$tis_id.";";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    
    my @results = $sth->fetchrow_array();
    my $TIS_calling_method = $results[0];
    
    return $TIS_calling_method;
}


##------ WRITE OUTPUT -------##
sub write_output {

	my $transcript 	= $_[0];
	my $output 		= $_[1];
	my $var_file 	= $_[2];

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
		
		if ($transcript->{$tr}->{'others'}) {
			my @Others = uniq(@{$transcript->{$tr}->{'others'}});
			$desc = $desc." [".join("#",@Others)."]";
		}
		my $seq_obj = Bio::Seq->new(-display_id => $id, -desc => $desc, -seq => $transcript->{$tr}->{'seq'});
		$seq_out->write_seq($seq_obj);

		$annotations{$transcript->{$tr}->{'anno'}}++;
  	}

		# write Summary of external mapping to file
	my $percent_overall = sprintf '%.1f', 100*($count_mapped/$total_tr);
	my $nonRedundantTrans = scalar(keys %$transcript);
	my $percent_mapped = sprintf '%.1f', 100*($count_mapped/($nonRedundantTrans+1));
	my $percent = sprintf '%.1f', 100*($nonRedundantTrans/$total_tr);

	print STDOUT "Total number of genes with at least one TIS called ", $total_gene, "\n";
	print STDOUT "Total number of possible TIS predicted ", $total_tr, "\n";
    if($mflag==5){
        print STDOUT "Total number of sequences in redundant fasta db ",$nonRedundantTrans,"\n";
    } else {
        print STDOUT "Total number of non redundant sequence in fasta db ",$nonRedundantTrans," ($percent%)\n";
    }

	if ($mflag == 1 or $mflag == 2 or $mflag == 2) {
		print STDOUT "Number of transcripts mapped to external refernece: $count_mapped ($percent_mapped%)\n\n";
	}

	foreach my $key (sort {$annotations{$b} <=> $annotations{$a} }keys %annotations) {
		if ($annotations_mapped{$key}) {
			print STDOUT "Number of ", $key," annotations ", $annotations{$key}." ($annotations_mapped{$key})","\n";
		} else {
			print STDOUT "Number of ", $key," annotations ", $annotations{$key}."\n";
		}
	}

	# Write SNP and indel Information to VAR file
	open(F, ">".$var_file) or die "Cannot create file $var_file \n";
	print F "transcript\tSNP_info\tindel_info\n";
	foreach my $tr (sort keys %$transcript) {
		next if ($transcript->{$tr}->{'snp'} eq "" and $transcript->{$tr}->{'indel'} eq "");
		print F "$tr\t$transcript->{$tr}->{'snp'}\t$transcript->{$tr}->{'indel'}\n";
	}
	close F;

	print STDOUT "Generated database written to file $output\n";

}



sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}


##------ SEQUENCE MAPPING -------##


	# Blast Search
sub blastp {

	my $transcript = $_[0];
	my $trans_to_blast = $_[1];
	
	print STDOUT "Running Local blastp search...\n";
	my ($blast_rpt) = $trans_to_blast =~ /^(.*)\./;
	$blast_rpt = $blast_rpt.".bls";
			
	my $command = "blastp -query ".$trans_to_blast." -db ".$blastdb." -out ".$blast_rpt." -num_alignments 1 -num_descriptions 1 -evalue ".$evalue." -num_threads ".$nr_of_cores." -gapopen ".$gapopen." -gapextend ".$gapextend." -matrix ".$matrix." > /dev/null 2>&1";
	system($command);

	return blastp_parser($transcript, $blast_rpt);
}

sub ublast {

	my $transcript = $_[0];
	my $file_to_blast = $_[1];
	
	print STDOUT "Running Local ublast search...\n";
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
	#print "number of files:$num_files \n";
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
	
	print STDOUT "Total number of transcripts mapped to a canonical ID.\n";
	return $transcript;
}

##------ REMOVE REDUNDANCY -------##

sub remove_redundancy {
    
    my $transcript = $_[0];
    
    my $transcript_seq;
    
    my $count = 0;
    foreach my $tr (keys %$transcript) {
        
        my $seq1 = $transcript->{$tr}->{'seq'};
        $transcript->{$tr}->{'len'} = length($seq1);
        
        if ($transcript_seq) {
            if (exists $transcript_seq->{$seq1}) {
                $transcript_seq->{$seq1}->{$tr} = 1;
            } else {
                my $flag = 1;
                my $seq1_tmp = substr($seq1, 1);
                foreach my $seq2 (keys %$transcript_seq) {
                    my $seq2_tmp = substr($seq2, 1);
                    
                    if (index($seq1_tmp, $seq2_tmp) >= 0) {
                        foreach my $tr_tmp (keys %{$transcript_seq->{$seq2}}) {
                            $transcript_seq->{$seq1}->{$tr_tmp} = 1;
                        }
                        $transcript_seq->{$seq1}->{$tr} = 1;
                        delete $transcript_seq->{$seq2};
                        $flag = 0;
                    } elsif (index($seq2_tmp, $seq1_tmp) >= 0) {
                        $transcript_seq->{$seq2}->{$tr} = 1;
                        $flag = 0;
                        last;
                    }
                }
                
                # seq not a subseq of any in non redundant hash
                if ($flag == 1) {
                    $transcript_seq->{$seq1}->{$tr} = 1;
                }
            }
        } else {
            $transcript_seq->{$seq1}->{$tr} = 1;
        }
    }
    
    # select representative accession for longest non redundant sequence
    my $non_red_trans;
    foreach my $seq (keys %$transcript_seq) {
        my @acessions = (keys %{$transcript_seq->{$seq}});
        my $selected_acc = pop @acessions;
        foreach my $tr (@acessions) {
            if ($transcript->{$tr}->{'len'} > $transcript->{$selected_acc}->{'len'}) {
                $selected_acc = $tr;
            } elsif ($transcript->{$tr}->{'len'} == $transcript->{$selected_acc}->{'len'}) {
                if ($transcript->{$selected_acc}->{'snp'} ne "" and $transcript->{$selected_acc}->{'indel'} ne "") {
                    # Selected acession have SNP and INDEL info hence keep
                } else {
                    if ($transcript->{$tr}->{'snp'} ne "" and $transcript->{$tr}->{'indel'} ne "") {
                        $selected_acc = $tr;    # current accession has SNP and INDEL  while select doesn't
                    } else {
                        if ($transcript->{$selected_acc}->{'indel'} eq "") {
                            if ($transcript->{$tr}->{'indel'} ne "") {
                                $selected_acc = $tr;    #INDEL takes precedence
                            } else {
                                if ($transcript->{$tr}->{'snp'} ne "" and $transcript->{$selected_acc}->{'snp'} eq "") {
                                    $selected_acc = $tr;    #INDEL takes precedence
                                }
                            }
                        }
                    }
                }
            }
        }
        
        $non_red_trans->{$selected_acc} = $transcript->{$selected_acc};
        
        foreach my $tr (keys %{$transcript_seq->{$seq}}) {
            next if ($selected_acc eq $tr); #Not necessary to throw this info away
            #next if ($transcript->{$selected_acc}->{'gene'} eq $transcript->{$tr}->{'gene'});
            push @{$non_red_trans->{$selected_acc}->{'others'}}, $tr;
        }
        
    }
    
    return $non_red_trans;
    
}




##------ GET TRANSCRIPTS -------##
sub get_transcripts_from_resultdb {

	# sub routine to extract Ribo-seq information from SQlite result database

    my ($dbh,$tbl,$tis,$TIS_calling_method,$ens_db) = @_;

	my $var_tracker = {};
	my $transcript = {};
	my $gene_transcript = {};

	my ($transcript2geneid,$annotated_tr) = transcript_gene_id($dbh,$tbl);
		
	print STDOUT "Extracting transcripts form SQLite database. Please wait ....\n";
	my $query = "SELECT DISTINCT tr_stable_id, chr, start, start_codon, dist_to_aTIS, aTIS_call, annotation, peak_shift, SNP, INDEL, aa_seq FROM ".$tbl.";";
 	my $sth = $dbh->prepare($query);
	$sth->execute();
	
	while ( my @row = $sth->fetchrow()) {
        
        #Change undefined values (because of empty SQLite values) to empty strings
        foreach (@row) {$_ = '' unless defined};
        my ($tr_stable_id, $chr, $start, $start_codon, $dist_to_aTIS, $aTIS_call, $annotation, $peak_shift, $snp, $indel, $aa_seq) = @row;
        
        
		# if instructed to not keep aTIS with not enough coverage to call TIS
		if (uc($tis_call) eq "N") {next if (uc($aTIS_call) eq 'NO_DATA' or uc($aTIS_call) eq 'FALSE')}
		$aa_seq =~ s/\*//g;
		next if (length($aa_seq) < $mslength);	# skip if sequence is less than minimum allowed amino acid length

		# Skip all non aTIS without SNP or indel information that corresponds to an annotated TIS i.e redundant non annotated TIS
		my $red_tis = 0;
		if ($annotation ne 'aTIS' and $snp eq "" and $indel eq "") {
			if ($transcript2geneid->{$tr_stable_id}) {
				my $gene = $transcript2geneid->{$tr_stable_id}->{'gene'};
				foreach my $start1 (keys %{$annotated_tr->{$gene}}) {
					if ($start1 == $start) {$red_tis = 1}
				}
			}
		}
		next if ($red_tis == 1);
        
        #Check for ORFs present in the ORF table but also in transcripts not in the translated transcripts table
        #PRICE and SPECtre call ORFs based on the gtf transcript annotation, so they can call ORFs in uncalled transcripts
        #These transcripts are not present in the tr_translation table and thus, not present in the transcript2geneid hash
        #For these transcripts, include gene id and biotype info from ensembl
        
        #if stable id not present in transcript2geneid:
        #   search in ensembl for gene id and biotype
        unless(defined($transcript2geneid->{$tr_stable_id})){
            ($transcript2geneid->{$tr_stable_id}->{'gene'}, $transcript2geneid->{$tr_stable_id}->{'biotype'}) = search_gene_ensembl($tr_stable_id, $ens_db);
        }

		# create unique transcript ID
		my $tr = $tr_stable_id."_".$chr."_".$start."_".$annotation;
		if ($snp ne "" or $indel ne "") {		# if snp or indel info exist for current record
			if ($var_tracker->{$tr}) {
				$var_tracker->{$tr}++;
				$tr = $tr."_".$var_tracker->{$tr};
			} else {
				$var_tracker->{$tr} = 1;
				$tr = $tr."_".$var_tracker->{$tr};
			}
		}

		$transcript->{$tr}->{'chr'} 		= $chr;
		$transcript->{$tr}->{'snp'} 		= $snp;
        $transcript->{$tr}->{'indel'}       = $indel;
		$transcript->{$tr}->{'codon'} 		= $start_codon;
		$transcript->{$tr}->{'aTIS_call'} 	= $aTIS_call;
		$transcript->{$tr}->{'anno'} 		= $annotation;
		$transcript->{$tr}->{'peak_shift'} 	= $peak_shift;
		$transcript->{$tr}->{'seq'} 		= $aa_seq;
		$transcript->{$tr}->{'biotype'} 	= $transcript2geneid->{$tr_stable_id}->{'biotype'};
		$transcript->{$tr}->{'gene'} 		= $transcript2geneid->{$tr_stable_id}->{'gene'};
		$transcript->{$tr}->{'tr'} 			= $tr_stable_id;
		$transcript->{$tr}->{'red'} 		= "N";			# Set all transcript to default redundant value of N
		
		push @{$gene_transcript->{$transcript2geneid->{$tr_stable_id}->{'gene'}}}, $tr;
	
	}
	$sth->finish();
	
	return $transcript,$gene_transcript;

}

#PRICE and SPECtre call ORFs based on the gtf transcript annotation, so they can call ORFs in uncalled transcripts
#These transcripts are not present in the tr_translation table and thus, not present in the transcript2geneid hash
#For these transcripts, include gene id and biotype info from ensembl
sub search_gene_ensembl {
    
    #Catch
    my $tr_stable_id = $_[0];
    my $ens_db = $_[1];
    
    #Init
    my $gene_id;
    my $biotype;
    
    #Connect to ens_db
    my $user = "";
    my $pw = "";
    
    # Connect to ensembl sqlite database
    my $dbh  = DBI->connect('DBI:SQLite:'.$ens_db,$user,$pw,
    { RaiseError => 1},) || die "Database connection not made: $DBI::errstr";
    
    my $query = "SELECT g.stable_id, g.biotype FROM gene AS g JOIN transcript AS t ON g.gene_id=t.gene_id WHERE t.stable_id='".$tr_stable_id."';";
    my $execute = $dbh->prepare($query);
    $execute->execute();
    while(my @result = $execute->fetchrow()){
        $gene_id = $result[0];
        $biotype = $result[1];
    }
    
    $execute->finish();
    
    #Disconnect
    $dbh->disconnect();
    
    return ($gene_id, $biotype);
}

sub convert_species {
    
    #Catch
    my $species = $_[0];
    my $tmp_folder = $_[1];
    
    #Conversion
    my $species_peff = (uc($species) eq "MOUSE") ? "Mus musculus"
    : (uc($species) eq "HUMAN") ? "Homo sapiens"
    : (uc($species) eq "ARABIDOPSIS") ? "Arabidopsis thaliana"
    : (uc($species) eq "FRUITFLY") ? "Drosophila melanogaster" : "";
    #Error check
    if ($species_peff eq ""){
        print "Error could not convert the species to peff format: ".$species."\n";
        die;
    }
    
    #Search tax ID
    #download tax ID lookup file
    if (!-e $tmp_folder."/names.dmp"){
        system("wget -q ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz");
        system("mv taxdump.tar.gz ".$tmp_folder);
        #Unpack
        system("gunzip ".$tmp_folder."/taxdump.tar.gz");
        system("tar -xf ".$tmp_folder."/taxdump.tar -C ".$tmp_folder);
    }
    
    #Read tax file
    my $taxid = "";
    my $file = $tmp_folder."/names.dmp";
    open(my $FR, '<', $file) or die "Could not find $file file";
    while (my $line = <$FR>){
        chomp $line;
        $line =~ m/^(.*?)\t\|\t(.*?)\t/;
        if (uc($2) eq uc($species_peff)){
            $taxid = $1;
            last;
        }
    }
    
    return ($species_peff, $taxid);
}

sub transcript_gene_id {

	my $dbh = $_[0];
	my $tbl = $_[1];

	my $annotated_tr = {};
	my $transcript2geneid = {};

	my $query = "SELECT a.stable_id, a.biotype, a.gene_stable_id, b.start, b.annotation, b.aTIS_call FROM tr_translation a, $tbl b WHERE a.stable_id == b.tr_stable_id";
 	my $sth = $dbh->prepare($query);
	$sth->execute();
	while ( my ($stable_id, $biotype, $gene_stable_id, $start, $annotation, $aTIS_call) = $sth->fetchrow()) {
        $transcript2geneid->{$stable_id}->{'gene'} = $gene_stable_id;
        $transcript2geneid->{$stable_id}->{'biotype'} = $gene_stable_id;
        
        if (uc($tis_call) eq "N") {next if (uc($aTIS_call) eq 'NO_DATA' or uc($aTIS_call) eq 'FALSE')}
        if ($annotation eq 'aTIS') {$annotated_tr->{$gene_stable_id}->{$start} = 1;}
	}
	$sth->finish();

	return $transcript2geneid, $annotated_tr;
}

sub transcript_gene_id_per_chr {
    
    my $dbh = $_[0];
    my $tbl = $_[1];
    my $chr = $_[2];
    my $ens_db = $_[3];
    
    my $annotated_tr = {};
    my $transcript2geneid = {};
    
    #Attach ens db to dbh
    $dbh->do("ATTACH '".$ens_db."' AS ens;");
    
    my $query = "SELECT a.stable_id, a.biotype, a.gene_stable_id, b.start, b.annotation, b.aTIS_call, g.description FROM tr_translation AS a JOIN $tbl AS b ON a.stable_id == b.tr_stable_id JOIN ens.gene AS g ON a.gene_stable_id=g.stable_id WHERE a.chr='".$chr."';";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    while ( my ($stable_id, $biotype, $gene_stable_id, $start, $annotation, $aTIS_call, $description) = $sth->fetchrow()) {
        #Capture only essential info of description
        if($description =~ m/(.*?)\[.*\]/){
            $description=$1;
        }
        #Parse
        unless (exists $transcript2geneid->{$stable_id}){
            $transcript2geneid->{$stable_id}->{'gene'} = $gene_stable_id;
            $transcript2geneid->{$stable_id}->{'biotype'} = $biotype;
            $transcript2geneid->{$stable_id}->{'description'} = $description;
        }
            
        if (uc($tis_call) eq "N") {next if (uc($aTIS_call) eq 'NO_DATA' or uc($aTIS_call) eq 'FALSE')}
        if ($annotation eq 'aTIS') {$annotated_tr->{$gene_stable_id}->{$start} = 1;}
    }
    $sth->finish();

    return $transcript2geneid, $annotated_tr;
}


##------ MAPPING TO CANONICAL -------##
sub remote_biomart_mapping {
	
	# This function remotely searches Biomart website for ensemble transcript to swissprot ID mapppings
	# input is the ensemble database name
	
	my $mapping_db = $_[0];
	my $version = $_[1];
	my $ref_attribute = $_[2];
	my $other_attribute = $_[3];
	
	print STDOUT "Retriving the list of Biomart Ensembl transcript to Swissprot ids mappings. Please wait...\n";
	
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

sub get_analysis_ids{
    
    # Catch
    my $dbh    = $_[0];
    my $ids_in = $_[1]; #Either comma separated list of identifiers or "all"
    my @idsref;
    
    if ($ids_in eq "all") {
        
        my $query = "select ID, SNP from TIS_overview";
        my $sth = $dbh->prepare($query);
        $sth->execute();
        
        while ( my ($id, $snp) = $sth->fetchrow()) {
            my $add = "";
            if (uc($snp) eq "NO") {
                push @idsref, $id;
            } else {
                    $add = "_".$snp;
                push @idsref, $id.$add;
            }
        }
        $sth->finish();
    } else {
        
        my @sel_ids  = split ',', $ids_in;
        foreach (@sel_ids) {
            
            my $query = "select ID, SNP from TIS_overview where ID = $_";
            my $sth = $dbh->prepare($query);
            $sth->execute();
            my ($id, $snp) = $sth->fetchrow();
            
            my $add = "";
            if (uc($snp) eq "NO") {
                push @idsref, $id;
            } else {
                $add = "_".$snp;
                push @idsref, $id.$add;
            }
            $sth->finish();
        }
    }
    
    return @idsref;
    
}

sub get_analysis_ids_indels {
    
	# Catch
	my $dbh    = $_[0];
	my $ids_in = $_[1]; #Either comma separated list of identifiers or "all"
	my @idsref;

	if ($ids_in eq "all") {

		my $query = "select ID, SNP, indel from TIS_overview";
	    my $sth = $dbh->prepare($query);
		$sth->execute();

		while ( my ($id, $snp, $indel) = $sth->fetchrow()) {
            my $add = "";
			if ((uc($snp) eq "NO") && (uc($indel) eq "NO")) {
				push @idsref, $id;
			} elsif (uc($snp) eq "NO") {
                my @indel_tools = split '_', $indel;
                foreach my $tool (@indel_tools){
                    $add = $add."_indel".$tool;
                }
				push @idsref, $id.$add;
            } elsif (uc($indel) eq "NO") {
                my @snp_tools = split '_', $snp;
                foreach my $tool (@snp_tools){
                    $add = $add."_snp".$tool;
                }
                push @idsref, $id.$add;
            } else {
                my @snp_tools = split '_', $snp;
                foreach my $tool (@snp_tools){
                    $add = $add."_snp".$tool;
                }
                my @indel_tools = split '_', $indel;
                foreach my $tool (@indel_tools){
                    $add = $add."_indel".$tool;
                }
                push @idsref, $id.$add;
            }
		}
		$sth->finish();
    } else {
        	
		my @sel_ids  = split ',', $ids_in;
		foreach (@sel_ids) {

			my $query = "select ID, SNP, indel from TIS_overview where ID = $_";
		    my $sth = $dbh->prepare($query);
			$sth->execute();
			my ($id, $snp, $indel) = $sth->fetchrow();
			
            my $add = "";
            if ((uc($snp) eq "NO") && (uc($indel) eq "NO")) {
                push @idsref, $id;
            } elsif (uc($snp) eq "NO") {
                my @indel_tools = split '_', $indel;
                foreach my $tool (@indel_tools){
                    $add = $add."_indel".$tool;
                }
                push @idsref, $id.$add;
            } elsif (uc($indel) eq "NO") {
                my @snp_tools = split '_', $snp;
                foreach my $tool (@snp_tools){
                    $add = $add."_snp".$tool;
                }
                push @idsref, $id.$add;
            } else {
                my @snp_tools = split '_', $snp;
                foreach my $tool (@snp_tools){
                    $add = $add."_snp".$tool;
                }
                my @indel_tools = split '_', $indel;
                foreach my $tool (@indel_tools){
                    $add = $add."_indel".$tool;
                }
                push @idsref, $id.$add;
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

sub get_arguments{
    
    # Catch
    my $dbh = $_[0];
    
    # Get input variables
    my $query = "select value from `arguments` where variable = \'ensembl_version\'";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    my $ensemblversion = $sth->fetch()->[0];
    
    $query = "select value from `arguments` where variable = \'species\'";
    $sth = $dbh->prepare($query);
    $sth->execute();
    my $species = $sth->fetch()->[0];
    
    $query = "select value from `arguments` where variable = \'igenomes_root\'";
    $sth = $dbh->prepare($query);
    $sth->execute();
    my $igenomes_root = $sth->fetch()->[0];
    
    $query = "select value from `arguments` where variable = \'nr_of_cores\'";
    $sth = $dbh->prepare($query);
    $sth->execute();
    my $nr_of_cores = $sth->fetch()->[0];
    
    $query = "select value from `arguments` where variable = \'ens_db\'";
    $sth = $dbh->prepare($query);
    $sth->execute();
    my $ens_db = $sth->fetch()->[0];
    
    
    # Return input variables
    return($ensemblversion,$ens_db,$species,$igenomes_root,$nr_of_cores);
    
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
		print STDOUT "Error creating file ", $file, ". Ensure you have the required permission.\n";
	} else {
		print STDOUT "Error reading file ", $file, ". Ensure the file exist.\n";
	}
	
}

sub timer {
	my $startRun = shift;

	my $endRun 	= time();
	my $runTime = $endRun - $startRun;

	printf("\nTotal running time: %02d:%02d:%02d\n\n", int($runTime / 3600), int(($runTime  % 3600) / 60), int($runTime % 60));
}
