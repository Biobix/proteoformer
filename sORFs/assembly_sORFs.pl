#!/usr/bin/perl -w

use strict;

use warnings;

use DBI;

use DBD::SQLite;

use Data::Dumper;

use Getopt::Long;

use v5.10;

use Parallel::ForkManager;

use Storable 'dclone';

use Cwd;



##############

##Command-line

# ./assembly_sORFs.pl  --sqlite_db SQLite/results2.db --tis_ids 1 --min_length 10 --max_length 100 --ex_overlap 0.5



# ./assembly_sORFs.pl  --sqlite_db SQLite/results_LTM.db --tis_ids 3 --min_length 7 --max_length 150 --coverage 0.25 --ex_overlap 0.5





# get the command line arguments

my ($tmpfolder,$snp,$min,$max,$coverage,$sqlite_db,$work_dir,$tis_ids,$out_sqlite,$R,$ex_overlap);



GetOptions(

"sqlite_db=s"=>\$sqlite_db,                 # The sqlite DB holding all RIBO-pipeline results,                     mandatory argument

"dir:s"=>\$work_dir,                        # Path to the working directory,                                        optional  argument

"tmp:s" =>\$tmpfolder,                      # Folder where temporary files are stored,                              optional  argument (default = $TMP env setting

"min_length=i"=>\$min,                      # Minimal sORF length in Amino Acids                                    mandatory argument

"max_length=i"=>\$max,                      # Maximal sORF length in Amino Acids                                    mandatory argument

#"coverage:f" => \$coverage,                 # The minimal coverage for sORF sequence based on total ribo_reads      optional argument (default .75)

"ex_overlap:f" => \$ex_overlap,             # The maximal exon overlap for a sORF sequence                          optional argument (default .75)

"snp:s"=>\$snp,                             # The snp calling algorithm applied                                     optional argument (default "NO", others can be "samtools", "GATK")

"tis_ids=s"  =>\$tis_ids,                   # list of analysis ids                                                  mandatory argument

"out_sqlite=s"=>\$out_sqlite,               # Galaxy specific history file location                                                         Galaxy specific

);



my $CWD             = getcwd;

my $HOME            = $ENV{'HOME'};

my $TMP             = ($ENV{'TMP'}) ? $ENV{'TMP'} : ($tmpfolder) ? $tmpfolder : "$CWD/tmp"; # (1) get the TMP environment variable, (2) get the $tmpfolder variable, (3) get current_working_dir/tmp

print "The following tmpfolder is used                          : $TMP\n";

print "The following results db folder is used                  : $sqlite_db\n";



#Check if tmpfolder exists, if not create it...

if (!-d "$TMP") {

    system ("mkdir ". $TMP);

}



#comment on these

if ($work_dir){

    print "Working directory                                        : $work_dir\n";

} else {

    #Choose default value

    $work_dir = $CWD;

    print "Working directory                                        : $CWD\n";

}

if ($min){

    print "Minimal sORF length (AA)                                 : $min\n";

} else {

    die "\nDon't forget to pass the minimal sORF length using the --min_length or -min argument!\n\n";

}

if ($max){

    print "Maximal sORF length (AA)                                 : $max\n";

} else {

    die "\nDon't forget to pass the maximal sORF length using the --max_length or -max argument!\n\n";

}

#if ($coverage){

#    print "The minimal sORF coverage based on ribosome profiles     : $coverage\n";

#} else {

#    #Choose default value for coverage

#    $coverage = 0.75;

#}

if ($ex_overlap){

    print "The maximal sORF exon overlap based on ribosome profiles : $ex_overlap\n";

} else {

    #Choose default value for exon overlap

    $ex_overlap = 0.5;

}

if ($snp){

    print "The snp calling algorithm used is                        : $snp\n";

} else {

    #Choose default value for mapper

    $snp = "NO";

    print "The snp calling algorithm used is                        : $snp\n";

}



# Create output files for command line script

if (!defined($out_sqlite))       {$out_sqlite          = $sqlite_db;}



# DB settings

# Sqlite Riboseq

my $db_results  = $sqlite_db;

my $dsn_results = "DBI:SQLite:dbname=$db_results";

my $us_results  = "";

my $pw_results  = "";



#Get the input variables (Check nog welke eruit mogen)

#Get arguments from arguments table

my ($ensemblversion,$species,$ens_db,$IGENOMES_ROOT,$cores,$Mreads,$mean_length_fastq1,$mean_length_fastq2) = get_arguments($dsn_results,$us_results,$pw_results);



# Sqlite Ensembl

my $db_ENS  = $ens_db;

my $dsn_ENS = "DBI:SQLite:dbname=$db_ENS";

my $us_ENS  = "";

my $pw_ENS  = "";



print "The igenomes_root folder used is                         : $IGENOMES_ROOT\n";

print "Number of cores to use for Mapping                       : $cores\n";

print "The following Ensembl db folder is used                  : $ens_db\n";

print "Total number of CHX reads in millions of reads           : $Mreads\n";



#Conversion for species terminology

my $spec = ($species eq "mouse") ? "Mus_musculus" : ($species eq "human") ? "Homo_sapiens" : ($species eq "arabidopsis") ? "Arabidopsis_thaliana" : ($species eq "fruitfly") ? "Drosophila_melanogaster" : "";

my $spec_short = ($species eq "mouse") ? "mmu" : ($species eq "human") ? "hsa" : ($species eq "arabidopsis") ? "ath" : ($species eq "fruitfly") ? "dme" : "";

#Old mouse assembly = NCBIM37, new one is GRCm38

my $assembly = ($species eq "mouse" && $ensemblversion >= 70 ) ? "GRCm38"

: ($species eq "mouse" && $ensemblversion < 70 ) ? "NCBIM37"

: ($species eq "human") ? "GRCh37"

: ($species eq "arabidopsis") ? "TAIR10"

: ($species eq "fruitfly") ? "BDGP5" : "";



#my $chrom_file = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Annotation/Genes/ChromInfo.txt";

my $BIN_chrom_dir = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Sequence/Chromosomes_BIN";



#################################################

##                                             ##

## ASSEMBLY OF INTERGENIC TRANSLATION PRODUCTS ##

##                                             ##

#################################################



# Start time

my $start = time;

## Get chromosome sizes

print "Getting chromosome sizes  ...\n";

my $chr_sizes = get_chr_sizes();

# Get the analysis_id that corresponds to the TIS-calling input parameters

my $idsref = get_analysis_ids($dsn_results,$us_results,$pw_results,$tis_ids);



print "\nGet chromosomes... \n";

## Get chromosomes based on seq_region_id ##

my $chrs = get_chrs($dsn_ENS,$us_ENS,$pw_ENS,get_chr_sizes(),$assembly);



# Create binary chromosomes if they don't exist

print "Checking/Creating binary chrom files ...\n";

if (!-d "$BIN_chrom_dir") {

    create_BIN_chromosomes($BIN_chrom_dir,$cores,$chrs,$TMP);

}



my $analysis_id;

#Loop over all selected analysis_ids

print "Assembly of translation products with RPKM and coverage analysis ...\n";

foreach $analysis_id (@$idsref) {

    

    print "Processing analysis_id $analysis_id ...\n";

    construct_sORFs($chrs,$TMP,$analysis_id,$min,$max,$Mreads);



    ## Store in DB

    print "   Storing all sORFs in DB \n";

    store_in_db($dsn_results,$us_results,$pw_results,$analysis_id,$work_dir,$chrs,$TMP);

}



#Move to galaxy history

system ("mv ".$sqlite_db." ".$out_sqlite);



# End time

print "   DONE! \n";

my $end = time - $start;

printf("runtime assembly: %02d:%02d:%02d\n\n",int($end/3600), int(($end % 3600)/60), int($end % 60));



############

# THE SUBS #

############



### Construct sORFs ###



sub construct_sORFs{

    

    #Catch

    my $chrs        =   $_[0];

    my $tmp         =   $_[1];

    my $analysis_id =   $_[2];

    my $min_length  =   $_[3];

    my $max_length  =   $_[4];

    my $Mreads      =   $_[5];

    

    # Get R from TIS_overview

    my $R = get_R($dsn_results,$us_results,$pw_results,$analysis_id);

    

    # Init multi core

    my $pm = new Parallel::ForkManager($cores);

    

    ## Loop over all chromosomes

    foreach my $chr (sort keys %{$chrs}){

        

        ### Start parallel process

        $pm->start and next;

        

        ### DBH per process

        my $dbh = dbh($dsn_results,$us_results,$pw_results);

        

        ### Output fasta and db_csv file per process

        open TMP_db, "+>>".$TMP."/".$analysis_id."_".$chr."_sORF_tmp.csv" or die $!;

        

        #Get CHX and LTM reads

        my ($CHX_for,$LTM_for,$CHX_rev,$LTM_rev) = get_reads($chr);

        

        ## Get intergenic_ids per chr from transcript calling

        my $ig_sORF_starts = get_intergenic_sORF_starts_per_chromosome($dbh,$analysis_id,$chr);



        # Init

        my ($sORF_seq,$tr_sORF_seq,$tmp_sORF_seq,$ig_sORF_id,$g_sORF_id,$TIS,$strand,$start_codon,$peak_shift,$count,$sORF_length,$FPKM,$coverage_sORF,$sorf_end,$R_sORF,$sorf_chr,$sorf_begin,$sorf_strand,$annotation,$biotype,$upstream_gene_distance,$downstream_gene_distance,$coverage_uniformity,$exon_overlap,$start_seq_region_start,$start_seq_region_end,$end_seq_region_start,$end_seq_region_end,$phase,$end_phase,$in_frame,$exon_id,$peptide_mass);

        my $ig_sorfs = {};

        my $g_sorfs = {};

        my $seq_region_id = $chrs->{$chr}{'seq_region_id'};

        my $max_DNA_length = ($max_length*3) -1;

        

        foreach my $ig_sORF_start ( keys %{$ig_sORF_starts}){

            #my %tr_SNPs;

            $sORF_seq = '';

            $tmp_sORF_seq = '';

            

            # Split transcript and TIS

            $ig_sORF_id                     =   $ig_sORF_starts->{$ig_sORF_start}{'intergene_id'};

            $TIS                            =   $ig_sORF_starts->{$ig_sORF_start}{'start'};

            $strand                         =   $ig_sORF_starts->{$ig_sORF_start}{'strand'};

            $start_codon                    =   $ig_sORF_starts->{$ig_sORF_start}{'start_codon'};

            $peak_shift                     =   $ig_sORF_starts->{$ig_sORF_start}{'peak_shift'};;

            $count                          =   $ig_sORF_starts->{$ig_sORF_start}{'count'};

            $downstream_gene_distance       =   $ig_sORF_starts->{$ig_sORF_start}{'downstream_gene_distance'};

            $upstream_gene_distance         =   $ig_sORF_starts->{$ig_sORF_start}{'upstream_gene_distance'};

            

            # Create max_length sequence from TIS (strand specific)

            $tmp_sORF_seq = ($strand eq '1') ? get_sequence($chr,$TIS,$TIS + $max_DNA_length) : revdnacomp(get_sequence($chr,$TIS - $max_DNA_length,$TIS));

            

            # Sequence ad the end of chromosome

            if ((length($tmp_sORF_seq) % 3) != 0){next;}

            

            # Create translated sequence

            ($tr_sORF_seq,$sORF_seq,$sORF_length,$peptide_mass) = translate($tmp_sORF_seq,$min,$max);

            

            #next if no sorf sequence

            if($sORF_seq eq 'X'){next;}

            

            #Save in ig_sorfs

            $sorf_end   =   ($strand eq '1') ? $TIS + (($sORF_length * 3) -1) : $TIS;

            $sorf_begin =   ($strand eq '1') ? $TIS : $TIS - (($sORF_length * 3) -1);

            

            $ig_sorfs->{$ig_sORF_start}{'intergene_id'}             =   $ig_sORF_id;

            $ig_sorfs->{$ig_sORF_start}{'sorf_chr'}                 =   $chr;

            $ig_sorfs->{$ig_sORF_start}{'sorf_begin'}               =   $sorf_begin;

            $ig_sorfs->{$ig_sORF_start}{'FPKM_begin'}               =   $sorf_begin - 15;

            $ig_sorfs->{$ig_sORF_start}{'sorf_end'}                 =   $sorf_end;

            $ig_sorfs->{$ig_sORF_start}{'FPKM_end'}                 =   $sorf_end +15;

            $ig_sorfs->{$ig_sORF_start}{'sorf_strand'}              =   $strand;

            $ig_sorfs->{$ig_sORF_start}{'start_codon'}              =   $start_codon;

            $ig_sorfs->{$ig_sORF_start}{'peak_shift'}               =   $peak_shift;

            $ig_sorfs->{$ig_sORF_start}{'count'}                    =   $count;

            $ig_sorfs->{$ig_sORF_start}{'sorf_length'}              =   $sORF_length;

            $ig_sorfs->{$ig_sORF_start}{'sorf_seq'}                 =   $sORF_seq;

            $ig_sorfs->{$ig_sORF_start}{'tr_sorf_seq'}              =   $tr_sORF_seq;

            $ig_sorfs->{$ig_sORF_start}{'downstream_gene_distance'} =   $downstream_gene_distance;

            $ig_sorfs->{$ig_sORF_start}{'upstream_gene_distance'}   =   $upstream_gene_distance;

            $ig_sorfs->{$ig_sORF_start}{'annotation'}               =   'intergenic';

            $ig_sorfs->{$ig_sORF_start}{'mass'}                     =   $peptide_mass;

            

            #For LTM reads start from sorf_begin -1 to calculate total LTM

            #Otherwize, if sORF has only LTM_reads on sorf_begin -1 position total_LTM = 0 -> illegal division

            $ig_sorfs->{$ig_sORF_start}{'LTM_begin'}   =    ($strand eq '1') ?  $sorf_begin - 1                          :   $ig_sorfs->{$ig_sORF_start}{'sorf_begin'};

            $ig_sorfs->{$ig_sORF_start}{'LTM_end'}     =    ($strand eq '1') ?  $ig_sorfs->{$ig_sORF_start}{'sorf_end'}  :   $sorf_end + 1;

        }

        

        # Match reads to intergenic sORFs

        $ig_sorfs = match_reads_to_intergenic_sORFs($ig_sorfs,$CHX_for,$CHX_rev,$LTM_for,$LTM_rev);

        

        # Loop over sORFs in $ig_sORFs

        foreach my $ig_sorf_id (keys %{$ig_sorfs}){

            

            # Get sORF specific reads

            my ($LTM_reads,$CHX_reads) = get_overlapping_reads($ig_sorf_id,$ig_sorfs,$chr,$CHX_for,$CHX_rev,$LTM_for,$LTM_rev);

            

            # Calculate FPKM, coverage and R

            ($FPKM,$coverage_sORF,$R_sORF,$coverage_uniformity) = calculate_FPKM_and_coverage($ig_sorfs,$ig_sorf_id,$Mreads,$LTM_reads,$CHX_reads,$mean_length_fastq1,$mean_length_fastq2);

            

            #Check if R value and coverage are higher or equal as $R and

            if($R_sORF < $R){next;}

            if($coverage_sORF <= 0){next;}

            

            # Save sORF in csv

            $ig_sORF_id                 =   $ig_sorfs->{$ig_sorf_id}{'intergene_id'};

            $sorf_chr                   =   $ig_sorfs->{$ig_sorf_id}{'sorf_chr'};

            $sorf_begin                 =   $ig_sorfs->{$ig_sorf_id}{'sorf_begin'};

            $sorf_end                   =   $ig_sorfs->{$ig_sorf_id}{'sorf_end'};

            $sorf_strand                =   $ig_sorfs->{$ig_sorf_id}{'sorf_strand'};

            $start_codon                =   $ig_sorfs->{$ig_sorf_id}{'start_codon'};

            $peak_shift                 =   $ig_sorfs->{$ig_sorf_id}{'peak_shift'};

            $count                      =   $ig_sorfs->{$ig_sorf_id}{'count'};

            $sORF_length                =   $ig_sorfs->{$ig_sorf_id}{'sorf_length'};

            $sORF_seq                   =   $ig_sorfs->{$ig_sorf_id}{'sorf_seq'};

            $tr_sORF_seq                =   $ig_sorfs->{$ig_sorf_id}{'tr_sorf_seq'};

            $downstream_gene_distance   =   $ig_sorfs->{$ig_sorf_id}{'downstream_gene_distance'};

            $upstream_gene_distance     =   $ig_sorfs->{$ig_sorf_id}{'upstream_gene_distance'};

            $annotation                 =   $ig_sorfs->{$ig_sorf_id}{'annotation'};

            $peptide_mass               =   $ig_sorfs->{$ig_sorf_id}{'mass'};

            

            # Check how to implement

            my $out_cnt        =   '1';

            # Should becom $_->{'SNP'}

            my $snp_tmp = '0';

            

            #Add extra columns to be compatible with genic_sORF_calling

            my $biotype         =   'NA';

            my $exon_overlap    =   0;

            my $in_frame        =   'NA';

            

            print TMP_db $ig_sORF_id.",".$biotype.",".$sorf_chr.",".$sorf_strand.",".$sorf_begin.",".$sorf_end.",".$sORF_length.",".$peptide_mass.",".$start_codon.",".$downstream_gene_distance.",".$upstream_gene_distance.",".$annotation.",".$exon_overlap.",".$in_frame.",".$peak_shift.",".$count.",".$R_sORF.",".$coverage_sORF.",".$coverage_uniformity.",".$FPKM.",".$snp_tmp.",".$sORF_seq.",".$tr_sORF_seq."\n";

            

        }

        

        ## Get gene_ids per chr from TIScalling

        my $genes = get_genes_per_chromosome($dbh,$analysis_id,$chr);

            

        #Run over genes

        foreach my $gene_id (keys %{$genes}){

        

            #Init

            my $g_sorfs_tmp = {};

            

            #get data

            ## Get TISses per gene_id

            my $g_sORF_starts = get_genic_sORF_starts_per_gene_id($dbh,$analysis_id,$gene_id);



            #Check if gene protein_coding

            if($genes->{$gene_id}{'biotype'} eq 'protein_coding') {

            

                ## Get transcripts from ENS db

                my $gene_id_transcripts = get_gene_id_transcripts($dsn_ENS,$us_ENS,$pw_ENS,$gene_id);

                

                #Create list of aTISses, exons and UTRs

                #Init

                my $aTIS_hash   =   {};

                my $exon_hash   =   {};

                my $UTR5_hash   =   {};

                my $UTR3_hash   =   {};

                my ($start,$stop,$strand,$gene_id_transcript_id,$exon);

                

                foreach $gene_id_transcript_id (keys %{$gene_id_transcripts}){

                

                    #Get translation data

                    my ($trans,$exons)  =   get_translation_data($dsn_ENS,$us_ENS,$pw_ENS,$gene_id_transcript_id,$gene_id_transcripts,$seq_region_id,$chr);

                    my $strand          =   $gene_id_transcripts->{$gene_id_transcript_id}->{'seq_region_strand'};

                

                    #Check if transcript has translation info

                    if($trans->[0]){

                        

                        #Get genomic aTIS from first exon and save in aTIS_array

                        $start = ($strand eq '1') ? $exons->{$trans->[0]}{'seq_region_start'} + ${$trans}[2] -1 : $exons->{$trans->[0]}{'seq_region_end'} - ${$trans}[2] + 1;

                        $stop =  ($strand eq '1') ? $exons->{$trans->[1]}{'seq_region_start'} + ${$trans}[3] -1 : $exons->{$trans->[1]}{'seq_region_end'} - ${$trans}[3] + 1;

                        $aTIS_hash->{$start}    =   $start;

                        

                        #Save 5UTR_start_exon and 3UTR_stop_exon in hash

                        $UTR5_hash->{$gene_id_transcript_id}->{$trans->[0]}->{'start'}         =   ($strand eq '1') ? $exons->{$trans->[0]}{'seq_region_start'} : $start + 1;

                        $UTR5_hash->{$gene_id_transcript_id}->{$trans->[0]}->{'end'}           =   ($strand eq '1') ? $start - 1 : $exons->{$trans->[0]}{'seq_region_end'};

                        $UTR5_hash->{$gene_id_transcript_id}->{$trans->[0]}->{'transcript_id'} =   $gene_id_transcript_id;

                        

                        $UTR3_hash->{$gene_id_transcript_id}->{$trans->[0]}->{'start'}         =   ($strand eq '1') ? $stop + 1 : $exons->{$trans->[1]}{'seq_region_start'};

                        $UTR3_hash->{$gene_id_transcript_id}->{$trans->[0]}->{'end'}           =   ($strand eq '1') ? $exons->{$trans->[1]}{'seq_region_end'} : $stop - 1;

                        $UTR3_hash->{$gene_id_transcript_id}->{$trans->[0]}->{'transcript_id'} =   $gene_id_transcript_id;

                        

                        #Delete UTRs from translation start and stop exon (work with new_seq_region_start for single exon transcripts)

                        $start_seq_region_start =   $exons->{$trans->[0]}{'seq_region_start'};

                        $start_seq_region_end   =   $exons->{$trans->[0]}{'seq_region_end'};

                        $end_seq_region_start   =   $exons->{$trans->[1]}{'seq_region_start'};

                        $end_seq_region_end     =   $exons->{$trans->[1]}{'seq_region_end'};

                        

                        $exons->{$trans->[0]}{'seq_region_start'}   =   ($strand eq '1') ? $start_seq_region_start + ${$trans}[2] -1 : $start_seq_region_start;

                        $exons->{$trans->[0]}{'seq_region_end'}     =   ($strand eq '1') ?  $start_seq_region_end : $start_seq_region_end - ${$trans}[2] + 1;

                        

                        if($trans->[0] eq $trans->[1]){

                            $exons->{$trans->[1]}{'seq_region_start'}   =   ($strand eq '1') ? $exons->{$trans->[0]}{'seq_region_start'} : $end_seq_region_end - ${$trans}[3] +1;

                            $exons->{$trans->[1]}{'seq_region_end'}     =   ($strand eq '1') ? $end_seq_region_start + ${$trans}[3] -1 : $exons->{$trans->[0]}{'seq_region_end'};

                        }else{

                            $exons->{$trans->[1]}{'seq_region_start'}   =   ($strand eq '1') ? $end_seq_region_start : $end_seq_region_end - ${$trans}[3] +1;

                            $exons->{$trans->[1]}{'seq_region_end'}     =   ($strand eq '1') ? $end_seq_region_start + ${$trans}[3] -1 : $end_seq_region_end;

                        }

                        

                        #Add exons to exon_hash | but not UTR exons (rank < as rank of $trans->[0] exon)

                        foreach $exon (keys %{$exons}){

                            

                            # but not UTR exons (rank < as rank of $trans->[0] exon)

                            if($exons->{$exon}{'rank'} < $exons->{$trans->[0]}{'rank'}){

                                $UTR5_hash->{$gene_id_transcript_id}->{$exon}{'start'} = $exons->{$exon}{'seq_region_start'};

                                $UTR5_hash->{$gene_id_transcript_id}->{$exon}{'end'} = $exons->{$exon}{'seq_region_end'};

                                $UTR5_hash->{$gene_id_transcript_id}->{$exon}{'transcript_id'} = $gene_id_transcript_id;

                                next;

                            }elsif($exons->{$exon}{'rank'} > $exons->{$trans->[1]}{'rank'}){

                                $UTR3_hash->{$gene_id_transcript_id}->{$exon}{'start'} = $exons->{$exon}{'seq_region_start'};

                                $UTR3_hash->{$gene_id_transcript_id}->{$exon}{'end'} = $exons->{$exon}{'seq_region_end'};

                                $UTR3_hash->{$gene_id_transcript_id}->{$exon}{'transcript_id'} = $gene_id_transcript_id;

                                next;

                            }

                            

                            #Check if exon already in exon_hash

                            if (exists $exon_hash->{$exon}){

                                next;

                            }else{

                                $exon_hash->{$exon} =   $exons->{$exon};

                            }

                        }

                    }else{

                        next;

                    }

                }

                

                #Run over TISses

                foreach my $g_sORF_start ( keys %{$g_sORF_starts}){

                    

                    ## Create sORF_seq

                    $sORF_seq = '';

                    $tmp_sORF_seq = '';

                    $g_sORF_id                      =   $g_sORF_starts->{$g_sORF_start}{'gene_id'};

                    $biotype                        =   $g_sORF_starts->{$g_sORF_start}{'biotype'};

                    $TIS                            =   $g_sORF_starts->{$g_sORF_start}{'start'};

                    $strand                         =   $g_sORF_starts->{$g_sORF_start}{'strand'};

                    $start_codon                    =   $g_sORF_starts->{$g_sORF_start}{'start_codon'};

                    $peak_shift                     =   $g_sORF_starts->{$g_sORF_start}{'peak_shift'};

                    $count                          =   $g_sORF_starts->{$g_sORF_start}{'count'};

                    

                    # Create max_length sequence from TIS (strand specific)

                    $tmp_sORF_seq = ($strand eq '1') ? get_sequence($chr,$TIS,$TIS + $max_DNA_length) : revdnacomp(get_sequence($chr,$TIS - $max_DNA_length,$TIS));

                    

                    # Sequence ad the end of chromosome

                    if ((length($tmp_sORF_seq) % 3) != 0){next;}

                    

                    # Create translated sequence

                    ($tr_sORF_seq,$sORF_seq,$sORF_length,$peptide_mass) = translate($tmp_sORF_seq,$min,$max);

                    

                    

                    #next if no sorf sequence

                    if($sORF_seq eq 'X'){

                        next;

                    }

                    

                    #Save in g_tmp_sorfs

                    $sorf_end   =   ($strand eq '1') ? $TIS + (($sORF_length * 3) -1) : $TIS;

                    $sorf_begin =   ($strand eq '1') ? $TIS : $TIS - (($sORF_length * 3) -1);

                    

                    $g_sorfs_tmp->{$g_sORF_start}{'gene_id'}                  =   $g_sORF_id;

                    $g_sorfs_tmp->{$g_sORF_start}{'biotype'}                  =   $biotype;

                    $g_sorfs_tmp->{$g_sORF_start}{'sorf_chr'}                 =   $chr;

                    $g_sorfs_tmp->{$g_sORF_start}{'sorf_begin'}               =   $sorf_begin;

                    $g_sorfs_tmp->{$g_sORF_start}{'sorf_end'}                 =   $sorf_end;

                    $g_sorfs_tmp->{$g_sORF_start}{'sorf_strand'}              =   $strand;

                    $g_sorfs_tmp->{$g_sORF_start}{'start_codon'}              =   $start_codon;

                    $g_sorfs_tmp->{$g_sORF_start}{'peak_shift'}               =   $peak_shift;

                    $g_sorfs_tmp->{$g_sORF_start}{'count'}                    =   $count;

                    $g_sorfs_tmp->{$g_sORF_start}{'sorf_length'}              =   $sORF_length;

                    $g_sorfs_tmp->{$g_sORF_start}{'sorf_seq'}                 =   $sORF_seq;

                    $g_sorfs_tmp->{$g_sORF_start}{'tr_sorf_seq'}              =   $tr_sORF_seq;

                    $g_sorfs_tmp->{$g_sORF_start}{'exon_overlap'}             =   '0';

                    $g_sorfs_tmp->{$g_sORF_start}{'exon_id'}                  =   '0';

                    $g_sorfs_tmp->{$g_sORF_start}{'in_frame'}                 =   'NA';

                    $g_sorfs_tmp->{$g_sORF_start}{'mass'}                     =   $peptide_mass;

                    

                    #For LTM reads start from sorf_begin -1 to calculate total LTM

                    #Otherwize, if sORF has only LTM_reads on sorf_begin -1 position total_LTM = 0 -> illegal division

                    $g_sorfs_tmp->{$g_sORF_start}{'LTM_begin'}   =    ($strand eq '1') ?  $sorf_begin - 1                        :   $g_sorfs_tmp->{$g_sORF_start}{'sorf_begin'};

                    $g_sorfs_tmp->{$g_sORF_start}{'LTM_end'}     =    ($strand eq '1') ?  $g_sorfs_tmp->{$g_sORF_start}{'sorf_end'}  :   $sorf_end + 1;

                    

                    ## Check sORF_start_position

                    my $annotation = 'False';

                    

            #If Tis overlaps aTIS -> aTIS annotation (To check presence of tal, positive control)

            if (exists $aTIS_hash->{$g_sORF_start}){

            $g_sorfs_tmp->{$g_sORF_start}{'annotation'} =   'aTIS';

            $annotation                 =   'True';

            }   

        

                    #If TIS overlaps exon -> keep track of exon_id

                    if ($annotation eq 'False'){

             foreach my $exon_id (keys %{$exon_hash}){

                               

                             if ($exon_hash->{$exon_id}{'seq_region_start'} <= $g_sORF_start && $g_sORF_start <= $exon_hash->{$exon_id}{'seq_region_end'}) {

                                 $g_sorfs_tmp->{$g_sORF_start}{'annotation'} =   'exonic';

                                 $g_sorfs_tmp->{$g_sORF_start}{'exon_id'}    =   $exon_hash->{$exon_id}{'exon_id'};

                                 $annotation                                 =   'True';

                             }

                     }

           }

                    

                    #If TIS overlaps 5UTR

                    if ($annotation eq 'False'){

                        foreach my $tr_UTR5 (keys %{$UTR5_hash}){

                            foreach my $ex_UTR5 (keys %{$UTR5_hash->{$tr_UTR5}}){

                                if($UTR5_hash->{$tr_UTR5}->{$ex_UTR5}->{'start'} <= $g_sORF_start && $g_sORF_start <= $UTR5_hash->{$tr_UTR5}->{$ex_UTR5}->{'end'}){

                            

                                    $g_sorfs_tmp->{$g_sORF_start}{'annotation'} =   '5UTR';

                                    $g_sorfs_tmp->{$g_sORF_start}{'5UTR_tr_id'} =   $UTR5_hash->{$tr_UTR5}->{$ex_UTR5}->{'transcript_id'};

                                    $annotation                                 =   'True';

                                }

                            }

                        }

                    }

                    

                    #If TIS overlaps 3UTR

                    if ($annotation eq 'False'){

                        foreach my $tr_UTR3 (keys %{$UTR3_hash}){

                            foreach my $ex_UTR3 (keys %{$UTR3_hash->{$tr_UTR3}}){

                                if($UTR3_hash->{$tr_UTR3}->{$ex_UTR3}->{'start'} <= $g_sORF_start && $g_sORF_start <= $UTR3_hash->{$tr_UTR3}->{$ex_UTR3}->{'end'}){

                            

                                    $g_sorfs_tmp->{$g_sORF_start}{'annotation'} =   '3UTR';

                                    $g_sorfs_tmp->{$g_sORF_start}{'3UTR_tr_id'} =   $UTR3_hash->{$tr_UTR3}->{$ex_UTR3}->{'transcript_id'};

                                    $annotation                                 =   'True';

                                }

                            }

                        }

                    }

                    

                    #Else TIS is intronic

                    if ($annotation eq 'False'){

                    

                        $g_sorfs_tmp->{$g_sORF_start}{'annotation'} =   'intronic';

                    }

                    

                    #For LTM reads start from sorf_begin -1 to calculate total LTM

                    #Otherwize, if sORF has only LTM_reads on sorf_begin -1 position total_LTM = 0 -> illegal division

                    $g_sorfs_tmp->{$g_sORF_start}{'LTM_begin'}   =    ($strand eq '1') ?  $sorf_begin - 1                        :   $g_sorfs_tmp->{$g_sORF_start}{'sorf_begin'};

                    $g_sorfs_tmp->{$g_sORF_start}{'LTM_end'}     =    ($strand eq '1') ?  $g_sorfs_tmp->{$g_sORF_start}{'sorf_end'}  :   $sorf_end + 1;

                }

                

                ## Check sORF overlap

                foreach my $g_sorf_tmp (keys %{$g_sorfs_tmp}){

                    

                    #Init

                    my $overlap = 'False';

                    my $overlap_hash = {};

                    my $overlap_length = 0;

                    my $sorf_length = $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'} - $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'} +1;

                    $exon_overlap = 0;

                    

                    if ($g_sorfs_tmp->{$g_sorf_tmp}{'annotation'} eq 'aTIS') {

                        next;

                    }

                    

                    if ($g_sorfs_tmp->{$g_sorf_tmp}{'annotation'} eq 'exonic'){



                        #If exonic -> Check frame and exon overlap (frame is more important as overlap)

                        foreach my $exon (keys %{$exon_hash}){

                            

                            if ($g_sorfs_tmp->{$g_sorf_tmp}{'sorf_strand'} eq '1'){

                                

                                if ($exon_hash->{$exon}{'seq_region_start'} <= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'} && $exon_hash->{$exon}{'seq_region_end'} >= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'} ){

                                    if ($exon_hash->{$exon}{'seq_region_end'} >= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'} ){

                                        $overlap = 'True';

                                        for (my $i = $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'}; $i <= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'}; $i++ ){

                                            $overlap_hash->{$i} = '1';

                                        }

                                    }elsif ($exon_hash->{$exon}{'seq_region_end'} <= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'} ){

                                        $overlap = 'True';

                                        for (my $i = $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'}; $i <= $exon_hash->{$exon}{'seq_region_end'}; $i++ ){

                                            $overlap_hash->{$i} = '1';

                                        }

                                    }

                                }elsif ($exon_hash->{$exon}{'seq_region_start'} >= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'} && $exon_hash->{$exon}{'seq_region_end'} <= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'} ){

                                    $overlap = 'True';

                                    for (my $i = $exon_hash->{$exon}{'seq_region_start'}; $i <= $exon_hash->{$exon}{'seq_region_end'}; $i++ ){

                                        $overlap_hash->{$i} = '1';

                                    }

                                }elsif ($exon_hash->{$exon}{'seq_region_start'} <= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'} && $exon_hash->{$exon}{'seq_region_end'} >= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'} ){

                                    $overlap = 'True';

                                    for (my $i = $exon_hash->{$exon}{'seq_region_start'}; $i <= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'}; $i++ ){

                                        $overlap_hash->{$i} = '1';

                                    }

                                }

                            }elsif ($g_sorfs_tmp->{$g_sorf_tmp}{'sorf_strand'} eq '-1'){

                                if ($exon_hash->{$exon}{'seq_region_start'} <= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'} && $exon_hash->{$exon}{'seq_region_end'} >= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'} ){

                                    if ($exon_hash->{$exon}{'seq_region_start'} <= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'} ){

                                        $overlap = 'True';

                                        for (my $i = $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'}; $i <= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'}; $i++ ){

                                            $overlap_hash->{$i} = '1';

                                        }

                                    }elsif ($exon_hash->{$exon}{'seq_region_start'} >= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'} ){

                                        $overlap = 'True';

                                        for (my $i = $exon_hash->{$exon}{'seq_region_start'}; $i <= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'}; $i++ ){

                                            $overlap_hash->{$i} = '1';

                                        }

                                    }

                                }elsif ($exon_hash->{$exon}{'seq_region_start'} >= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'} && $exon_hash->{$exon}{'seq_region_end'} <= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'} ){

                                    $overlap = 'True';

                                    for (my $i = $exon_hash->{$exon}{'seq_region_start'}; $i <= $exon_hash->{$exon}{'seq_region_end'}; $i++ ){

                                        $overlap_hash->{$i} = '1';

                                    }

                                }elsif ($exon_hash->{$exon}{'seq_region_start'} <= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'} && $exon_hash->{$exon}{'seq_region_end'} >= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'} ){

                                    $overlap = 'True';

                                    for (my $i = $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'}; $i <= $exon_hash->{$exon}{'seq_region_end'}; $i++ ){

                                        $overlap_hash->{$i} = '1';

                                    }

                                }

                            }

                        }

                        # percentage overlap

                        if ($overlap eq 'True'){

                            for (my $i= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'}; $i<=$g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'}; $i++){

                                if(defined $overlap_hash->{$i}){

                                    if($overlap_hash->{$i} == 1){

                                        $overlap_length++;

                                    }

                                }

                            }

                            $g_sorfs_tmp->{$g_sorf_tmp}{'exon_overlap'} = ($overlap_length / $sorf_length);

                            $exon_overlap = ($overlap_length / $sorf_length);

                        }

                    

                        ##Check phase/frame

                        #Init

                        $sorf_begin =   $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'};

                        $sorf_end   =   $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'};

                        $exon_id    =   $g_sorfs_tmp->{$g_sorf_tmp}{'exon_id'};

                        $phase      =   $exon_hash->{$exon_id}{'phase'};

                        $end_phase  =   $exon_hash->{$exon_id}{'end_phase'};

       

                        if ($g_sorfs_tmp->{$g_sorf_tmp}{'sorf_strand'} eq '1'){

                            

                            # Phase is -1 when UTR is part of exon. We deleted UTRs from exons, so phase -1 is equal to phase 0

                            if($phase eq '0' || $phase eq '-1'){

                                

                                $phase = '0';

                            

                                ## Rest van ((sorf_begin - exon_start ) / 3 ) should be equal to the phase -> in phase.

                                if (($sorf_begin -  $exon_hash->{$exon_id}{'seq_region_start'}) % 3 == $phase){

                                    $g_sorfs_tmp->{$g_sorf_tmp}{'in_frame'} = 'Yes';

                                }else{

                                    $g_sorfs_tmp->{$g_sorf_tmp}{'in_frame'} = 'No';

                                }

                            }else{

                                

                                ## 3 - Rest van ((sorf_begin - exon_start ) / 3 ) should be equal to the phase -> in phase.

                                if(3 -(($sorf_begin -  $exon_hash->{$exon_id}{'seq_region_start'}) % 3) == $phase){

                                    $g_sorfs_tmp->{$g_sorf_tmp}{'in_frame'} = 'Yes';

                                }else{

                                    $g_sorfs_tmp->{$g_sorf_tmp}{'in_frame'} = 'No';

                                }

                            }

                        }elsif ($g_sorfs_tmp->{$g_sorf_tmp}{'sorf_strand'} eq '-1'){

                        

                            # Phase is -1 when UTR is part of exon. We deleted UTRs from exons, so phase -1 is equal to phase 0

                            if($phase eq '0' || $phase eq '-1'){

                                

                                $phase = '0';

                                ## Rest van ((exon_end - sorf_end ) / 3 ) should be equal to the phase -> in phase.

                                if (($exon_hash->{$exon_id}{'seq_region_end'} - $sorf_end) % 3 == $phase){

                                    $g_sorfs_tmp->{$g_sorf_tmp}{'in_frame'} = 'Yes';

                                }else{

                                    $g_sorfs_tmp->{$g_sorf_tmp}{'in_frame'} = 'No';

                                }

                            }else{

                                

                                ## 3 - Rest van ((exon_end - sorf_end ) / 3 ) should be equal to the phase -> in phase.

                                if(3 -(($exon_hash->{$exon_id}{'seq_region_end'} - $sorf_end) % 3) == $phase){

                                    $g_sorfs_tmp->{$g_sorf_tmp}{'in_frame'} = 'Yes';

                                }else{

                                    $g_sorfs_tmp->{$g_sorf_tmp}{'in_frame'} = 'No';

                                }

                            }

                         }

                    }elsif($g_sorfs_tmp->{$g_sorf_tmp}{'annotation'} eq 'intronic'){

                        

                        #If intronic -> Check if completely intronic or overlap with exon

                        foreach my $exon (keys %{$exon_hash}){

                            

                            if ($g_sorfs_tmp->{$g_sorf_tmp}{'sorf_strand'} eq '1'){

                                if ($exon_hash->{$exon}{'seq_region_start'} >= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'} && $exon_hash->{$exon}{'seq_region_start'} <= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'} ){

                                    if ($exon_hash->{$exon}{'seq_region_end'} <= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'} ){

                                        

                                        $overlap = 'True';

                                        #Create hash with overlapping positions

                                        for (my $i = $exon_hash->{$exon}{'seq_region_start'}; $i <= $exon_hash->{$exon}{'seq_region_end'}; $i++ ){

                                            $overlap_hash->{$i} = '1';

                                        }

                                    }elsif ($exon_hash->{$exon}{'seq_region_end'} >= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'} ){

                                        

                                        $overlap = 'True';

                                        #Create hash with overlapping positions

                                        for (my $i = $exon_hash->{$exon}{'seq_region_start'}; $i <= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'}; $i++ ){

                                            $overlap_hash->{$i} = '1';

                                        }

                                    }

                                }

                            }elsif ($g_sorfs_tmp->{$g_sorf_tmp}{'sorf_strand'} eq '-1'){

                                if ($exon_hash->{$exon}{'seq_region_end'} <= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'} && $exon_hash->{$exon}{'seq_region_end'} >= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'} ){

                                    if ($exon_hash->{$exon}{'seq_region_start'} >= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'} ){

                                        

                                        $overlap = 'True';

                                        #Create hash with overlapping positions

                                        for (my $i = $exon_hash->{$exon}{'seq_region_start'}; $i <= $exon_hash->{$exon}{'seq_region_end'}; $i++ ){

                                            $overlap_hash->{$i} = '1';

                                        }

                                    }elsif ($exon_hash->{$exon}{'seq_region_start'} <= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'} ){

                                        

                                        $overlap = 'True';

                                        #Create hash with overlapping positions

                                        for (my $i = $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'}; $i <= $exon_hash->{$exon}{'seq_region_end'}; $i++ ){

                                            $overlap_hash->{$i} = '1';

                                        }

                                    }

                                }

                            }

                        }

                        

                        #percentage overlap

                        if ($overlap eq 'True'){

                            for (my $i= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'}; $i<=$g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'}; $i++){

                                if(defined $overlap_hash->{$i}){

                                    if($overlap_hash->{$i} == 1){

                                        $overlap_length++;

                                    }

                                }

                            }

                            $g_sorfs_tmp->{$g_sorf_tmp}{'exon_overlap'} = ($overlap_length / $sorf_length);

                            $exon_overlap = ($overlap_length / $sorf_length);

                        }

                    }elsif($g_sorfs_tmp->{$g_sorf_tmp}{'annotation'} eq '5UTR'){

                        

                        #If 5UTR -> Check if completely UTR or overlap with exon

                        foreach my $exon (keys %{$exon_hash}){

                            

                            if ($g_sorfs_tmp->{$g_sorf_tmp}{'sorf_strand'} eq '1'){

                                if ($exon_hash->{$exon}{'seq_region_start'} >= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'} && $exon_hash->{$exon}{'seq_region_start'} <= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'} ){

                                    if ($exon_hash->{$exon}{'seq_region_end'} <= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'} ){

                                        

                                        $overlap = 'True';

                                        #Create hash with overlapping positions

                                        for (my $i = $exon_hash->{$exon}{'seq_region_start'}; $i <= $exon_hash->{$exon}{'seq_region_end'}; $i++ ){

                                            $overlap_hash->{$i} = '1';

                                        }

                                    }elsif ($exon_hash->{$exon}{'seq_region_end'} >= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'} ){

                                        

                                        $overlap = 'True';

                                        #Create hash with overlapping positions

                                        for (my $i = $exon_hash->{$exon}{'seq_region_start'}; $i <= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'}; $i++ ){

                                            $overlap_hash->{$i} = '1';

                                        }

                                    }

                                }

                            }elsif ($g_sorfs_tmp->{$g_sorf_tmp}{'sorf_strand'} eq '-1'){

                                

                                if ($exon_hash->{$exon}{'seq_region_end'} <= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'} && $exon_hash->{$exon}{'seq_region_end'} >= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'} ){

                                    if ($exon_hash->{$exon}{'seq_region_start'} >= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'} ){

                                        

                                        $overlap = 'True';

                                        #Create hash with overlapping positions

                                        for (my $i = $exon_hash->{$exon}{'seq_region_start'}; $i <= $exon_hash->{$exon}{'seq_region_end'}; $i++ ){

                                            $overlap_hash->{$i} = '1';

                                        }

                                    }elsif ($exon_hash->{$exon}{'seq_region_start'} <= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'} ){

                                        

                                        $overlap = 'True';

                                        #Create hash with overlapping positions

                                        for (my $i = $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'}; $i <= $exon_hash->{$exon}{'seq_region_end'}; $i++ ){

                                            $overlap_hash->{$i} = '1';

                                        }

                                    }

                                }

                            }

                        }

                        

                        #percentage overlap

                        if ($overlap eq 'True'){

                            for (my $i= $g_sorfs_tmp->{$g_sorf_tmp}{'sorf_begin'}; $i<=$g_sorfs_tmp->{$g_sorf_tmp}{'sorf_end'}; $i++){

                                if(defined $overlap_hash->{$i}){

                                    if($overlap_hash->{$i} == 1){

                                        $overlap_length++;

                                    }

                                }

                            }

                            $g_sorfs_tmp->{$g_sorf_tmp}{'exon_overlap'} = ($overlap_length / $sorf_length);

                            $exon_overlap = ($overlap_length / $sorf_length);

                        }

                    }else{

                        

                        #3UTR and no other downstream overlap possible

                        $g_sorfs_tmp->{$g_sorf_tmp}{'exon_overlap'} = 0;

                        $exon_overlap = '0';

                    }

                    

                    #Save in g_sorfs

                    $g_sorfs->{$g_sorf_tmp} = $g_sorfs_tmp->{$g_sorf_tmp};

                }

            #Else gene is non-coding

            }else{

                

                #Run over TISses

                foreach my $g_sORF_start ( keys %{$g_sORF_starts}){

                  

                    $sORF_seq = '';

                    $tmp_sORF_seq = '';

                    $g_sORF_id                      =   $g_sORF_starts->{$g_sORF_start}{'gene_id'};

                    $biotype                        =   $g_sORF_starts->{$g_sORF_start}{'biotype'};

                    $TIS                            =   $g_sORF_starts->{$g_sORF_start}{'start'};

                    $strand                         =   $g_sORF_starts->{$g_sORF_start}{'strand'};

                    $start_codon                    =   $g_sORF_starts->{$g_sORF_start}{'start_codon'};

                    $peak_shift                     =   $g_sORF_starts->{$g_sORF_start}{'peak_shift'};;

                    $count                          =   $g_sORF_starts->{$g_sORF_start}{'count'};

                    

                    # Create max_length sequence from TIS (strand specific)

                    $tmp_sORF_seq = ($strand eq '1') ? get_sequence($chr,$TIS,$TIS + $max_DNA_length) : revdnacomp(get_sequence($chr,$TIS - $max_DNA_length,$TIS));

                    

                    # Sequence ad the end of chromosome

                    if ((length($tmp_sORF_seq) % 3) != 0){next;}

                    

                    # Create translated sequence

                    ($tr_sORF_seq,$sORF_seq,$sORF_length,$peptide_mass) = translate($tmp_sORF_seq,$min,$max);

                    

                    #next if no sorf sequence

                    if($sORF_seq eq 'X'){next;}

                    

                    #Save in g_sorfs

                    $sorf_end   =   ($strand eq '1') ? $TIS + (($sORF_length * 3) -1) : $TIS;

                    $sorf_begin =   ($strand eq '1') ? $TIS : $TIS - (($sORF_length * 3) -1);

                    

                    $g_sorfs->{$g_sORF_start}{'gene_id'}                  =   $g_sORF_id;

                    $g_sorfs->{$g_sORF_start}{'biotype'}                  =   $biotype;

                    $g_sorfs->{$g_sORF_start}{'sorf_chr'}                 =   $chr;

                    $g_sorfs->{$g_sORF_start}{'sorf_begin'}               =   $sorf_begin;

                    $g_sorfs->{$g_sORF_start}{'sorf_end'}                 =   $sorf_end;

                    $g_sorfs->{$g_sORF_start}{'sorf_strand'}              =   $strand;

                    $g_sorfs->{$g_sORF_start}{'start_codon'}              =   $start_codon;

                    $g_sorfs->{$g_sORF_start}{'peak_shift'}               =   $peak_shift;

                    $g_sorfs->{$g_sORF_start}{'count'}                    =   $count;

                    $g_sorfs->{$g_sORF_start}{'sorf_length'}              =   $sORF_length;

                    $g_sorfs->{$g_sORF_start}{'sorf_seq'}                 =   $sORF_seq;

                    $g_sorfs->{$g_sORF_start}{'tr_sorf_seq'}              =   $tr_sORF_seq;

                    $g_sorfs->{$g_sORF_start}{'annotation'}               =   'ncRNA';

                    $g_sorfs->{$g_sORF_start}{'exon_overlap'}             =   0;

                    $g_sorfs->{$g_sORF_start}{'in_frame'}                 =   'NA';

                    $g_sorfs->{$g_sORF_start}{'mass'}                     =   $peptide_mass;

                    

                    #For LTM reads start from sorf_begin -1 to calculate total LTM

                    #Otherwize, if sORF has only LTM_reads on sorf_begin -1 position total_LTM = 0 -> illegal division

                    $g_sorfs->{$g_sORF_start}{'LTM_begin'}   =    ($strand eq '1') ?  $sorf_begin - 1                        :   $g_sorfs->{$g_sORF_start}{'sorf_begin'};

                    $g_sorfs->{$g_sORF_start}{'LTM_end'}     =    ($strand eq '1') ?  $g_sorfs->{$g_sORF_start}{'sorf_end'}  :   $sorf_end + 1;

                }

            }

        }

        

        # Match reads to genic sORFs

        $g_sorfs = match_reads_to_genic_sORFs($g_sorfs,$CHX_for,$CHX_rev,$LTM_for,$LTM_rev);

        

        # Loop over sORFs in $g_sORFs

        foreach my $g_sorf_id (keys %{$g_sorfs}){

        

            # Get sORF specific reads

            my ($LTM_reads,$CHX_reads) = get_overlapping_reads($g_sorf_id,$g_sorfs,$chr,$CHX_for,$CHX_rev,$LTM_for,$LTM_rev);

            

            # Calculate FPKM, coverage and R

            ($FPKM,$coverage_sORF,$R_sORF,$coverage_uniformity) = calculate_FPKM_and_coverage($g_sorfs,$g_sorf_id,$Mreads,$LTM_reads,$CHX_reads,$mean_length_fastq1,$mean_length_fastq2);

        

            #Check if R value and coverage are higher or equal as $R and $coverage

            if($R_sORF < $R){next;}

            if($coverage_sORF <= 0){next;}

            

            # Save sORF in csv

            $g_sORF_id                  =   $g_sorfs->{$g_sorf_id}{'gene_id'};

            $sorf_chr                   =   $g_sorfs->{$g_sorf_id}{'sorf_chr'};

            $sorf_begin                 =   $g_sorfs->{$g_sorf_id}{'sorf_begin'};

            $sorf_end                   =   $g_sorfs->{$g_sorf_id}{'sorf_end'};

            $sorf_strand                =   $g_sorfs->{$g_sorf_id}{'sorf_strand'};

            $start_codon                =   $g_sorfs->{$g_sorf_id}{'start_codon'};

            $peak_shift                 =   $g_sorfs->{$g_sorf_id}{'peak_shift'};

            $count                      =   $g_sorfs->{$g_sorf_id}{'count'};

            $sORF_length                =   $g_sorfs->{$g_sorf_id}{'sorf_length'};

            $sORF_seq                   =   $g_sorfs->{$g_sorf_id}{'sorf_seq'};

            $tr_sORF_seq                =   $g_sorfs->{$g_sorf_id}{'tr_sorf_seq'};

            $annotation                 =   $g_sorfs->{$g_sorf_id}{'annotation'};

            $biotype                    =   $g_sorfs->{$g_sorf_id}{'biotype'};

            $exon_overlap               =   $g_sorfs->{$g_sorf_id}{'exon_overlap'};

            $in_frame                   =   $g_sorfs->{$g_sorf_id}{'in_frame'};

            $peptide_mass               =   $g_sorfs->{$g_sorf_id}{'mass'};

            

            if($exon_overlap >= $ex_overlap){next;}

            if($annotation eq 'aTIS'){next;}

            

            $exon_overlap = sprintf("%.2f",$exon_overlap);

            

            # Check how to implement

            my $out_cnt        =   '1';

            # Should becom $_->{'SNP'}

            my $snp_tmp = '0';

            

            #Add extra columns to be compatible with intergenic_sORF_calling

            my $downstream_gene_distance    =   'NA';

            my $upstream_gene_distance      =   'NA';

            

            print TMP_db $g_sORF_id.",".$biotype.",".$sorf_chr.",".$sorf_strand.",".$sorf_begin.",".$sorf_end.",".$sORF_length.",".$peptide_mass.",".$start_codon.",".$downstream_gene_distance.",".$upstream_gene_distance.",".$annotation.",".$exon_overlap.",".$in_frame.",".$peak_shift.",".$count.",".$R_sORF.",".$coverage_sORF.",".$coverage_uniformity.",".$FPKM.",".$snp_tmp.",".$sORF_seq.",".$tr_sORF_seq."\n";

        }

        ### Finish childs

        print "     * Finished translating chromosome ".$chr."\n";

        $dbh->disconnect();

        $pm->finish;

    }

    #Waiting for all childs to finish

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

    

    # Init

    my $dbh     =   dbh($dsn,$us,$pw);

    my $table   =   "TIS_sORFs_".$id."_transcripts";

    

    # Drop table if exists transcript table

    my $query = "DROP TABLE IF EXISTS `".$table."`";

    $dbh->do($query);

    

    # Create table tmp

    $query = "CREATE TABLE IF NOT EXISTS `".$table."_tmp` (

    `id` varchar(128) NOT NULL default '',

    `biotype` varchar(128) NOT NULL default '',

    `chr` char(50) NOT NULL default '',

    `strand` int(2) NOT NULL default '',

    `sorf_begin` int(10) NOT NULL default '',

    `sorf_end` int(10) NOT NULL default '',

    `sorf_length` int(10) NOT NULL default '',

    `mass` decimal(10,5) NOT NULL default '0',

    `start_codon` varchar(128) NOT NULL default '',

    `downstream_gene_distance` int(10) NOT NULL default '',

    `upstream_gene_distance` int(10) NOT NULL default 'NA',

    `annotation` varchar(128) NOT NULL default 'NA',

    `exon_overlap` varchar(128) NOT NULL default 'NA',

    `in_frame` varchar(128) NOT NULL default '',

    `peak_shift` int(2) NOT NULL default '',

    `count` float default NULL,

    `Rltm_min_Rchx` decimal(11,8) NOT NULL default '0',

    `coverage` decimal(11,8) NOT NULL default '0',

    `coverage_uniformity` decimal(11,8) NOT NULL default '0',

    `RPKM` decimal(11,8) NOT NULL default '0',

    `SNP_` decimal(11,8) NOT NULL default '0',

    `tr_seq` TEXT NOT NULL default '',

    `aa_seq` TEXT NOT NULL default '' )"  ;

    $dbh->do($query);

    

    # Create table

    $query = "CREATE TABLE IF NOT EXISTS `".$table."` (

    `sorf_id` INTEGER PRIMARY KEY AUTOINCREMENT,

    `id` varchar(128) NOT NULL default '',

    `biotype` varchar(128) NOT NULL default '',

    `chr` char(50) NOT NULL default '',

    `strand` int(2) NOT NULL default '',

    `sorf_begin` int(10) NOT NULL default '',

    `sorf_end` int(10) NOT NULL default '',

    `sorf_length` int(10) NOT NULL default '',

    `mass` decimal(10,5) NOT NULL default '0',

    `start_codon` varchar(128) NOT NULL default '',

    `downstream_gene_distance` int(10) NOT NULL default '',

    `upstream_gene_distance` int(10) NOT NULL default 'NA',

    `annotation` varchar(128) NOT NULL default 'NA',

    `exon_overlap` varchar(128) NOT NULL default 'NA',

    `in_frame` varchar(128) NOT NULL default '',

    `peak_shift` int(2) NOT NULL default '',

    `count` float default NULL,

    `Rltm_min_Rchx` decimal(11,8) NOT NULL default '0',

    `coverage` decimal(11,8) NOT NULL default '0',

    `coverage_uniformity` decimal(11,8) NOT NULL default '0',

    `RPKM` decimal(11,8) NOT NULL default '0',

    `SNP_` decimal(11,8) NOT NULL default '0',

    `tr_seq` TEXT NOT NULL default '',

    `aa_seq` TEXT NOT NULL default '' )"  ;

    $dbh->do($query);  

    

    # Store

    foreach my $chr (sort keys %{$chrs}){

        system("sqlite3 -separator , ".$sqlite_db." \".import ".$TMP."/".$id."_".$chr."_sORF_tmp.csv ".$table."_tmp\"")== 0 or die "system failed: $?";

    }

    

    #Save in transcripts with auto_increment and drop tmp

    $query = "INSERT INTO `".$table."` (id,biotype,chr,strand,sorf_begin,sorf_end,sorf_length,mass,start_codon,downstream_gene_distance,upstream_gene_distance,annotation,exon_overlap,in_frame,peak_shift,count,Rltm_min_Rchx,coverage,coverage_uniformity,RPKM,SNP_,tr_seq,aa_seq) SELECT * FROM `".$table."_tmp`";

    $dbh->do($query);

    $query = "drop table `".$table."_tmp` ";

    $dbh->do($query);

    

    #Disconnect dbh

    $dbh->disconnect();

    

    # Unlink tmp csv files

    foreach my $chr (sort keys %{$chrs}){

        unlink $TMP."/".$id."_".$chr."_sORF_tmp.csv";

    }



}



### Calculate FPKM and coverage ###



sub calculate_FPKM_and_coverage{

    

    # Catch

    my $sorfs               =   $_[0];

    my $sorf_id             =   $_[1];

    my $Mreads              =   $_[2];

    my $LTM_reads           =   $_[3];

    my $CHX_reads           =   $_[4];

    my $mean_length_fastq1  =   $_[5];

    my $mean_length_fastq2  =   $_[6];

    

    # Init

    my $sorf_end    =   $sorfs->{$sorf_id}{'sorf_end'};

    my $sorf_begin  =   $sorfs->{$sorf_id}{'sorf_begin'};

    my $count       =   $sorfs->{$sorf_id}{'count'};

    my ($FPKM,$coverage_sORF,$read_begin,$read_end,$Rltm,$Rchx,$R_sORF,$coverage_uniformity);

    my $sorf_length =   $sorf_end-$sorf_begin +1;

    my $dbh = dbh($dsn_results,$us_results,$pw_results);

    my $hits = 0;

    my $hits_total = 0;

    my $covered_positions = 0;

    my $position = 0;

    my $coverage = {};

    

    #Start on 1 to avoid illegal division by 0.

    my $counts_1 = 1;

    my $counts_2 = 1;

    

    # Total number of LTM and CHX reads in gene

    my $total_LTM = 0;

    my $total_CHX = 0;



    foreach my $key (keys %{$CHX_reads}){

        $total_CHX += $CHX_reads->{$key}{'count'};

    }

    

    foreach my $key (keys %{$LTM_reads}){

        $total_LTM += $LTM_reads->{$key}{'count'};

    }

    

    # Run over CHX ribo_profiles

    foreach my $key (keys %{$CHX_reads}) {

        $coverage->{$key} = '1';

    }

    

    #Calculate percent coverage and FPKM for sORF

    for (my $i= $sorf_begin; $i<$sorf_end+1; $i++){

       if(defined $coverage->{$i}){

           if($coverage->{$i} == 1){

               $covered_positions++;

           }

       }

    }

    $coverage_sORF = $covered_positions/$sorf_length;

    $coverage_sORF = sprintf("%.2f",$coverage_sORF);

    

    # Calculate coverage uniformity

    my $half_sorf_length = int($sorf_length/2);



    for (my $i= $sorf_begin; $i<$sorf_begin + $half_sorf_length +1; $i++){

        if(defined $coverage->{$i}){

            if($coverage->{$i} == 1){

                $counts_1 += $CHX_reads->{$i}{'count'};

            }

        }

    }

    for (my $i= $sorf_begin + $half_sorf_length; $i<$sorf_end +1; $i++){

        if(defined $coverage->{$i}){

            if($coverage->{$i} == 1){

                $counts_2 += $CHX_reads->{$i}{'count'};

            }

        }

    }

    $coverage_uniformity = $counts_2/$counts_1;

    $coverage_uniformity = sprintf("%.2f",$coverage_uniformity);

    $FPKM = $total_CHX / (($sorf_length / 1000) * $Mreads);

    $FPKM = sprintf("%.2f",$FPKM);

    

    #Calculate RLTM - RCHX

    #Rk = (Xk/Nk) x 10 (k = LTM, CHX), Xk number of reads on that position in data k, Nk total number of reads for transcript.

    

    # Calculate Rltm

    $Rltm = ($count / ($total_LTM/$mean_length_fastq2)) * 10;

    if(exists $CHX_reads->{$sorf_begin}{'count'}){

        $Rchx = ($CHX_reads->{$sorf_begin}{'count'}/ ($total_CHX/$mean_length_fastq1));

        $Rchx = $Rchx*10;

    }else{

        $Rchx=0;

    }

    $R_sORF = $Rltm - $Rchx;

    $R_sORF = sprintf("%.2f",$R_sORF);

    

    #Return

    return($FPKM,$coverage_sORF,$R_sORF,$coverage_uniformity);

}



### translation subroutine ###

    

sub translate {

    

    # Add SNPs if available 

    #       -> Watch sequence length (if indel, add sequence)

    #       -> Watch positions, sort SNP Hash, if dels, downstream position changes

    

    #Catch

    my $tmp_sORF_seq    =   $_[0];

    my $min             =   $_[1];

    my $max             =   $_[2];

    

    #Init

    my ($i,$triplet,$tr_sORF_seq,$sORF_seq);

    my $AA_seq      =   '';

    my $DNA_seq     =   '';

    my $AA          =   '';

    my $sORF_length =   0;

       



    for($i=0;$i<=length($tmp_sORF_seq)-2;$i+=3){

        $sORF_length++;

        $triplet = substr($tmp_sORF_seq,$i,3);

        

        #Stop translation if DNAseq contains N

        if ($triplet =~ /N/){last;}

        

        #Translate DNA to AA

        $AA = DNAtoAA($triplet);

        $AA_seq     =   $AA_seq . $AA;

        $DNA_seq    =   $DNA_seq . $triplet;

        

        #Stop translation if STOP-codon is reached

        if ($AA eq '*'){last;}

    }

    

    # Check length of sORF_AA_seq

    if($sORF_length >= $min){

        if ($AA_seq =~ /\*$/ && $AA_seq ne '') {

           

            #Replace near-cognate start to cognate methionine...

            $AA_seq = (substr($AA_seq,0,1) ne 'M') ? 'M'.substr($AA_seq,1) : $AA_seq;

            

            $tr_sORF_seq    =   $AA_seq;

            $sORF_seq       =   $DNA_seq;

            

            

            #Calculate monoIsotopic mass of peptides

            my $peptide_mass = 0;

            

            ## @peptide_AAsplit peptide on AA

            my @peptide_AA = split('',$AA_seq);

            

            my %AA_mass = (

            'A'=>71.03711,  'R'=>156.10111,  'D'=>115.02694,  'N'=>114.04293,

            'C'=>103.00919,  'E'=>129.04259,  'Q'=>128.05858,  'G'=>57.02146,

            'H'=>137.05891,  'I'=>113.08406,  'L'=>113.08406,  'K'=>128.09496,

            'M'=>131.04049,  'F'=>147.06841,  'P'=>97.05276,  'S'=>87.03203,

            'T'=>101.04768,  'W'=>186.07931,  'Y'=>163.06333,  'V'=>99.06841, '*'=>0,

            );

            

            for my $AA( @peptide_AA ) {

                $peptide_mass += $AA_mass{$AA};

            }

            

            ## add mass of N-terminus H and C-terminus OH (1.00783 & 17.00274)

            $peptide_mass += 18.01057;

            

            #Return

            return($tr_sORF_seq,$sORF_seq,$sORF_length,$peptide_mass);

        }

    }

    

    #Bad sORF sequence (too short or no stop codon)

    $tr_sORF_seq    =   'X';

    $sORF_seq       =   'X';

    

    

    #Return

    return($tr_sORF_seq,$sORF_seq,$sORF_length,0);

}

    

### Get overlapping reads ###



sub get_overlapping_reads{

    

    # Catch

    my $sorf_id         =   $_[0];

    my $sorfs           =   $_[1];

    my $chr             =   $_[2];

    my $CHX_for         =   $_[3];

    my $CHX_rev         =   $_[4];

    my $LTM_for         =   $_[5];

    my $LTM_rev         =   $_[6];

    

    #Init

    my %LTM_reads = ();

    my %CHX_reads = ();

    

    my $sorf = $sorfs->{$sorf_id};

    my $strand = $sorf->{'sorf_strand'};

    my ($CHX_all,$LTM_all);

    

    # Sense specific

    if ($strand eq '1') { 

        $CHX_all = $CHX_for;

        $LTM_all = $LTM_for;

    }

    if ($strand eq '-1') {

        $CHX_all = $CHX_rev;

        $LTM_all = $LTM_rev;

    }

    

    #Get gene specific reads if exist

    if($sorfs->{$sorf_id}{'CHX'}){

        my @keys = @{$sorfs->{$sorf_id}{'CHX'}};

        @CHX_reads{@keys} = @{$CHX_all}{@keys};

    }

    

    if($sorfs->{$sorf_id}{'LTM'}){

        my @keys = @{$sorfs->{$sorf_id}{'LTM'}};

        @LTM_reads{@keys} = @{$LTM_all}{@keys};

    }

    

    #Return

    return(dclone \%LTM_reads,\%CHX_reads);

}



### Match reads to genes ##



sub match_reads_to_genic_sORFs{

    

    #Catch

    my $g_sorfs     =   $_[0];

    my $CHX_for     =   $_[1];

    my $CHX_rev     =   $_[2];

    my $LTM_for     =   $_[3];

    my $LTM_rev     =   $_[4];

    

    #Init

    my @window = ();

    my ($g_sorf_for,$g_sorf_rev,$g_sorf_for_LTM,$g_sorf_rev_LTM);

    

    #Split g_sorfs in forward and reverse arrays

    foreach my $g_sorf_id (sort { $g_sorfs->{$a}{'sorf_begin'} <=> $g_sorfs->{$b}{'sorf_begin'} } keys %{$g_sorfs}){

        if ($g_sorfs->{$g_sorf_id}{'sorf_strand'} eq '1'){

            push (@$g_sorf_for,$g_sorf_id);

            push (@$g_sorf_for_LTM,$g_sorf_id);

        }else{

            push (@$g_sorf_rev,$g_sorf_id);

            push (@$g_sorf_rev_LTM,$g_sorf_id);

        }

    }

    

    # Loop over CHX_forward

    foreach my $key (sort {$a <=> $b} keys %{$CHX_for}){

        

        # Push all g_sorf_ids to @window where g_sorf_start < window_pos

        foreach my $g_sorf_for_id (@$g_sorf_for){

            if($g_sorfs->{$g_sorf_for_id}{'sorf_begin'} <= $key){

                push(@window,$g_sorf_for_id);

            }else{last;}

        }

        # Get rid of g_sorf_for elements already in @window

        @$g_sorf_for = grep { $g_sorfs->{$_}{'sorf_begin'} > $key} @$g_sorf_for;

        

        # Get rid of g_sorf_ids in @$window where g_end < window_pos

        @window = grep { $g_sorfs->{$_}{'sorf_end'} >= $key} @window;

        

        # Loop over window and add read position to window_g_sorfs

        foreach my $window_id (@window){

            push(@{$g_sorfs->{$window_id}{'CHX'}},$key);

        }

    }

    

    #Empty @window

    @window = ();

    

    # Loop over CHX_reverse

    foreach my $key (sort {$a <=> $b} keys %{$CHX_rev}){

        

        # Push all g_sorf_ids to @window where g_sorf_start < window_pos

        foreach my $g_sorf_rev_id (@$g_sorf_rev){

            if($g_sorfs->{$g_sorf_rev_id}{'sorf_begin'} <= $key){

                push(@window,$g_sorf_rev_id);

            }else{last;}

        }

        # Get rid of g_sorf_for elements already in @window

        @$g_sorf_rev = grep { $g_sorfs->{$_}{'sorf_begin'} > $key} @$g_sorf_rev;

        

        # Get rid of g_sorf_ids in @$window where g_sorf_end < window_pos

        @window = grep { $g_sorfs->{$_}{'sorf_end'} >= $key} @window;

        

        # Loop over window and add read position to window_g_sorfs

        foreach my $window_id (@window){

            push(@{$g_sorfs->{$window_id}{'CHX'}},$key);

        }

    }

    

    #Empty @window

    @window = ();

    

    # Loop over LTM_forward

    foreach my $key (sort {$a <=> $b} keys %{$LTM_for}){

        

        # Push all g_sorf_ids to @window where g_sorf_start < window_pos

        foreach my $g_sorf_for_LTM_id (@$g_sorf_for_LTM){

            if($g_sorfs->{$g_sorf_for_LTM_id}{'LTM_begin'} <= $key){

                push(@window,$g_sorf_for_LTM_id);

            }else{last;}

        }

        # Get rid of g_sorf_for elements already in @window

        @$g_sorf_for_LTM = grep { $g_sorfs->{$_}{'LTM_begin'} > $key} @$g_sorf_for_LTM;

        

        # Get rid of g_sorf_ids in @$window where g_sorf_end < window_pos

        @window = grep { $g_sorfs->{$_}{'LTM_end'} >= $key} @window;

        

        # Loop over window and add read position to window_g_sorfs

        foreach my $window_id (@window){

            push(@{$g_sorfs->{$window_id}{'LTM'}},$key);

        }

    }

    

    #Empty @window

    @window = ();

    

    # Loop over LTM_reverse

    foreach my $key (sort {$a <=> $b} keys %{$LTM_rev}){

        

        # Push all g_sorf_ids to @window where g_sorf_start < window_pos

        foreach my $g_sorf_rev_LTM_id (@$g_sorf_rev_LTM){

            if($g_sorfs->{$g_sorf_rev_LTM_id}{'LTM_begin'} <= $key){

                push(@window,$g_sorf_rev_LTM_id);

            }else{last;}

        }

        # Get rid of g_sorf_for elements already in @window

        @$g_sorf_rev_LTM = grep { $g_sorfs->{$_}{'LTM_begin'} > $key} @$g_sorf_rev_LTM;

        

        # Get rid of g_sorf_ids in @$window where g_sorf_end < window_pos

        @window = grep { $g_sorfs->{$_}{'LTM_end'} >= $key} @window;

        

        # Loop over window and add read position to window_g_sorfs

        foreach my $window_id (@window){

            push(@{$g_sorfs->{$window_id}{'LTM'}},$key);

        }

    }

    

    #Return

    return($g_sorfs);

}



### Match reads to intergenes ##



sub match_reads_to_intergenic_sORFs{

    

    #Catch

    my $ig_sorfs    =   $_[0];

    my $CHX_for     =   $_[1];

    my $CHX_rev     =   $_[2];

    my $LTM_for     =   $_[3];

    my $LTM_rev     =   $_[4];

    

    #Init

    my @window = ();

    my ($ig_sorf_for,$ig_sorf_rev,$ig_sorf_for_LTM,$ig_sorf_rev_LTM);

    

    #Split ig_sorfs in forward and reverse arrays

    foreach my $ig_sorf_id (sort { $ig_sorfs->{$a}{'sorf_begin'} <=> $ig_sorfs->{$b}{'sorf_begin'} } keys %{$ig_sorfs}){

        if ($ig_sorfs->{$ig_sorf_id}{'sorf_strand'} eq '1'){

            push (@$ig_sorf_for,$ig_sorf_id);

            push (@$ig_sorf_for_LTM,$ig_sorf_id);

        }else{

            push (@$ig_sorf_rev,$ig_sorf_id);

            push (@$ig_sorf_rev_LTM,$ig_sorf_id);

        }

    }

    

    # Loop over CHX_forward

    foreach my $key (sort {$a <=> $b} keys %{$CHX_for}){

        

        # Push all ig_sorf_ids to @window where ig_sorf_start < window_pos

        foreach my $ig_sorf_for_id (@$ig_sorf_for){

            if($ig_sorfs->{$ig_sorf_for_id}{'sorf_begin'} <= $key){

                push(@window,$ig_sorf_for_id);

            }else{last;}

        }

        # Get rid of ig_sorf_for elements already in @window

        @$ig_sorf_for = grep { $ig_sorfs->{$_}{'sorf_begin'} > $key} @$ig_sorf_for;

        

        # Get rid of ig_sorf_ids in @$window where ig_end < window_pos

        @window = grep { $ig_sorfs->{$_}{'sorf_end'} >= $key} @window;

        

        # Loop over window and add read position to window_ig_sorfs

        foreach my $window_id (@window){

            push(@{$ig_sorfs->{$window_id}{'CHX'}},$key);

        }

    }

    

    #Empty @window

    @window = ();

    

    # Loop over CHX_reverse

    foreach my $key (sort {$a <=> $b} keys %{$CHX_rev}){

        

        # Push all ig_sorf_ids to @window where ig_sorf_start < window_pos

        foreach my $ig_sorf_rev_id (@$ig_sorf_rev){

            if($ig_sorfs->{$ig_sorf_rev_id}{'sorf_begin'} <= $key){

                push(@window,$ig_sorf_rev_id);

            }else{last;}

        }

        # Get rid of ig_sorf_for elements already in @window

        @$ig_sorf_rev = grep { $ig_sorfs->{$_}{'sorf_begin'} > $key} @$ig_sorf_rev;

        

        # Get rid of ig_sorf_ids in @$window where ig_sorf_end < window_pos

        @window = grep { $ig_sorfs->{$_}{'sorf_end'} >= $key} @window;

        

        # Loop over window and add read position to window_ig_sorfs

        foreach my $window_id (@window){

            push(@{$ig_sorfs->{$window_id}{'CHX'}},$key);

        }

    }

    

    #Empty @window

    @window = ();

    

    # Loop over LTM_forward

    foreach my $key (sort {$a <=> $b} keys %{$LTM_for}){

        

        # Push all ig_sorf_ids to @window where ig_sorf_start < window_pos

        foreach my $ig_sorf_for_LTM_id (@$ig_sorf_for_LTM){

            if($ig_sorfs->{$ig_sorf_for_LTM_id}{'LTM_begin'} <= $key){

                push(@window,$ig_sorf_for_LTM_id);

            }else{last;}

        }

        # Get rid of ig_sorf_for elements already in @window

        @$ig_sorf_for_LTM = grep { $ig_sorfs->{$_}{'LTM_begin'} > $key} @$ig_sorf_for_LTM;

        

        # Get rid of ig_sorf_ids in @$window where ig_sorf_end < window_pos

        @window = grep { $ig_sorfs->{$_}{'LTM_end'} >= $key} @window;

        

        # Loop over window and add read position to window_ig_sorfs

        foreach my $window_id (@window){

            push(@{$ig_sorfs->{$window_id}{'LTM'}},$key);

        }

    }

    

    #Empty @window

    @window = ();

    

    # Loop over LTM_reverse

    foreach my $key (sort {$a <=> $b} keys %{$LTM_rev}){

        

        # Push all ig_sorf_ids to @window where ig_sorf_start < window_pos

        foreach my $ig_sorf_rev_LTM_id (@$ig_sorf_rev_LTM){

            if($ig_sorfs->{$ig_sorf_rev_LTM_id}{'LTM_begin'} <= $key){

                push(@window,$ig_sorf_rev_LTM_id);

            }else{last;}

        }

        # Get rid of ig_sorf_for elements already in @window

        @$ig_sorf_rev_LTM = grep { $ig_sorfs->{$_}{'LTM_begin'} > $key} @$ig_sorf_rev_LTM;

        

        # Get rid of ig_sorf_ids in @$window where ig_sorf_end < window_pos

        @window = grep { $ig_sorfs->{$_}{'LTM_end'} >= $key} @window;

        

        # Loop over window and add read position to window_ig_sorfs

        foreach my $window_id (@window){

            push(@{$ig_sorfs->{$window_id}{'LTM'}},$key);

        }

    }

    

    #Return

    return($ig_sorfs);

}



### Get reads ###



sub get_reads {

    

    # Catch

    my $chr             =   $_[0];

    

    # Init

    my $dbh         = dbh($dsn_results,$us_results,$pw_results);

    my $CHX_for     = {};

    my $LTM_for     = {};

    my $CHX_rev     = {};

    my $LTM_rev     = {};

    

    # Get Reads

    my $query = "SELECT * FROM count_fastq1 WHERE chr = '$chr' and strand = '1'";

    my $sth = $dbh->prepare($query);

    $sth->execute();

    $CHX_for = $sth->fetchall_hashref('start');

    

    $query = "SELECT * FROM count_fastq2 WHERE chr = '$chr' and strand = '1'";

    $sth = $dbh->prepare($query);

    $sth->execute();

    $LTM_for = $sth->fetchall_hashref('start');

    

    $query = "SELECT * FROM count_fastq1 WHERE chr = '$chr' and strand = '-1'";

    $sth = $dbh->prepare($query);

    $sth->execute();

    $CHX_rev = $sth->fetchall_hashref('start');

    

    $query = "SELECT * FROM count_fastq2 WHERE chr = '$chr' and strand = '-1'";

    $sth = $dbh->prepare($query);

    $sth->execute();

    $LTM_rev = $sth->fetchall_hashref('start');

    

    #Disconnect DBH

    $dbh->disconnect();

    

    # Return

    return($CHX_for,$LTM_for,$CHX_rev,$LTM_rev);    

}



### Get transcript translation data ###



sub get_translation_data {



    # Catch

    my $dsn_ENS         =   $_[0];

    my $us_ENS          =   $_[1];

    my $pw_ENS          =   $_[2];

    my $tr_id           =   $_[3];

    my $trs             =   $_[4];

    my $seq_region_id   =   $_[5];

    my $chr             =   $_[6];

    

    #Init

    my $dbh_ENS = dbh($dsn_ENS,$us_ENS,$pw_ENS);

    my $exons = {};

    my $trans = [];

    

    # Get translation data for transcrip_id if available

    my $query = "SELECT start_exon_id,end_exon_id,seq_start,seq_end FROM translation where transcript_id = '".$tr_id."'";

    my $sth = $dbh_ENS->prepare($query);

    $sth->execute();

    $trans = $sth->fetchrow_arrayref();

    

    # Get exons from ENS DB

    $query = "SELECT a.exon_id,a.rank,b.seq_region_start,b.seq_region_end,b.phase,b.end_phase FROM exon_transcript a join exon b on a.exon_id = b.exon_id where a.transcript_id = '".$tr_id."' AND b.seq_region_id = '".$seq_region_id."'";

    $sth = $dbh_ENS->prepare($query);

    $sth->execute();

    $exons = $sth->fetchall_hashref('exon_id');

    

    #Return

    return($trans,$exons);



}



### get gene_id transcripts



sub get_gene_id_transcripts{



    #Catch

    my $dsn_ENS     =   $_[0];

    my $us_ENS      =   $_[1];

    my $pw_ENS      =   $_[2];

    my $gene_id     =   $_[3];

    

    #Init

    my $dbh_ENS = dbh($dsn_ENS,$us_ENS,$pw_ENS);

    my $trs = {};

    

    #Get transcripts

    my $query = "SELECT transcript_id,gene_id,seq_region_id,seq_region_strand,seq_region_start,seq_region_end,biotype,stable_id,canonical_translation_id from transcript where gene_id = $gene_id ";

    my $sth = $dbh_ENS->prepare($query);

    $sth->execute();

    $trs = $sth->fetchall_hashref('transcript_id');

    

    #Return

    return($trs);

    

}



### get genes from TIS-calling table ###



sub get_genes_per_chromosome {

    

    # Catch

    my $dbh         =   $_[0];

    my $id          =   $_[1];

    my $chr         =   $_[2];

    

    # Init

    my $genes = {};

    

    # Get genes

    my $query = "SELECT DISTINCT id as gene_id,biotype from TIS_sORFs_".$id." WHERE chr = '".$chr."' and annotation != 'intergenic'";

    my $sth = $dbh->prepare($query);

    $sth->execute();

    $genes = $sth->fetchall_hashref('gene_id');

    

    # Return

    return($genes);

}



### get genic sORF start positions from TIS-calling table ###



sub get_genic_sORF_starts_per_gene_id {

    

    # Catch

    my $dbh         =   $_[0];

    my $id          =   $_[1];

    my $gene_id     =   $_[2];

    

    # Init

    my $g_sORF_starts = {};

    

    # Get TISses

    my $query = "SELECT start,id as gene_id,biotype,chr,strand,annotation,start_codon,peak_shift,count from TIS_sORFs_".$id." WHERE id = '".$gene_id."' and annotation != 'intergenic'";

    my $sth = $dbh->prepare($query);

    $sth->execute();

    $g_sORF_starts = $sth->fetchall_hashref('start');

    

    # Return

    return($g_sORF_starts);

}



### get intergenic sORF start positions from TIS-calling table ###



sub get_intergenic_sORF_starts_per_chromosome {

    

    # Catch

    my $dbh         =   $_[0];

    my $id          =   $_[1];

    my $chr         =   $_[2];

    

    # Init

    my $ig_sORF_starts = {};

    

    # Get transcripts

    my $query = "SELECT id||'_'||start as intergene_start,id as intergene_id,chr,strand,start,annotation,upstream_gene_distance,downstream_gene_distance,start_codon,peak_shift,count from TIS_sORFs_".$id." WHERE chr = '".$chr."' and annotation = 'intergenic'";

    my $sth = $dbh->prepare($query);

    $sth->execute();

    $ig_sORF_starts = $sth->fetchall_hashref('intergene_start');

    

    # Return

    return($ig_sORF_starts);

}



### get snp information from SQLite into hash ###



sub get_SNPs_per_chromosome {

    

    # Catch

    my $dbh         = $_[0];

    my $chr         = $_[1];

    my $snp         = $_[2];

    my $CHX_lane    = $_[3];

    

    # Init

    my $SNPs = {};

    

    # Get SNPs based on chr

    my $query = "select chr,pos,ref,alt,af from snp_".$snp."_lane".$CHX_lane." where chr = '".$chr."'";

    

    my $sth = $dbh->prepare($query);

    $sth->execute();

    $SNPs = $sth->fetchall_hashref('pos');

    

    # Return

    return($SNPs);

}



### Fetch SNPs per sORF



sub fetch_SNPs_sORF {

    

    # Catch

    my $sorf_begin  =   $_[0];

    my $sorf_end    =   $_[1];

    my %SNPs        =   %{$_[2]};

    

    # Init

    my (@slice_list,%sORF_SNPs_all,%sORF_SNPs_def,$key,$value);

    

    # Take hash slice based on sORF_begin/sORF_end

    @slice_list = ($sorf_begin..$sorf_end);

    @sORF_SNPs_all{@slice_list} = @SNPs{@slice_list};

    foreach $key (keys %sORF_SNPs_all) {

        # Only print defined keys

        if ($sORF_SNPs_all{$key}) { $sORF_SNPs_def{$key} =  $sORF_SNPs_all{$key} }

    }

    

    #Return

    return(%sORF_SNPs_def);

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



### get sequence from binary chromosomes ###



sub get_sequence {

    

    #Catch

    my $chr   =     $_[0];

    my $start =     $_[1];

    my $end   =     $_[2];

    

    open(IN, "< ".$BIN_chrom_dir."/".$chr.".fa");

    binmode(IN);

    

    # Always forward strand sequence

    my $buffer; # Buffer to get binary data

    my $length = $end - $start + 1; # Length of string to read

    my $offset = $start-1; # Offset where to start reading

    

    seek(IN,$offset,0);

    read(IN,$buffer,$length);

    close(IN);

    

    # Return

    return($buffer);

}

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





### Get R



sub get_R{



    # Catch

    my $dsn                 =   $_[0];

    my $us                  =   $_[1];

    my $pw                  =   $_[2];

    my $analysis_id         =   $_[3];

    

    # Init

    my $dbh = dbh($dsn,$us,$pw);

    

    # Get R for analysis_id

    my $query = "select R from `TIS_sORFs_overview` where ID = '".$analysis_id."'";

    my $sth = $dbh->prepare($query);

    $sth->execute();

    my $R = $sth->fetch()->[0];

    

    # Return R

    return($R);

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

    

    $query = "select value from `arguments` where variable = \'nr_of_cores\'";

    $sth = $dbh->prepare($query);

    $sth->execute();

    my $nr_of_cores = $sth->fetch()->[0];

    

    $query = "select total from statistics where sample like \'%fastq1%\' and type = \'genomic\'";

    $sth = $dbh->prepare($query);

    $sth->execute();

    my $Mreads = $sth->fetch()->[0];

    $Mreads = $Mreads/1000000;

    

    $query = "select value from `arguments` where variable = \'mean_length_fastq1\'";

    $sth = $dbh->prepare($query);

    $sth->execute();

    my $mean_length_fastq1 = $sth->fetch()->[0];

    

    $query = "select value from `arguments` where variable = \'mean_length_fastq2\'";

    $sth = $dbh->prepare($query);

    $sth->execute();

    my $mean_length_fastq2 = $sth->fetch()->[0];

    

    # Return input variables

    return($ensemblversion,$species,$ens_db,$igenomes_root,$nr_of_cores,$Mreads,$mean_length_fastq1,$mean_length_fastq2);

    

}



### Get the analysis ids that need to be processed



sub get_analysis_ids {

    

    # Catch

    my $dsn     =   $_[0];

    my $us      =   $_[1];

    my $pw      =   $_[2];

    my $ids_in  =   $_[3]; #Either comma separated list of identifiers or "all"

    

    #Init

    my $dbh = dbh($dsn,$us,$pw);

    

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



### DNA to AA ###



sub DNAtoAA{

    

    #Catch

    my $triplet = $_[0];

    

    #Init

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

    

    my $AA = $AA1{uc($triplet)};

    

    #Return

    return($AA);

    

}



### get reverse complement sequence ###



sub revdnacomp {

    my $dna = shift;

    my $revcomp = reverse($dna);

    $revcomp =~ tr/ACGTacgt/TGCAtgca/;

    return $revcomp;

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

