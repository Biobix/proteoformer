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
# ./1_mapping.pl --name mESC --species mouse --ensembl 72 --cores 20 --readtype ribo --unique N --inputfile1 file1 --inputfile2 file2 --igenomes_root IGENOMES_ROOT (--mapper STAR --adaptor CTGTAGGCACCATCAAT --readlength 36 --truseq Y --firstrankmultimap N --out_bg_s_untr bg_s_untr --out_bg_as_untr bg_as_untr --out_bg_s_tr bg_s_tr --out_bg_as_tr bg_as_tr --out_sam_untr sam_untr --out_sam_tr sam_tr --out_sqlite sqliteDBName --work_dir getcwd --tmpfolder $TMP)

#For GALAXY
#1_mapping.pl --name "${experimentname}" --species "${organism}" --ensembl "${ensembl}" --cores "${cores}" --readtype $readtype.riboSinPair --unique "${unique}" --mapper "${mapper}" --readlength $readtype.readlength --adaptor $readtype.adaptor --inputfile1 $readtype.input_file1 --inputfile2 $readtype.input_file2 --out_bg_s_untr "${untreat_s_bg}"  --out_bg_as_untr "${untreat_as_bg}" --out_bg_s_tr "${treat_s_bg}" --out_bg_as_tr "${treat_as_bg}" --out_sam_untr "${untreat_sam}" --out_sam_tr "${treat_sam}" --out_sqlite "${out_sqlite}" --igenomes_root "${igenomes_root}"

# get the command line arguments
my ($work_dir,$run_name,$species,$ensemblversion,$cores,$mapper,$readlength,$readtype,$truseq,$tmpfolder,$adaptorSeq,$unique,$seqFileName1,$seqFileName2,$fastqName,$out_bg_s_untr,$out_bg_as_untr,$out_bg_s_tr,$out_bg_as_tr,$out_sam_untr,$out_sam_tr,$out_sqlite,$IGENOMES_ROOT,$ref_loc,$clipper,$phix,$rRNA,$snRNA,$tRNA,$tr_coord,$maxmultimap,$mismatch,$out_bam_tr_untr,$out_bam_tr_tr,$splicing,$FirstRankMultiMap);

GetOptions(
"inputfile1=s"=>\$seqFileName1,         	# the fastq file of the untreated data for RIBO-seq (no,CHX,EMT) or the 1st fastq for single/paired-end RNA-seq                  mandatory argument
"inputfile2=s"=>\$seqFileName2,         	# the fastq file of the treated data for RIBO-seq (PUR,LTM,HARR) or the 2nd fastq for paired-end RNA-seq                         mandatory argument
"name=s"=>\$run_name,                   	# Name of the run,                                                  			mandatory argument
"species=s"=>\$species,                 	# Species, eg mouse/human/fruitfly/arabidopsis                            mandatory argument
"ensembl=i"=>\$ensemblversion,          	# Ensembl annotation version, eg 66 (Feb2012),                      			mandatory argument
"cores=i"=>\$cores,                     	# Number of cores to use for Bowtie Mapping,                        			mandatory argument
"readtype=s"=>\$readtype,              		# The readtype (ribo, PE_polyA, SE_polyA, PE_total, SE_total)       			mandatory argument (default = ribo)
"mapper:s"=>\$mapper,                   	# The mapper used for alignment (STAR,TopHat2)       			                optional  argument (default = STAR)
"readlength:i"=>\$readlength,           	# The readlength (if RiboSeq take 50 bases),                        			optional  argument (default = 50)
"adaptor:s"=>\$adaptorSeq,              	# The adaptor sequence that needs to be clipped with fastx_clipper, 			optional  argument (default = CTGTAGGCACCATCAAT) => Ingolia paper (for ArtSeq = AGATCGGAAGAGCACAC)
"unique=s" =>\$unique,                  	# Retain the uniquely (and multiple) mapping reads (Y or N),        			mandatory argument
"tmp:s" =>\$tmpfolder,                  	# Folder where temporary files are stored,                          			optional  argument (default = $TMP or $CWD/tmp env setting)
"work_dir:s" =>\$work_dir,              	# Working directory ,                                               			optional  argument (default = $CWD env setting)
"out_bg_s_untr:s" =>\$out_bg_s_untr,    	# Output file for sense untreated count data (bedgraph)             			optional  argument (default = untreat_sense.bedgraph)
"out_bg_as_untr:s" =>\$out_bg_as_untr,  	# Output file for antisense untreated count data (bedgraph)         			optional  argument (default = untreat_antisense.bedgraph)
"out_bg_s_tr:s" =>\$out_bg_s_tr,        	# Output file for sense treated count data (bedgraph)               			optional  argument (default = treat_sense.bedgraph)
"out_bg_as_tr:s" =>\$out_bg_as_tr,      	# Output file for antisense treated count data (bedgraph)           			optional  argument (default = treat_antisense.bedgraph)
"out_sam_untr:s" =>\$out_sam_untr,      	# Output file for alignments of untreated data (sam)                			optional  argument (default = untreat.sam)
"out_sam_tr:s" =>\$out_sam_tr,         	  # Output file for alignments of treated data (sam)                  			optional  argument (default = treat.sam)
"out_bam_tr_untr:s" =>\$out_bam_tr_untr,	# Output file for alignments on transcript coordinates for untreated data (bam) optional  argument (default = untreat_tr.bam)
"out_bam_tr_tr:s" =>\$out_bam_tr_tr,   	  # Output file for alignments on transcript coordinates for treated data (bam)   optional  argument (default = treat_tr.bam)
"out_sqlite:s" =>\$out_sqlite,          	# sqlite DB output file                                             			optional  argument (default = results.db)
"igenomes_root=s" =>\$IGENOMES_ROOT,    	# IGENOMES ROOT FOLDER                                              			mandatory argument
"clipper:s" =>\$clipper,                	# what clipper is used (none or STAR or fastx)                     	      optional argument (default = none) or STAR or fastx
"phix:s" =>\$phix,                      	# map to phix DB prior to genomic mapping (Y or N)                  			optional argument (default = N)
"rRNA:s" =>\$rRNA,                      	# map to rRNA DB prior to genomic mapping (Y or N)                  			optional argument (default = Y)
"snRNA:s" =>\$snRNA,                    	# map to snRNA DB prior to genomic mapping (Y or N)                				optional argument (default = N)
"tRNA:s" =>\$tRNA,                      	# map to tRNA DB prior to genomic mapping (Y or N)                   			optional argument (default = N)
"tr_coord=s" =>\$tr_coord,					      # Generate alignment file based on transcript coordinates (Y or N)				optional argument (default = N)
"truseq=s" =>\$truseq,                    # If strands (+ and -) are assigned as in TruSeq or not (Y or N)          optional argument (default = Y)
"mismatch=i" =>\$mismatch,	              # Alignment will be output only if it has fewer mismatches than this value		optional argument (default = 2)
"maxmultimap=i" =>\$maxmultimap,		      # Alignments will be output only if the read maps fewer than this value		    optional argument (default = 16)
"splicing=s" =>\$splicing,                # Allow splicing for genome alignment for eukaryotic species (Y or N)         optional argument (default = Y)
"firstrankmultimap=s" =>\$FirstRankMultiMap  # Only retain the first ranked alignment of multimapper (Y or N)           optional argument (default = N)
);



###########################################################################
#Check all input variable and/or get default values and set extra variables
###########################################################################

my $CWD             = getcwd;
my $HOME            = $ENV{'HOME'};
$IGENOMES_ROOT      = ($ENV{'IGENOMES_ROOT'}) ? $ENV{'IGENOMES_ROOT'} : $IGENOMES_ROOT;
print "The following igenomes folder is used                    : $IGENOMES_ROOT\n";
my $TMP             = ($ENV{'TMP'}) ? $ENV{'TMP'} : ($tmpfolder) ? $tmpfolder : "$CWD/tmp"; # (1) get the TMP environment variable, (2) get the $tmpfolder variable, (3) get current_working_dir/tmp
print "The following tmpfolder is used                          : $TMP\n";

#Check if tmpfolder exists, if not create it...
if (!-d "$TMP") {
    system ("mkdir ". $TMP);
}

# comment on these
if ($work_dir){
    print "Working directory                                        : $work_dir\n";
} else {
    $work_dir = $CWD;
    print "Working directory                                        : $work_dir\n";
}
if ($seqFileName1){
    print "the fastq file of the untreated data for RIBO-seq (no,CHX,EMT) or the 1st fastq for single/paired-end RNA-seq                               : $seqFileName1\n";
} else {
    die "\nDon't forget to pass the FastQ file for untreated RIBO-seq or single or first paired-end RNA-seq using the --file or -f argument!\n\n";
}
if ($seqFileName2){
    print "the fastq file of the treated data for RIBO-seq (PUR,LTM,HARR) or the 2nd fastq for paired-end RNA-seq                                      : $seqFileName2\n";
} elsif (!defined($seqFileName2) && uc($readtype) eq 'RIBO') {
    die "\nDon't forget to pass the FastQ file for treated RIBO-seq (PUR,LTM,HARR) or 2nd fastq for paired-end RNA-seq using the --file or -f argument!\n\n";
}
if ($IGENOMES_ROOT){
    print "the igenomes_root folder used is                         : $IGENOMES_ROOT\n";
} else {
    die "\nDon't forget to pass the igenomes_root folder --igenomes_root or -ig argument!\n\n";
}
if ($run_name){
    print "Run name                                                 : $run_name\n";
} else {
    die "\nDon't forget to pass the Run Name using the --name or -n argument!\n\n";
}
if ($species){
    print "Species                                                  : $species\n";
} else {
    die "\nDon't forget to pass the Species name using the --species or -sp argument!\n\n";
}
if ($ensemblversion){
    print "Ensembl Version                                          : $ensemblversion\n";
} else {
    die "\nDon't forget to pass the Ensembl Version using the --ensembl or -ens argument!\n\n";
}
if ($cores){
    print "Number of cores to use for  mapping                      : $cores\n";
} else {
    die "\nDon't forget to pass number of cores to use for mapping using the --cores or -c argument!\n\n";
}
if ($readtype){
    print "readtype                                                 : $readtype\n";
} else {
    die "\nDon't forget to pass the read type using the --readtype or -r argument!\n\n";
}
if ($readlength){
    print "The readLength (for RiboSeq it should be set to 36)      : $readlength\n";
} else {
    #Choose default value for readlength
    $readlength = 36;
    print "The readLength (for RiboSeq it should be set to 36)      : $readlength\n";
}
if ($mismatch){
    print "The number of allowed mismatches during mapping          : $mismatch\n";
} else {
    #Choose default value for mismatch
    $mismatch = 2;
    print "The number of allowed mismatches during mapping          : $mismatch\n";
}
if ($mapper){
    print "The mapper used is                                       : $mapper\n";
} else {
    #Choose default value for mapper
    $mapper = "STAR";
    print "The mapper used is                                       : $mapper\n";
}
if ($unique){
    print "Unique mapped reads                                      : $unique\n";
} else {
    die "\nDon't forget to pass the unique or multiple read retention parameter --unique or -u argument!\n\n";
}
if ($maxmultimap){
    print "Maximun number of loci for reads to be acceptable        : $maxmultimap\n";
} else {
	$maxmultimap = 16;
    print "Maximun number of loci for reads to be acceptable        : $maxmultimap\n";
}
if ($FirstRankMultiMap){
    print "Only first ranked loci for multimap read is accepted     : $FirstRankMultiMap\n";
} else {
	$FirstRankMultiMap = 'N';
    print "Only first ranked loci for multimap read is accepted     : $FirstRankMultiMap\n";
}
if ($mismatch){
    print "Maximum number of mismatches to allow in an alignment    : $mismatch\n";
} else {
	$mismatch = 2;
    print "Maximum number of mismatches to allow in an alignment    : $mismatch\n";
}
if ($splicing){
    print "Splicing is set to                                       : $splicing\n";
} else {
    $splicing = 'Y';
    print "Splicing is set to                                       : $splicing\n";
}
if ($adaptorSeq){
    print "The adaptor sequence to be clipped with fastx_clipper/STAR    : $adaptorSeq\n";
} else {
    #Choose default value for AdaptorSeq
    $adaptorSeq = "CTGTAGGCACCATCAAT";
    print "The adaptor sequence to be clipped                       : $adaptorSeq\n";
}
if ($clipper){
    if (uc($clipper) eq 'STAR' && uc($mapper) eq 'TOPHAT2') {
        die "\nfastx adaptor clipper should be chosen when mapping with TopHat2! (--clipper fastx)\n\n";
    }
	print "The clipper used is                                      : $clipper\n";
} else {
		#Choose default value for clipper
        $clipper = 'none';
		print "The clipper used is                                      : $clipper\n";
}
if ($truseq){
    print "TruSeq strand assignment                                 : $truseq\n";
} else {
    $truseq = 'Y';
    print "TruSeq strand assignment                                 : $truseq\n";
}

if ($phix){
    print "Phix mapping prior to genomic                            : $phix\n";
} else {
    #Choose default value for phix
    $phix = 'N';
    print "Phix mapping prior to genomic                            : $phix\n";
}
if ($rRNA){
    print "rRNA mapping prior to genomic                            : $rRNA\n";
} else {
    #Choose default value for rRNA
    $rRNA = 'Y';
    print "rRNA mapping prior to genomic                            : $rRNA\n";
}
if ($snRNA){
    print "snRNA mapping prior to genomic                           : $snRNA\n";
} else {
    #Choose default value for snRNA
    $snRNA = 'N';
    print "snRNA mapping prior to genomic                           : $snRNA\n";
}
if ($tRNA){
    print "tRNA mapping prior to genomic                            : $tRNA\n";
} else {
    #Choose default value for snRNA
    $tRNA = 'N';
    print "tRNA mapping prior to genomic                            : $tRNA\n";
}
if ($tr_coord){
    print "generate alignment file based on transcriptome coord     : $tr_coord\n";
} else {
    #Choose default value for tr_coord
    $tr_coord = 'N';
    print "generate alignment file based on transcriptome coord     : $tr_coord\n";
}

# Create output directory
system "mkdir -p ".$work_dir."/output/";
system "mkdir -p ".$work_dir."/fastq/";

if (!defined($out_bg_s_untr))  		{$out_bg_s_untr     = $work_dir."/output/untreat_sense.bedgraph";}
if (!defined($out_bg_as_untr)) 		{$out_bg_as_untr    = $work_dir."/output/untreat_antisense.bedgraph";}
if (!defined($out_bg_s_tr))    		{$out_bg_s_tr       = $work_dir."/output/treat_sense.bedgraph";}
if (!defined($out_bg_as_tr))   		{$out_bg_as_tr      = $work_dir."/output/treat_antisense.bedgraph";}
if (!defined($out_sam_untr))   		{$out_sam_untr      = $work_dir."/".$mapper."/fastq1/untreat.sam";}
if (!defined($out_sam_tr))     		{$out_sam_tr        = $work_dir."/".$mapper."/fastq2/treat.sam";}
if (!defined($out_bam_tr_untr)) 	{$out_bam_tr_untr    = $work_dir."/".$mapper."/fastq1/untreat_tr.bam";}
if (!defined($out_bam_tr_tr))		{$out_bam_tr_tr     = $work_dir."/".$mapper."/fastq2/treat_tr.bam";}
if (!defined($out_sqlite))     		{$out_sqlite        = $work_dir."/SQLite/results.db";}

#ADDED FOR TEST ISSUES
my $clusterPhoenix = "N";
my $bowtie2Setting = "local"; # Or "end-to-end"
my $prev_mapper = (uc($mapper) eq "STAR") ? 'STAR' : "bowtie2"; # For STAR runs     => prev_mapper = "STAR",
                                                            # For TopHat2 runs  => prev_mapper = "bowtie2" (prev_mapper is used for Phix and rRNA mapping)

print "Phix/rRNA/tRNA/snRNA mapper                              : $prev_mapper\n";
#Set program run_name
my $run_name_short = $run_name;
$run_name = $run_name."_".$mapper."_".$unique."_".$ensemblversion;

#Conversion for species terminology
my $spec = (uc($species) eq "MOUSE") ? "Mus_musculus" : (uc($species) eq "HUMAN") ? "Homo_sapiens" : (uc($species) eq "ARABIDOPSIS") ? "Arabidopsis_thaliana" : (uc($species) eq "FRUITFLY") ? "Drosophila_melanogaster" : "";
my $spec_short = (uc($species) eq "MOUSE") ? "mmu" : (uc($species) eq "HUMAN") ? "hsa" : (uc($species) eq "ARABIDOPSIS") ? "ath" : (uc($species) eq "FRUITFLY") ? "dme" : "";
#Old mouse assembly = NCBIM37, new one is GRCm38. Old human assembly = GRCh37, the new one is GRCh38
my $assembly = (uc($species) eq "MOUSE" && $ensemblversion >= 70 ) ? "GRCm38"
: (uc($species) eq "MOUSE" && $ensemblversion < 70 ) ? "NCBIM37"
: (uc($species) eq "HUMAN" && $ensemblversion >= 76) ? "GRCh38"
: (uc($species) eq "HUMAN" && $ensemblversion < 76) ? "GRCh37"
: (uc($species) eq "ARABIDOPSIS") ? "TAIR10"
: (uc($species) eq "FRUITFLY" && $ensemblversion < 79) ? "BDGP5"
: (uc($species) eq "FRUITFLY" && $ensemblversion >= 79) ? "BDGP6" : "";

#Names for STAR/Bowtie2/Bowtie Indexes
#rRNA
my $IndexrRNA = $spec_short."_rRNA";
#snRNA
my $IndexsnRNA = $spec_short."_sn-o-RNA";
#tRNA
my $IndextRNA = $spec_short."_tRNA";
#Phix
my $IndexPhix = "phix";
#Genome
my $IndexGenome = "genome";
#cDNA
my $IndexCDNA = $spec.".".$assembly.".".$ensemblversion.".cdna.all";

#For STAR-Genome
my $readlengthMinus = $readlength - 1;
#my $ensemblversionforStar = ($species eq "mouse" && $ensemblversion >= 70) ? "70" : $ensemblversion;
my $STARIndexGenomeFolder = $spec_short.".".$assembly.".".$ensemblversion.".".$IndexGenome.".".$readlengthMinus."bpOverhang";


#Get executables
my ($bowtie_loc,$bowtie2_loc,$tophat2_loc,$STAR_loc,$sqlite_loc,$samtools_loc,$fastx_clip_loc,$fastx_trim_loc,$python_loc);
 $bowtie_loc = "bowtie";
 $bowtie2_loc = "bowtie2";
 $tophat2_loc = "tophat2";
 $STAR_loc = "STAR";
 $sqlite_loc = "sqlite3";
 $samtools_loc = "samtools";
 $fastx_clip_loc = "fastx_clipper";
 #$fastx_trim_loc = "fastx_trimmer";
 $fastx_trim_loc = "fastq_quality_trimmer";
 $python_loc = "python";

my $mapper_loc = (uc($mapper) eq "BOWTIE") ? $bowtie_loc
: (uc($mapper) eq "BOWTIE2") ? $bowtie2_loc
: (uc($mapper) eq "TOPHAT2") ? $tophat2_loc
: (uc($mapper) eq "STAR") ? $STAR_loc
: "" ;


#####Locations of rRNA databases and STAR indexes######
my $STAR_ref_loc =  $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Sequence/STARIndex/";
my $rRNA_fasta   =  $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Sequence/AbundantSequences/".$IndexrRNA.".fa";
my $snRNA_fasta  =  $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Sequence/AbundantSequences/".$IndexsnRNA.".fa";
my $phix_fasta   =  $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Sequence/AbundantSequences/".$IndexPhix.".fa";
my $cDNA_fasta   =  $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Sequence/AbundantSequences/".$IndexCDNA.".fa";
my $tRNA_fasta   =  $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Sequence/AbundantSequences/".$IndextRNA.".fa";
my $genome_fasta =  $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Sequence/WholeGenomeFasta/".$IndexGenome.".fa";

my $chromosome_sizes = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Annotation/Genes/ChromInfo.txt";
#print "chromosome_sizes_file=".$IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Annotation/Genes/ChromInfo.txt\n";

########################
#### Initiate results DB
########################

# Create SQLite DB for runName and create connection
my $db_sqlite_results  = $out_sqlite;

system "mkdir -p $work_dir/SQLite";
my $cmd = $sqlite_loc." ".$db_sqlite_results." \"PRAGMA auto_vacuum=1\"";
if (! -e $db_sqlite_results){system($cmd);}

# Sqlite Ensembl
my $dsn_sqlite_results = "DBI:SQLite:dbname=$db_sqlite_results";
my $us_sqlite_results  = "";
my $pw_sqlite_results  = "";

# Store input arguments
store_input_vars($dsn_sqlite_results,$us_sqlite_results,$pw_sqlite_results,$run_name_short,$ensemblversion,$species,$mapper,$unique,$adaptorSeq,$readlength,$readtype,$IGENOMES_ROOT,$cores);


############
## MAPPING
############
# Print/stati
print "Mapping sequences!\n";

my $stat_file;
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
    # Init statistics
    $stat_file=$run_name.".".$fastqName.".statistics.txt";
    system("touch ".$stat_file);

    if (uc($mapper) eq "BOWTIE") {
        map_bowtie($_,$fastqName,$mismatch);
    }
    elsif (uc($mapper) eq "BOWTIE2") {
        map_bowtie2($_,$fastqName);
    }
    elsif (uc($mapper) eq "TOPHAT2") {

        print "Mapping with TopHat2\n";
        my $start = time;
        if (uc($readtype) eq "RIBO") {
            map_topHat2($_,$fastqName,$clipper,$mismatch);
            RIBO_parse_store($_,$fastqName, 'Y'); # Only A-site parsing if RIBO-seq
            if ($unique eq "N" && $FirstRankMultiMap eq "N") {RIBO_parse_store($_,$fastqName, $unique)}
        }
        if (uc($readtype) eq "SE_POLYA") {
            map_topHat2($_,$fastqName,$clipper,$mismatch);
            RNA_parse_store($_,$fastqName, $unique, $truseq);
        }
        if (uc($readtype) =~ m/PE/) {
            map_topHat2($_,$fastqName,$clipper,$mismatch);
            RNA_parse_store($_,$fastqName, $unique, $truseq);
        }

        my $end = time - $start;
        printf("runtime TopHat against genomic: %02d:%02d:%02d\n\n",int($end/3600), int(($end % 3600)/60), int($end % 60));
    }
    elsif (uc($mapper) eq "STAR") {

		print "Mapping with STAR\n";
        my $start = time;
        if (uc($readtype) eq "RIBO") {
            map_STAR($_,$fastqName,$clipper,$mismatch);
            RIBO_parse_store($_,$fastqName, 'Y'); # Only A-site parsing if RIBO-seq
			if ($unique eq "N" && $FirstRankMultiMap eq "N") {RIBO_parse_store($_,$fastqName, $unique)}
        }
        if (uc($readtype) eq "SE_POLYA") {
            map_STAR($_,$fastqName,$clipper,$mismatch);
			RNA_parse_store($_,$fastqName, $unique, $truseq);
        }
        if (uc($readtype) =~ m/PE/) {
            map_STAR($_,$fastqName,$clipper,$mismatch);
            RNA_parse_store($_,$fastqName , $unique, $truseq);
        }

        my $end = time - $start;
        printf("runtime STAR against genomic: %02d:%02d:%02d\n\n",int($end/3600), int(($end % 3600)/60), int($end % 60));
    }
    # Store statistics in DB
    store_statistics($stat_file,$dsn_sqlite_results,$us_sqlite_results,$pw_sqlite_results);

}



############
# THE SUBS #
############

sub map_bowtie {

    # Catch
    my $seqFile = $_[0];
    my $seqFileName = $_[1];
    my $mismatch = $_[2];

    my $seed = 25;      # From Nature Protocol Paper Ingolia (for Bowtie1)

    # Print
    print "   BOWTIE used for mapping----------------------------------\n";
    print "   Mapping $seqFileName\n";

    # Prepare
    my $directory = $work_dir."/".$mapper."/".$seqFileName."/";
    system "mkdir -p $directory";
    my ($fasta,$command,$full_command);

    ############################################################################
    #Map to rrna first, cDNA second (non-unique), afterwards to genomic (unique)
    ############################################################################

    # Get input fastq file
    #$fasta = $work_dir."/fastq/$seqFileName".".fastq";
    $fasta = $seqFile;

    ###### Clip and Trim sequence
    print "     Clipping $seqFileName"."\n";

    my $adapter = $adaptorSeq;
    my $clip_command = $fastx_clip_loc." -Q33 -a ".$adapter." -l 20 -c -n –v -i ".$fasta." -o ".$work_dir."/fastq/$seqFileName"."_clip.fastq";

    system ($clip_command);
    print "     Trimming $seqFileName"."\n";
    system ($fastx_trim_loc." -Q33 -v -t 28 -l 26 -i ".$work_dir."/fastq/$seqFileName"."_clip.fastq -o ".$work_dir."/fastq/$seqFileName"."_trim.fastq");

    ###### Map against rRNA db
    # Check rRNA Index
    check_Bowtie_index($IndexrRNA,$mapper);
    # Build command
    $ref_loc = get_ref_loc($mapper);
    print "     Mapping against rRNA $seqFileName"."\n";
    $command = $mapper_loc." -p ".$cores." --sam --seedlen 23 --al ".$work_dir."/fastq/$seqFileName"."_rrna.fq --un ".$work_dir."/fastq/$seqFileName"."_norrna.fq ".$ref_loc."".$IndexrRNA." ".$work_dir."/fastq/$seqFileName"."_trim.fastq>/dev/null";
    # Run
    system($command);
    # Print
    print "   Finished rRNA mapping $seqFileName"." with seedlength 23\n";
    # Process statistics
    system("wc ".$work_dir."/fastq/$seqFileName"."_trim.fastq >> ".$stat_file);
    system("wc ".$work_dir."/fastq/$seqFileName"."_rrna.fq >> ".$stat_file);
    system("wc ".$work_dir."/fastq/$seqFileName"."_norrna.fq >> ".$stat_file);

    my $fasta_norrna = $work_dir."/fastq/$seqFileName"."_norrna.fq";


    ###### Map against cDNA
    # Check cDNA Index
    check_Bowtie_index($IndexCDNA,$mapper);
    # Build command
    print "     Mapping against cDNA $seqFileName"."\n";
    $command = "-l ".$seed." -n ".$mismatch." -p ".$cores." -m 255 --norc --best --strata --phred33-quals --quiet --al ".$directory."$seqFileName"."_".$seed."_cDNA_hit --un ".$directory."$seqFileName"."_".$seed."_cDNA_unhit --max ".$directory."$seqFileName"."_".$seed."_cDNA_max";
    $full_command = $mapper_loc." ".$command." ".$ref_loc."".$IndexCDNA." ".$fasta_norrna." > ".$directory."$seqFileName"."_".$seed."_cDNA_mapped";
    # Run
    system($full_command);
    # Print
    print "   Finished cDNA mapping $seqFileName"." with seedlength ".$seed."\n";
    # Process statistics
    system("wc ".$directory."$seqFileName"."_".$seed."_cDNA_hit  >> ".$stat_file);
    system("wc ".$directory."$seqFileName"."_".$seed."_cDNA_unhit >> ".$stat_file);

    my $fasta_nocDNA = $directory."$seqFileName"."_".$seed."_cDNA_unhit";

    ###### Map against genomic
    # Check genome Index
    #check_Bowtie_index($IndexGenome,$mapper);
    # Build command
    print "     Mapping against genomic $seqFileName"."\n";
    $command = "-l ".$seed." -n ".$mismatch." -p ".$cores." -m 1 --phred33-quals --quiet --al ".$directory."$seqFileName"."_".$seed."_genomic_hit --un ".$directory."$seqFileName"."_".$seed."_genomic_unhit";
    $full_command = $mapper_loc." ".$command." ".$ref_loc."".$IndexGenome." ".$fasta_nocDNA." > ".$directory."$seqFileName"."_".$seed."_genomic_mapped";
    # Run
    system($full_command);
    # Print
    print "   Finished genomic mapping $seqFileName"." with seedlength ".$seed."\n";
    #Process statistics
    system("wc ".$directory."$seqFileName"."_".$seed."_genomic_hit >> ".$stat_file);
    system("wc ".$directory."$seqFileName"."_".$seed."_genomic_unhit >> ".$stat_file);

}

sub map_bowtie2 {

    # Catch
    my $seqFile = $_[0];
    my $seqFileName = $_[1];

    # Print
    print "   BOWTIE2 used for mapping----------------------------------\n";
    print "   Mapping $seqFileName"."\n";

    # Prepare
    my $directory = $work_dir."/".$mapper.$bowtie2Setting."/".$seqFileName.""."/";
    system "mkdir -p $directory";
    my ($fasta,$command,$full_command);

    ############################################################################
    #Map to rrna first, cDNA second (non-unique), afterwards to genomic (unique)
    ############################################################################

    # Get input fastq file
    $fasta = $seqFile;

    my $fasta_to_rRNA;
    my $adapter = $adaptorSeq;

    ######  Clip and Trim sequence
    if ($bowtie2Setting eq 'local') {
        print "     Bowtie2 local version is selected, no trimming/clipping performed, feed file to rRNA mapping\n"; #Bowtie local allows mismatches at beginning or end of read
        $fasta_to_rRNA = $fasta;
    } elsif ($bowtie2Setting eq 'end-to-end') {
        print "     Clipping $seqFileName".", Bowtie2 end-to-end version selected\n";
        my $clip_command = $fastx_clip_loc." -Q33 -a ".$adapter." -l 20 -c -n –v -i ".$work_dir."/fastq/$seqFileName".".fastq -o ".$work_dir."/fastq/$seqFileName"."_clip.fastq";
        system ($clip_command);
        print "     Trimming $seqFileName".", Bowtie2 end-to-end version selected\n";
        system ($fastx_trim_loc." -Q33 -v -t 28 -l 26 -i ".$work_dir."/fastq/$seqFileName"."_clip.fastq -o ".$work_dir."/fastq/$seqFileName"."_trim.fastq");
        $fasta_to_rRNA = $work_dir."/fastq/$seqFileName"."_trim.fastq";
    }


    ###### Map against rRNA db
    # Check rRNA Index
    check_Bowtie_index($IndexrRNA,$mapper);
    # Build command
    print "     Mapping against rRNA $seqFileName"."\n";
    $ref_loc = get_ref_loc($mapper);
    $command = $mapper_loc." -p ".$cores." --".$bowtie2Setting." --sensitive-local --al ".$work_dir."/fastq/$seqFileName"."_rrna.fq --un ".$work_dir."/fastq/$seqFileName"."_norrna.fq  -x ".$ref_loc."".$IndexrRNA." -U ".$fasta_to_rRNA.">/dev/null";
    # Run
    system($command);
    # Print
    print "   Finished rRNA multiseed alignment of $seqFileName"." with seedlength 20\n";
    # Process statistics
    system("wc ".$fasta_to_rRNA." >> ".$stat_file);
    system("wc ".$work_dir."/fastq/$seqFileName"."_rrna.fq >> ".$stat_file);
    system("wc ".$work_dir."/fastq/$seqFileName"."_norrna.fq >> ".$stat_file);

    my $fasta_norrna = $work_dir."/fastq/$seqFileName"."_norrna.fq";


    ###### Map against cDNA
    # Check cDNA Index
    check_Bowtie_index($IndexCDNA,$mapper);
    # Build command
    print "     Mapping against cDNA $seqFileName"."\n";
    $command = " -p ".$cores." --norc --".$bowtie2Setting." --phred33 --quiet --al ".$directory."$seqFileName"."_20_cDNA_hit --un ".$directory."$seqFileName"."_20_cDNA_unhit ";
    $full_command = $mapper_loc." ".$command." -x ".$ref_loc."".$IndexCDNA." -U ".$fasta_norrna." -S ".$directory."$seqFileName"."_20_cDNA_mapped";
    # Run
    system($full_command);
    # Print
    print "   Finished cDNA mapping $seqFileName"." with seedlength 20"."\n";
    # Process statistics
    system("wc ".$directory."$seqFileName"."_20_cDNA_hit  >> ".$stat_file);
    system("wc ".$directory."$seqFileName"."_20_cDNA_unhit >> ".$stat_file);

    my $fasta_nocDNA = $directory."$seqFileName"."_20_cDNA_unhit";

    ###### Map against genomic
    # Check genome Index
    #check_Bowtie_index($IndexGenome,$mapper);
    # Build command
    print "     Mapping against genomic $seqFileName"."\n";
    $command = "-p ".$cores." --".$bowtie2Setting." --phred33 --quiet --al ".$directory."$seqFileName"."_20_genomic_hit --un ".$directory."$seqFileName"."_20_genomic_unhit";
    $full_command = $mapper_loc." ".$command." -x ".$ref_loc."".$IndexGenome." -U  ".$fasta_nocDNA." -S ".$directory."$seqFileName"."_20_genomic_mapped";
    # Run
    system($full_command);
    # Print
    print "   Finished genomic mapping $seqFileName"." with seedlength 20"."\n";
    #Process statistics
    system("wc ".$directory."$seqFileName"."_20_genomic_hit >> ".$stat_file);
    system("wc ".$directory."$seqFileName"."_20_genomic_unhit >> ".$stat_file);

}

sub map_STAR {

    # Catch
    my $seqFiles = $_[0];
    my $seqFileName = $_[1];
    my $clipper = $_[2];
    my $mismatch = $_[3];
    my @splitSeqFiles = split(/,/,$seqFiles);
    my $seqFile  = $splitSeqFiles[0];
    my $seqFile2 = $splitSeqFiles[1];

    my ($fasta,$fasta1,$fasta2,$fasta_norrna,$command,$full_command,$outSAMmultNmax);
    my $ref_loc = get_ref_loc($mapper);
    # Check if main STAR index folder exists
    if (!-d $ref_loc) { system("mkdir -p ".$ref_loc); print "main STAR folder has been created\n";}

    # Print
    print "   STAR used for mapping----------------------------------\n";
    print "   Mapping $seqFileName"."\n";

    # Get input fastq file
    if (uc($readtype) eq "RIBO" || uc($readtype) =~ m/SE/) {
        $fasta = $seqFile;
    } elsif (uc($readtype) =~ m/PE/) {
        $fasta1 = $seqFile;
        $fasta2 = $seqFile2;
    }

    # Prepare
    my $directory = $work_dir."/".$mapper."/".$seqFileName.""."/";
    system "mkdir -p $directory";

    # We need to first clip the adapter with fastx_clipper?
    if (uc($clipper) eq "FASTX") {

        print "     Clipping $seqFileName"." using fastx_clipper tool\n";
        # Without length cut-off and adaptor presence
        my $clip_command = $fastx_clip_loc." -Q33 -a ".$adaptorSeq." -l 26 -n –v -i ".$fasta." -o ".$work_dir."/fastq/".$seqFileName."_clip.fastq";
        print "     ".$clip_command."\n";
        system ($clip_command);
        $fasta = $work_dir."/fastq/$seqFileName"."_clip.fastq";
        print "clipfasta= $fasta\n";

	    print "     Trimming $seqFileName"." using fastq_quality_trimmer tool\n";
		my $trim_command = $fastx_trim_loc." -Q33 -v -t 28 -l 26  -i ".$fasta." -o ".$work_dir."/fastq/".$seqFileName."_clip_trim.fastq";
	    system ($trim_command);
		$fasta = $work_dir."/fastq/$seqFileName"."_clip_trim.fastq";
		print "trimfasta= $fasta\n";

    }

    my $clip_stat = (uc($clipper) eq "FASTX" || uc($clipper) eq "NONE") ? " " : "--clip3pAdapterSeq ".$adaptorSeq." --clip3pAdapterMMp 0.1 ";

    #GO FOR STAR-pHIX mapping
    if ($phix eq "Y") {
        #Check if pHIX STAR index exists?
        check_STAR_index($IndexPhix);

        print "     Mapping against Phix $seqFileName"."\n";
        ###### Map against phix db using STAR

        # Build command
        system("mkdir  -p ".$work_dir."/fastq/nophix");
        $command = $STAR_loc." --genomeLoad NoSharedMemory --outFilterMismatchNmax 2  --seedSearchStartLmaxOverLread .5 ".$clip_stat." --alignIntronMax 1 --genomeDir ".$ref_loc.$IndexPhix." --readFilesIn ".$fasta."  --outFileNamePrefix ".$work_dir."/fastq/nophix/ --runThreadN ".$cores." --outReadsUnmapped Fastx";
        # Run
        system($command);
        # Print
        print "   Finished rRNA multiseed mapping $seqFileName"."\n";
        # Process statistics
        open (STATS, ">>".$stat_file);
        my ($inReads, $mappedReadsU,$mappedReadsM, $unmappedReads) = parseLogSTAR($work_dir."/fastq/nophix/");
        my $mappedReads = $mappedReadsU + $mappedReadsM;
        print STATS "STAR ".$run_name." ".$seqFileName." "."fastq phix ".$inReads."\n";
        print STATS "STAR ".$run_name." ".$seqFileName." "."hit phix ".$mappedReads."\n";
        print STATS "STAR ".$run_name." ".$seqFileName." "."unhit phix ".$unmappedReads."\n";
        close(STATS);

        # Rename unmapped.out.mate1 for further analysis
        system("mv ".$work_dir."/fastq/nophix/Unmapped.out.mate1 ".$work_dir."/fastq/".$seqFileName."_nophix.fq"); #For further processing against genomic!!
        $fasta = $work_dir."/fastq/".$seqFileName."_nophix.fq";

    }

    # If RIBO-SEQ: prior rRNA and/or snRNA and/or trRNA mapping is NECESSARY
    # If RNA-SEQ:  prior filtering on rRNA snRNA and/or trRNA is OPTIONAL
    #if (uc($readtype) eq "RIBO") {
        #GO FOR rRNA mapping
        if ($rRNA eq "Y") {
			if ($prev_mapper eq 'STAR') {

				#Check if rRNA STAR index exists?
				check_STAR_index($IndexrRNA);

				print "     Mapping against rRNA $seqFileName"."\n";
				###### Map against rRNA db using STAR (includes adapter clipping or exluding adapter clipping, dependent on previous clipping process)

				$ref_loc = get_ref_loc($prev_mapper);
				# Build command
				$command = $STAR_loc." --genomeLoad NoSharedMemory --seedSearchStartLmaxOverLread .5 ".$clip_stat." --genomeDir ".$ref_loc.$IndexrRNA." --readFilesIn ".$fasta." --outFilterMultimapNmax 1000 --outFilterMismatchNmax 2 --outFileNamePrefix ".$work_dir."/fastq/ --runThreadN ".$cores." --outReadsUnmapped Fastx";

				#print "     $command\n";
				# Run
				system($command);
				# Print
				print "   Finished rRNA multiseed mapping $seqFileName"."\n";
				# Process statistics
				open (STATS, ">>".$stat_file);
				my ($inReads, $mappedReadsU,$mappedReadsM, $unmappedReads) = parseLogSTAR($work_dir."/fastq/");
				my $mappedReads = $mappedReadsU + $mappedReadsM;
				print STATS "STAR ".$run_name." ".$seqFileName." "."fastq rRNA ".$inReads."\n";
				print STATS "STAR ".$run_name." ".$seqFileName." "."hit rRNA ".$mappedReads."\n";
				print STATS "STAR ".$run_name." ".$seqFileName." "."unhit rRNA ".$unmappedReads."\n";
				close(STATS);

				# Rename unmapped.out.mate1 for further analysis
				system("mv ".$work_dir."/fastq/Unmapped.out.mate1 ".$work_dir."/fastq/".$seqFileName."_norrna.fq"); #For further processing against genomic!!

			}
			#GO FOR CLIP-TRIM-BOWTIE to get clipped, trimmed, norRNA reads
			elsif ($prev_mapper eq 'bowtie') {
				###### Clip and Trim sequence
				print "     Clipping $seqFileName"."\n";

				my $adapter = $adaptorSeq;
				my $clip_command = $fastx_clip_loc." -Q33 -a ".$adapter." -l 20 -c -n –v -i ".$work_dir."/fastq/$seqFileName".".fastq -o ".$work_dir."/fastq/$seqFileName"."_clip.fastq";
				#print "     ".$clip_command."\n";
				system ($clip_command);
				print "     Trimming $seqFileName"."\n";
				system ($fastx_trim_loc." -Q33 -v -t 28 -l 26 -i ".$work_dir."/fastq/$seqFileName"."_clip.fastq -o ".$work_dir."/fastq/$seqFileName"."_trim.fastq");

				###### Map against rRNA db using Bowtie
				# Check rRNA Index
				check_Bowtie_index($IndexrRNA,$prev_mapper);
				# Build command
				print "     Mapping against rRNA $seqFileName"."\n";
				$ref_loc = get_ref_loc($prev_mapper);
				$command = $bowtie_loc." -p ".$cores." --sam --seedlen 23 --al ".$work_dir."/fastq/$seqFileName"."_rrna.fq --un ".$work_dir."/fastq/$seqFileName"."_norrna.fq ".$ref_loc."".$IndexrRNA." ".$work_dir."/fastq/$seqFileName"."_trim.fastq>/dev/null";
				# Run
				system($command);
				# Print
				print "   Finished rRNA mapping $seqFileName"." with seedlength 23\n";
				# Process statistics
				system("wc ".$work_dir."/fastq/$seqFileName"."_trim.fastq >> ".$stat_file);
				system("wc ".$work_dir."/fastq/$seqFileName"."_rrna.fq >> ".$stat_file);
				system("wc ".$work_dir."/fastq/$seqFileName"."_norrna.fq >> ".$stat_file);

			}
			# GO FOR BOWTIE2_LOCAL_VERY-FAST to get unclipped/untrimmed norRNA reads
			elsif ($prev_mapper eq 'bowtie2') {
				print "     Bowtie2 local version is selected, no trimming performed, feed clipped file to rRNA mapping\n"; #Bowtie local allows mismatches at beginning or end of read
				###### Map against rRNA db
				# Check rRNA Index
				check_Bowtie_index($IndexrRNA,$prev_mapper);
				# Build command
				$ref_loc = get_ref_loc($prev_mapper);
				print "     Mapping against rRNA $seqFileName"."\n";
				$command = $bowtie2_loc." -p ".$cores." --".$bowtie2Setting." --sensitive-local --al ".$work_dir."/fastq/$seqFileName"."_rrna.fq --un ".$work_dir."/fastq/$seqFileName"."_norrna.fq  -x ".$ref_loc."".$IndexrRNA." -U ".$fasta.">/dev/null";
				# Run
				system($command);
				# Print
				print "   Finished rRNA multiseed alignment of $seqFileName"." with seedlength 20\n";
				# Process statistics
				system("wc ".$fasta." >> ".$stat_file);
				system("wc ".$work_dir."/fastq/".$seqFileName."_rrna.fq >> ".$stat_file);
				system("wc ".$work_dir."/fastq/".$seqFileName."_norrna.fq >> ".$stat_file);

			}
		$fasta = $work_dir."/fastq/$seqFileName"."_norrna.fq";
        }

        #GO FOR snRNA mapping
        if ($snRNA eq "Y") {
            if ($prev_mapper eq 'STAR') {

                #Check if snRNA STAR index exists?
                check_STAR_index($IndexsnRNA);

                print "     Mapping against snRNA $seqFileName"."\n";
                ###### Map against snRNA db using STAR (includes adapter clipping or exluding adapter clipping, dependent on previous clipping process)

                $ref_loc = get_ref_loc($prev_mapper);
                # Build command
                $command = $STAR_loc." --genomeLoad NoSharedMemory --seedSearchStartLmaxOverLread .5 ".$clip_stat." --genomeDir ".$ref_loc.$IndexsnRNA." --readFilesIn ".$fasta." --outFilterMultimapNmax 1000 --outFilterMismatchNmax 2 --outFileNamePrefix ".$work_dir."/fastq/ --runThreadN ".$cores." --outReadsUnmapped Fastx";

                #print "     $command\n";
                # Run
                system($command);
                # Print
                print "   Finished snRNA multiseed mapping $seqFileName"."\n";
                # Process statistics
                open (STATS, ">>".$stat_file);
                my ($inReads, $mappedReadsU,$mappedReadsM, $unmappedReads) = parseLogSTAR($work_dir."/fastq/");
                my $mappedReads = $mappedReadsU + $mappedReadsM;
                print STATS "STAR ".$run_name." ".$seqFileName." "."fastq snRNA ".$inReads."\n";
                print STATS "STAR ".$run_name." ".$seqFileName." "."hit snRNA ".$mappedReads."\n";
                print STATS "STAR ".$run_name." ".$seqFileName." "."unhit snRNA ".$unmappedReads."\n";
                close(STATS);

                # Rename unmapped.out.mate1 for further analysis
                system("mv ".$work_dir."/fastq/Unmapped.out.mate1 ".$work_dir."/fastq/".$seqFileName."_norrna_nosnrna.fq"); #For further processing against genomic!!
            }
        $fasta = $work_dir."/fastq/$seqFileName"."_norrna_nosnrna.fq";
        }

        #GO FOR tRNA mapping
        if ($tRNA eq "Y") {
            if ($prev_mapper eq 'STAR') {

                #Check if tRNA STAR index exists?
                check_STAR_index($IndextRNA);

                print "     Mapping against tRNA $seqFileName"."\n";
                ###### Map against tRNA db using STAR (includes adapter clipping or exluding adapter clipping, dependent on previous clipping process)

                $ref_loc = get_ref_loc($prev_mapper);
                # Build command
                $command = $STAR_loc." --genomeLoad NoSharedMemory --seedSearchStartLmaxOverLread .5 ".$clip_stat." --genomeDir ".$ref_loc.$IndextRNA." --readFilesIn ".$fasta." --outFilterMultimapNmax 1000 --outFilterMismatchNmax 2 --outFileNamePrefix ".$work_dir."/fastq/ --runThreadN ".$cores." --outReadsUnmapped Fastx";

                #print "     $command\n";
                # Run
                system($command);
                # Print
                print "   Finished tRNA multiseed mapping $seqFileName"."\n";
                # Process statistics
                open (STATS, ">>".$stat_file);
                my ($inReads, $mappedReadsU,$mappedReadsM, $unmappedReads) = parseLogSTAR($work_dir."/fastq/");
                my $mappedReads = $mappedReadsU + $mappedReadsM;
                print STATS "STAR ".$run_name." ".$seqFileName." "."fastq tRNA ".$inReads."\n";
                print STATS "STAR ".$run_name." ".$seqFileName." "."hit tRNA ".$mappedReads."\n";
                print STATS "STAR ".$run_name." ".$seqFileName." "."unhit tRNA ".$unmappedReads."\n";
                close(STATS);

                # Rename unmapped.out.mate1 for further analysis
                system("mv ".$work_dir."/fastq/Unmapped.out.mate1 ".$work_dir."/fastq/".$seqFileName."_norrna_nosnrna_notrna.fq"); #For further processing against genomic!!

            }
        $fasta = $work_dir."/fastq/$seqFileName"."_norrna_nosnrna_notrna.fq";
        }
        #}

    # Check genomic STAR index
    # If it doesn't exist, it's created from data within the iGenome directory
    check_STAR_index($STARIndexGenomeFolder);

    # Map vs genomic
    print "     Mapping against genomic $seqFileName"."\n";

    if (uc($readtype) eq "RIBO") {

		    if ($unique eq 'Y') {$maxmultimap = 1;}	# to ensure unique mapping
        $ref_loc = get_ref_loc($mapper);
        $command = $STAR_loc." --outSAMattributes Standard --outSAMtype BAM SortedByCoordinate --genomeLoad NoSharedMemory ".$clip_stat." --seedSearchStartLmaxOverLread .5 --outFilterMultimapNmax ".$maxmultimap." --outMultimapperOrder Random --outFilterMismatchNmax ".$mismatch." --genomeDir ".$ref_loc.$STARIndexGenomeFolder." --runThreadN ".$cores." --outFileNamePrefix ".$directory." --readFilesIn ".$fasta;
		    if (uc($tr_coord) eq 'Y') {$command .= " --quantMode TranscriptomeSAM "}		# Output transcriptome coordinate to bam file
        if (uc($splicing) eq 'N') {$command .= " --alignIntronMax 1 --alignIntronMin 2 "}        # Don't allow splicing if splicing is set to 'N'
        if ($FirstRankMultiMap eq 'Y') {$command .= " --outSAMmultNmax 1 "}   # To ensure to output only first hit to SAM (either the best scoring one or a random chosen one from a list of equally scoring hits)
    } elsif (uc($readtype) =~ m/PE/) {

		    if ($unique eq 'Y') {$maxmultimap = 1;}	# to ensure unique mapping
        $command = $STAR_loc." --outSAMattributes Standard --outSAMtype BAM SortedByCoordinate --genomeLoad NoSharedMemory ".$clip_stat." --seedSearchStartLmaxOverLread .5 --outFilterMultimapNmax ".$maxmultimap." --outMultimapperOrder Random --outFilterMismatchNmax ".$mismatch." --genomeDir ".$ref_loc.$STARIndexGenomeFolder." --runThreadN ".$cores." --outFileNamePrefix ".$directory." --readFilesIn ".$fasta1." ".$fasta2." --chimSegmentMin 15 --chimJunctionOverhangMin 15 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMunmapped Within --outReadsUnmapped Fastx --outFilterType BySJout";
		    if (uc($tr_coord) eq 'Y') {$command .= " --quantMode TranscriptomeSAM "}		# Add transcript coordinate to bam output
        if (uc($splicing) eq 'N') {$command .= " --alignIntronMax 1 --alignIntronMin 2 "}        # Don't allow splicing if splicing is set to 'N'
        if ($FirstRankMultiMap eq 'Y') {$command .= " --outSAMmultNmax 1 "}   # To ensure to output only first hit to SAM (either the best scoring one or a random chosen one from a list of equally scoring hits)
    } elsif (uc($readtype) =~ m/SE/) {

		    if ($unique eq 'Y') {$maxmultimap = 1;}	# to ensure unique mapping
        $command = $STAR_loc." --outSAMattributes Standard --outSAMtype BAM SortedByCoordinate --genomeLoad NoSharedMemory ".$clip_stat." --seedSearchStartLmaxOverLread .5 --outFilterMultimapNmax ".$maxmultimap." --outMultimapperOrder Random --outFilterMismatchNmax ".$mismatch." --genomeDir ".$ref_loc.$STARIndexGenomeFolder." --runThreadN ".$cores." --outFileNamePrefix ".$directory." --readFilesIn ".$fasta;
		    if (uc($tr_coord) eq 'Y') {$command .= " --quantMode TranscriptomeSAM "}		# Add transcript coordinate to bam output
        if (uc($splicing) eq 'N') {$command .= " --alignIntronMax 1 --alignIntronMin 2 "}        # Don't allow splicing if splicing is set to 'N'
        if ($FirstRankMultiMap eq 'Y') {$command .= " --outSAMmultNmax 1 "}   # To ensure to output only first hit to SAM (either the best scoring one or a random chosen one from a list of equally scoring hits)
    }


    print "     ".$command."\n";
    system($command);
    systemError("STAR",$?,$!);

    # Samtools sam_to_bam and sort not necessary, STAR already includes this with "--outSAMtype BAM SortedByCoordinate" option
    # # convert SAM output file to BAM file
    # print "converting SAM output to BAM...\n";
    # system($samtools_loc." view -bS -o ".$directory."Aligned.out.bam ".$directory."Aligned.out.sam > /dev/null 2>&1");
    # systemError("Samtools view",$?,$!);
    # # sort the BAM file
    # print "sorting STAR hits...\n";
    # system($samtools_loc." sort -@ ". $cores. " -m 1000M ".$directory."Aligned.out.bam ".$directory."Aligned.sorted 2>&1" );
    # systemError("Samtools sort",$?,$!);

    # #  convert BAM back to SAM file
     print "converting BAM back to SAM...\n";
     system($samtools_loc." view -h -o ".$directory."Aligned.sorted.sam ".$directory."Aligned.sortedByCoord.out.bam > /dev/null 2>&1");
     systemError("Samtools view",$?,$!);

    # rename SAM output file
    print "renaming SAM output file...\n";

    # Bam file depends on what fastq file is processed (fastq1 = untreated, fastq2 = treaeted; that is for RIBO-seq experiments)
    my $bamf = ($seqFileName  eq 'fastq1') ? $out_sam_untr : $out_sam_tr;
	system("mv ".$directory."Aligned.sorted.sam ".$bamf);

	# remove redundant files
	system("rm ".$directory."Aligned.out.bam > /dev/null 2>&1");
	system("rm ".$directory."Aligned.out.sam > /dev/null 2>&1");
	system("rm ".$directory."Aligned.sorted.bam > /dev/null 2>&1");

	# Handle transcript coordinates if required
	if (uc($tr_coord) eq "Y" ) {

		print "sorting STAR for transcripts coordinates hits...\n";
		system($samtools_loc." sort -@ ". $cores. " -m 1000M ".$directory."Aligned.toTranscriptome.out.bam ".$directory."Aligned.toTranscriptome.out.sorted 2>&1" );
		systemError("Samtools sort",$?,$!);

		#  convert BAM back to SAM file
		#print "converting BAM to SAM...\n";
		#system($samtools_loc." view -h -o ".$directory."Aligned.toTranscriptome.out.sorted.sam ".$directory."Aligned.toTranscriptome.out.sorted.bam > /dev/null 2>&1");
		#systemError("Samtools view",$?,$!);

		# rename SAM output file
		#print "renaming SAM output file...\n";

		# Bam file depends on what fastq file is processed (fastq1 = untreated, fastq2 = treated; that is for RIBO-seq experiments)
		my $bamf = ($seqFileName  eq 'fastq1') ? $out_bam_tr_untr : $out_bam_tr_tr;
		system("mv ".$directory."Aligned.toTranscriptome.out.sorted.bam ".$bamf);

		print "transcript coordinate $bamf\n";
		system("rm ".$directory."Aligned.toTranscriptome.out.bam > /dev/null 2>&1");
		#system("rm ".$directory."Aligned.toTranscriptome.out.sorted.bam");
		#system("rm ".$directory."Aligned.toTranscriptome.out.sorted.sam");
	}


    # Process statistics
    open (STATS, ">>".$stat_file);
    my ($inReads,$mappedReadsU,$mappedReadsM, $unmappedReads) = parseLogSTAR($work_dir."/".$mapper."/".$seqFileName."/");
    print STATS "STAR ".$run_name." ".$seqFileName." "."fastq genomic ".$inReads."\n";
    print STATS "STAR ".$run_name." ".$seqFileName ." "."hitU genomic ".$mappedReadsU."\n";
    print STATS "STAR ".$run_name." ".$seqFileName ." "."hitM genomic ".$mappedReadsM."\n";
    print STATS "STAR ".$run_name." ".$seqFileName." "."unhit genomic ".$unmappedReads."\n";
    close(STATS);

}

sub check_STAR_index {

    # Catch
    my $starIndexDir = $_[0];

    $ref_loc = get_ref_loc($mapper);
    my $starIndexDirComplete = $ref_loc.$starIndexDir;

    #print "$starIndexDirComplete\n";
    print "     ----checking for STAR index folder...\n";
    if (!-d $starIndexDirComplete){
        if ($starIndexDir =~ /rRNA/ || $starIndexDir =~ /sn-o-RNA/ || $starIndexDir =~ /tRNA/ ) {
            print "no STAR directory ". $starIndexDir ." found\ncreating STAR index without annotation ...\n";
            system("mkdir -p ".$starIndexDirComplete);

            #Get rRNA fasta and from iGenome folder (rRNA sequence are located in /Sequence/AbundantSequences folder)
            #GenomeSAIndexNbases needs to be adapted because small "rRNA-genome" size +/- 8000bp. => [GenomeSAIndexNbases = log2(size)/2 - 1]
            my $PATH_TO_FASTA = ($starIndexDir =~ /rRNA/) ? $rRNA_fasta : ($starIndexDir =~ /tRNA/) ? $tRNA_fasta : $snRNA_fasta;
            system($STAR_loc." --runMode genomeGenerate --genomeSAindexNbases 6  --genomeDir ".$starIndexDirComplete." --runThreadN ". $cores." --genomeFastaFiles ".$PATH_TO_FASTA);
            print $STAR_loc." --runMode genomeGenerate --genomeSAindexNbases 6  --genomeDir ".$starIndexDirComplete." --runThreadN ". $cores." --genomeFastaFiles ".$PATH_TO_FASTA. "\n";
            systemError("STAR genome",$?,$!);

        }
        elsif ($starIndexDir =~ /phix/) {
            print "no STAR directory". $starIndexDir ."found\ncreating STAR index without annotation ...\n";
            system("mkdir -p ".$starIndexDirComplete);

            #Get phix fasta and from iGenome folder (phix sequence(s) are located in /Sequence/AbundantSequences folder)
            #GenomeSAIndexNbases needs to be adapted because small "rRNA-genome" size +/- 8000bp. => [GenomeSAIndexNbases = log2(size)/2 - 1]
            my $PATH_TO_FASTA = $phix_fasta;
            system($STAR_loc." --runMode genomeGenerate --genomeSAindexNbases 6  --genomeDir ".$starIndexDirComplete." --runThreadN ". $cores." --genomeFastaFiles ".$PATH_TO_FASTA);
            systemError("STAR genome",$?,$!);

        }
        elsif ($starIndexDir =~ /genome/) {
            print "no STAR ".$spec." genome directory found\ncreating STAR index with annotation ...\n";
            system("mkdir -p ".$starIndexDirComplete);

            #Get Genome fasta and GTF file from iGenome folder (genome sequence is located in /Sequence/WholeGenomeFasta folder)
            my $PATH_TO_FASTA = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Sequence/WholeGenomeFasta/genome.fa";
            my $PATH_TO_GTF   = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Annotation/Genes/genes_".$ensemblversion".gtf";
            system($STAR_loc." --runMode genomeGenerate --sjdbOverhang ".$readlengthMinus." --genomeDir ".$starIndexDirComplete." --genomeFastaFiles ".$PATH_TO_FASTA." --sjdbGTFfile ".$PATH_TO_GTF." --runThreadN ".$cores);
            systemError("STAR genome",$?,$!);
        }
    }

}

sub check_Bowtie_index {

    # Catch
    my $bowtieIndex = $_[0];
    my $mapper      = $_[1];

    # Set
    $ref_loc = get_ref_loc($mapper);
    if (!-d $ref_loc){ system("mkdir -p ".$ref_loc); }

    my $bowtieIndexFile = $ref_loc.$bowtieIndex;
    my $bowtieIndexFileCheck = ($mapper eq "bowtie") ? $ref_loc.$bowtieIndex.".1.ebwt" : $ref_loc.$bowtieIndex.".1.bt2";
    my $bowtieloc = ($mapper eq "bowtie") ? $bowtie_loc : $bowtie2_loc;

    # Check if index exists
    print "     ----checking for Bowtie index file...\n";
    print "         $bowtieIndexFileCheck\n";
    if (-e $bowtieIndexFileCheck){
        print "     ----ok, file found\n\n";
    } else {
        if ($bowtieIndex =~ /rRNA/) {
            print "no Bowtie ". $bowtieIndexFile ." found\ncreating bowtie index ...\n";

            system("mkdir -p ".$bowtieIndexFile);

            #Get rRNA fasta from iGenome folder (rRNA sequence are located in /Sequence/AbundantSequences folder)
            my $PATH_TO_FASTA = $rRNA_fasta;
            print $bowtieloc."-build ".$PATH_TO_FASTA." ".$bowtieIndexFile."\n";
            system($bowtieloc."-build ".$PATH_TO_FASTA." ".$bowtieIndexFile);
            systemError("Bowtie index creation error",$?,$!);
            print "done!\n\n";

        }
        if ($bowtieIndex =~ /tRNA/) {
            print "no Bowtie ". $bowtieIndexFile ." found\ncreating bowtie index ...\n";

            #Get tRNA fasta from iGenome folder (tRNA sequence are located in /Sequence/AbundantSequences folder)
            my $PATH_TO_FASTA = $tRNA_fasta;
            print $bowtieloc."-build ".$PATH_TO_FASTA." ".$bowtieIndexFile."\n";
            system($bowtieloc."-build ".$PATH_TO_FASTA." ".$bowtieIndexFile);
            systemError("Bowtie index creation error",$?,$!);
            print "done!\n\n";

        }
        if ($bowtieIndex =~ /sn-o-RNA/) {
            print "no Bowtie ". $bowtieIndexFile ." found\ncreating bowtie index ...\n";

            #Get snRNA fasta from iGenome folder (snRNA sequence are located in /Sequence/AbundantSequences folder)
            my $PATH_TO_FASTA = $snRNA_fasta;
            print $bowtieloc."-build ".$PATH_TO_FASTA." ".$bowtieIndexFile."\n";
            system($bowtieloc."-build ".$PATH_TO_FASTA." ".$bowtieIndexFile);
            systemError("Bowtie index creation error",$?,$!);
            print "done!\n\n";

        }
        elsif ($bowtieIndex =~ /phix/) {
            print "no Bowtie ". $bowtieIndexFile ."found\ncreating bowtie index ...\n";

            #Get phix fasta from iGenome folder (phix sequence(s) are located in /Sequence/AbundantSequences folder)
            my $PATH_TO_FASTA = $phix_fasta;
            print $bowtieloc."-build ".$PATH_TO_FASTA." ".$bowtieIndexFile."\n";
            system($bowtieloc."-build ".$PATH_TO_FASTA." ".$bowtieIndexFile);
            systemError("Bowtie index creation error",$?,$!);
            print "done!\n\n";

        }
        elsif ($bowtieIndex =~ /cdna/) {
            print "no Bowtie ". $bowtieIndexFile ."found\ncreating bowtie index ...\n";

            #Get phix fasta from iGenome folder (cdna sequence(s) are located in /Sequence/AbundantSequences folder)
            my $PATH_TO_FASTA = $cDNA_fasta;
            print $bowtieloc."-build ".$PATH_TO_FASTA." ".$bowtieIndexFile."\n";
            system($bowtieloc."-build ".$PATH_TO_FASTA." ".$bowtieIndexFile);
            systemError("Bowtie index creation error",$?,$!);
            print "done!\n\n";

        }
        elsif ($bowtieIndex =~ /genome/) {
            print "no Bowtie ". $bowtieIndexFile ."found\ncreating bowtie index ...\n";

            #Get genome fasta from iGenome folder (genome is located in .../WholeGenomeFasta/genome.fa)
            my $PATH_TO_FASTA = $genome_fasta;
            print $bowtieloc."-build ".$PATH_TO_FASTA." ".$bowtieIndexFile."\n";
            system($bowtieloc."-build ".$PATH_TO_FASTA." ".$bowtieIndexFile);
            systemError("Bowtie index creation error",$?,$!);
            print "done!\n\n";

        }
    }

}

sub map_topHat2 {

    # Catch
    my $seqFiles = $_[0];
    my $seqFileName = $_[1];
    my $clipper = $_[2];
    my $mismatch = $_[3];

    my @splitSeqFiles = split(/,/,$seqFiles);
    my $seqFile  = $splitSeqFiles[0];
    my $seqFile2 = $splitSeqFiles[1];
    $ref_loc = get_ref_loc($mapper);

    # Prepare
    my $directory = $work_dir."/".$mapper."/".$seqFileName."/";
    system "mkdir -p $directory";

    my ($fasta,$fasta1,$fasta2,$command,$full_command);

    # Get input fastq file
    if (uc($readtype) eq "RIBO" || uc($readtype) =~ m/SE/) {
        $fasta = $seqFile;
    } elsif (uc($readtype) =~ m/PE/) {
        $fasta1 = $seqFile;
        $fasta2 = $seqFile2;
    }

    # We need to first clip the adapter with fastx_clipper?
    if (uc($clipper) eq "FASTX") {

        print "     Clipping $seqFileName"." using fastx_clipper tool\n";
        # Without length cut-off and adaptor presence
        my $clip_command = $fastx_clip_loc." -Q33 -a ".$adaptorSeq." -l 26 -n –v -i ".$fasta." -o ".$work_dir."/fastq/$seqFileName"."_clip.fastq";
        print "     ".$clip_command."\n";
        system ($clip_command);
        $fasta = $work_dir."/fastq/$seqFileName"."_clip.fastq";
        print "clipfasta= $fasta\n";

        print "     Trimming $seqFileName"." using fastq_quality_trimmer tool\n";
        my $trim_command = $fastx_trim_loc." -Q33 -v -t 28 -l 26  -i ".$fasta." -o ".$work_dir."/fastq/$seqFileName"."_clip_trim.fastq";
        system ($trim_command);
        $fasta = $work_dir."/fastq/$seqFileName"."_clip_trim.fastq";
        print "trimfasta= $fasta\n";

    }

    #GO FOR Bowtie2-pHIX mapping
    if ($phix eq "Y") {
        #Check if pHIX STAR index exists?
        check_Bowtie_index($IndexPhix,$prev_mapper);
        print "     Mapping against Phix $seqFileName"."\n";

        # Build command
        system("mkdir  -p ".$work_dir."/fastq/nophix");

        #$command = $STAR_loc." --genomeLoad NoSharedMemory --outFilterMismatchNmax 2  --seedSearchStartLmaxOverLread .5 ".$clip_stat." --alignIntronMax 1 --genomeDir ".$ref_loc.$IndexPhix." --readFilesIn ".$fasta."  --outFileNamePrefix ".$work_dir."/fastq/nophix/ --runThreadN ".$cores." --outReadsUnmapped Fastx";
        $command = $prev_mapper." -p ".$cores." --".$bowtie2Setting." --sensitive-local --al ".$work_dir."/fastq/nophix/$seqFileName"."_phix.fq --un ".$work_dir."/fastq/nophix/$seqFileName"."_nophix.fq  -x ".$ref_loc."".$IndexPhix." -U ".$fasta.">/dev/null 2>".$work_dir."/fastq/nophix/bowtie.log";
        print "$command\n";
        # Run
        system($command);
        # Print
        print "   Finished phix multiseed mapping $seqFileName"."\n";
        # Process statistics
        open (STATS, ">>".$stat_file);
        my ($inReads, $mappedReadsU,$mappedReadsM, $unmappedReads) = parseLogBowtie($work_dir."/fastq/nophix/");
        my $mappedReads = $mappedReadsU + $mappedReadsM;
        print STATS "bowtie2 ".$run_name." ".$seqFileName." "."fastq phix ".$inReads."\n";
        print STATS "bowtie2 ".$run_name." ".$seqFileName." "."hit phix ".$mappedReads."\n";
        print STATS "bowtie2 ".$run_name." ".$seqFileName." "."unhit phix ".$unmappedReads."\n";
        close(STATS);

        $fasta = $work_dir."/fastq/nophix/".$seqFileName."_nophix.fq";

    }

    # If RIBO-SEQ: prior rRNA and/or snRNA and/or trRNA mapping is NECESSARY
    # if RNA-seq:  prior rRNA, snRNA and/or trRNA filtering is OPTIONAL
    #if (uc($readtype) eq "RIBO") {
        #GO FOR rRNA mapping
        if ($rRNA eq "Y") {
            # GO FOR BOWTIE2_LOCAL_SENSITIVE to get norRNA reads
            if ($prev_mapper eq 'bowtie2') {
                print "     Bowtie2 local version is selected for rRNA mapping\n"; #Bowtie local allows mismatches at beginning or end of read
                ###### Map against rRNA db
                # Check rRNA Index
                check_Bowtie_index($IndexrRNA,$prev_mapper);
                # Build command
                $ref_loc = get_ref_loc($prev_mapper);
                print "     Mapping against rRNA $seqFileName"."\n";
                $command = $bowtie2_loc." -p ".$cores." --".$bowtie2Setting." --sensitive-local --al ".$work_dir."/fastq/".$seqFileName."_rrna.fq --un ".$work_dir."/fastq/".$seqFileName."_norrna.fq  -x ".$ref_loc."".$IndexrRNA." -U ".$fasta." >/dev/null  2>".$work_dir."/fastq/bowtie.log";
                # Run
                print "$command\n";
                system($command);
                # Print
                print "   Finished rRNA multiseed alignment of $seqFileName\n";
                # Process statistics
                open (STATS, ">>".$stat_file);
                my ($inReads, $mappedReadsU,$mappedReadsM, $unmappedReads) = parseLogBowtie($work_dir."/fastq/");
                my $mappedReads = $mappedReadsU + $mappedReadsM;
                print STATS "bowtie2 ".$run_name." ".$seqFileName." "."fastq rRNA ".$inReads."\n";
                print STATS "bowtie2 ".$run_name." ".$seqFileName." "."hit rRNA ".$mappedReads."\n";
                print STATS "bowtie2 ".$run_name." ".$seqFileName." "."unhit rRNA ".$unmappedReads."\n";
                close(STATS);
            }
        $fasta = $work_dir."/fastq/$seqFileName"."_norrna.fq";
        }

        #GO FOR snRNA mapping
        if ($snRNA eq "Y") {
            # GO FOR BOWTIE2_LOCAL_SENSITIVE to get no snoRNA reads
            if ($prev_mapper eq 'bowtie2') {
                print "     Bowtie2 local version is selected for snoRNA mapping\n"; #Bowtie local allows mismatches at beginning or end of read
                ###### Map against snoRNA db
                # Check snoRNA Index
                check_Bowtie_index($IndexsnRNA,$prev_mapper);
                # Build command
                $ref_loc = get_ref_loc($prev_mapper);
                print "     Mapping against snoRNA $seqFileName"."\n";
                $command = $bowtie2_loc." -p ".$cores." --".$bowtie2Setting." --sensitive-local --al ".$work_dir."/fastq/".$seqFileName."_snrna.fq --un ".$work_dir."/fastq/".$seqFileName."_norrna_nosnrna.fq  -x ".$ref_loc."".$IndexsnRNA." -U ".$fasta.">/dev/null 2>".$work_dir."/fastq/bowtie.log";
                # Run
                system($command);
                # Print
                print "   Finished snoRNA multiseed alignment of $seqFileName\n";
                # Process statistics
                open (STATS, ">>".$stat_file);
                my ($inReads, $mappedReadsU,$mappedReadsM, $unmappedReads) = parseLogBowtie($work_dir."/fastq/");
                my $mappedReads = $mappedReadsU + $mappedReadsM;
                print STATS "bowtie2 ".$run_name." ".$seqFileName." "."fastq snRNA ".$inReads."\n";
                print STATS "bowtie2 ".$run_name." ".$seqFileName." "."hit snRNA ".$mappedReads."\n";
                print STATS "bowtie2 ".$run_name." ".$seqFileName." "."unhit snRNA ".$unmappedReads."\n";
                close(STATS);

            }
        $fasta = $work_dir."/fastq/$seqFileName"."_norrna_nosnrna.fq";
        }

        #GO FOR tRNA mapping
        if ($tRNA eq "Y") {
            # GO FOR BOWTIE2_LOCAL_SENSITIVE to get no tRNA reads
            if ($prev_mapper eq 'bowtie2') {
                print "     Bowtie2 local version is selected for tRNA mapping\n"; #Bowtie local allows mismatches at beginning or end of read
                ###### Map against tRNA db
                # Check tRNA Index
                check_Bowtie_index($IndextRNA,$prev_mapper);
                # Build command
                $ref_loc = get_ref_loc($prev_mapper);
                print "     Mapping against tRNA $seqFileName"."\n";
                $command = $bowtie2_loc." -p ".$cores." --".$bowtie2Setting." --sensitive-local --al ".$work_dir."/fastq/".$seqFileName."_trna.fq --un ".$work_dir."/fastq/".$seqFileName."_norrna_nosnrna_notrna.fq  -x ".$ref_loc."".$IndextRNA." -U ".$fasta.">/dev/null 2>".$work_dir."/fastq/bowtie.log";
                # Run
                system($command);
                # Print
                print "   Finished tRNA multiseed alignment of $seqFileName\n";
                # Process statistics
                open (STATS, ">>".$stat_file);
                my ($inReads, $mappedReadsU,$mappedReadsM, $unmappedReads) = parseLogBowtie($work_dir."/fastq/");
                my $mappedReads = $mappedReadsU + $mappedReadsM;
                print STATS "bowtie2 ".$run_name." ".$seqFileName." "."fastq tRNA ".$inReads."\n";
                print STATS "bowtie2 ".$run_name." ".$seqFileName." "."hit tRNA ".$mappedReads."\n";
                print STATS "bowtie2 ".$run_name." ".$seqFileName." "."unhit tRNA ".$unmappedReads."\n";
                close(STATS);
            }
        $fasta = $work_dir."/fastq/$seqFileName"."_norrna_nosnrna_notrna.fq";
        }
    #}

    # Run TopHat2
    $ref_loc = get_ref_loc($mapper);
    my $PATH_TO_GTF   = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Annotation/Genes/genes_".$ensemblversion".gtf";
    my $PATH_TO_GENOME_INDEX = $ref_loc."".$IndexGenome;

    check_Bowtie_index($IndexGenome,$prev_mapper);

    # alignment dependent on read type
    # Tophat parameters:
    #   --max-multihits: only ouput that many hits if read maps to multiple locations, if more equally scoring options are available than set with this value, a limited number is randomly selected)
    #   --report-secondary-alignments: this parameter shouldn't be set, suboptimal (i.e. with lower score) are not retained
    if (uc($readtype) eq "RIBO") {

        if ($unique eq 'Y' || $FirstRankMultiMap eq 'Y') {$maxmultimap = 1}  # to ensure unique mapping
        $command = $tophat2_loc." -N ". $mismatch." --max-multihits ".$maxmultimap." --no-coverage-search --segment-length 15 --output-dir ".$directory." -p ".$cores." --GTF ".$PATH_TO_GTF." ". $PATH_TO_GENOME_INDEX ." ". $fasta. " 2>&1";

    } elsif (uc($readtype) =~ m/PE/) {

        if ($unique eq 'Y' || $FirstRankMultiMap eq 'Y') {$maxmultimap = 1}  # to ensure unique mapping
        $command = $tophat2_loc." -N ". $mismatch." --max-multihits ".$maxmultimap." --no-coverage-search --segment-length 15 --output-dir ".$directory." -p ".$cores." --GTF ".$PATH_TO_GTF." ". $PATH_TO_GENOME_INDEX ." ". $fasta1. " ".$fasta2." 2>&1";
    } elsif (uc($readtype) =~ m/SE/) {

        if ($unique eq 'Y' || $FirstRankMultiMap eq 'Y') {$maxmultimap = 1}  # to ensure unique mapping
        $command = $tophat2_loc." -N ". $mismatch." --max-multihits ".$maxmultimap." --no-coverage-search --segment-length 15 --output-dir ".$directory." -p ".$cores." --GTF ".$PATH_TO_GTF." ". $PATH_TO_GENOME_INDEX ." ". $fasta. " 2>&1";
    }

    print "     ".$command."\n";
    system($command);
    systemError("STAR",$?,$!);

    # #  convert BAM back to SAM file
     print "converting BAM back to SAM...\n";
     system($samtools_loc." view -h -o ".$directory."Aligned.sorted.sam ".$directory."accepted_hits.bam > /dev/null 2>&1");
     systemError("Samtools view",$?,$!);

    # rename SAM output file
    print "renaming SAM output file...\n";
    # Bam file depends on what fastq file is processed (fastq1 = untreated, fastq2 = treated; that is for RIBO-seq experiments)
    my $bamf = ($seqFileName  eq 'fastq1') ? $out_sam_untr : $out_sam_tr;
    system("mv ".$directory."Aligned.sorted.sam ".$bamf);

    # remove redundant files
    system("rm ".$directory."Aligned.out.bam> /dev/null 2>&1");
    system("rm ".$directory."Aligned.out.sam> /dev/null 2>&1");
    system("rm ".$directory."Aligned.sorted.bam> /dev/null 2>&1");

    open (STATS, ">>".$stat_file);
    my ($inReads,$mappedReadsU,$mappedReadsM, $unmappedReads) = parseLogTopHat($seqFileName);
    print "$inReads,$mappedReadsU,$mappedReadsM, $unmappedReads\n";
    print STATS "TopHat2 ".$run_name." ".$seqFileName." "."fastq genomic ".$inReads."\n";
    print STATS "TopHat2 ".$run_name." ".$seqFileName." "."hitU genomic ".$mappedReadsU."\n";
    print STATS "TopHat2 ".$run_name." ".$seqFileName." "."hitM genomic ".$mappedReadsM."\n";
    print STATS "TopHat2 ".$run_name." ".$seqFileName." "."unhit genomic ".$unmappedReads."\n";
    close(STATS);

}

sub store_statistics {

    # DBH
    my ($file,$dsn,$us,$pw) = @_;
    my $dbh_sqlite_results = dbh($dsn,$us,$pw);


    ###############
    ## STATISTICS
    ##

    # Print
    print "Processing statistics!\n";

    # Gather statistics
    open (IN,"<".$stat_file) || die "ERROR";

    my $stat; my $ref; my $size; my $total;
    $total = "";

    while(my $line = <IN>){
        $line =~ s/^\s+//;
        my @a = split(/\s+/,$line);

        if ($line =~ m/^STAR/ ||  $line =~ m/^TopHat2/ || $line =~ m/^bowtie/) {
            my $run_name = $a[1];
            my $lane = $a[2];
            my $ref = $a[4];#($line =~ m/(hit|unhit) genomic/) ? "genomic" : "rRNA";
            my $type = $a[3];
            $stat->{$run_name.".".$lane}->{$ref}->{$type} = ($a[5]);
        }
        else {

            # Get some params
            $a[3] =~ /(fastq|norrna\.fq|rrna\.fq|hit|unhit)$/;
            my $type = ($1 eq "norrna.fq") ? "unhit" : ($1 eq "rrna.fq") ? "hit" : $1;

            # Get ref
            if ($a[3] =~ /\/fastq\//) { $ref = "rRNA";}
            elsif ($a[3] =~ /_(\d*)_(cDNA_hit|cDNA_unhit)/) { $ref = "cDNA"; }
            elsif ($a[3] =~ /_(\d*)_(genomic_hit|genomic_unhit)/) { $ref = "genomic";}

            # Get lane number
            $a[3] =~ /$fastqName(\d)/;
            my $lane = $1;

            # Fill hashref
            $stat->{$run_name."-".$lane}->{$ref}->{$type} = ($a[0] / 4);
        }

    }

    close(IN);
    print Dumper ($stat);

    my $query_table = "CREATE TABLE IF NOT EXISTS `statistics` (
    `sample` varchar(200) default NULL,
    `type` varchar(200) default NULL,
    `total` int(11) default NULL,
    `mapped_U` int(11) default NULL,
    `mapped_M` int(11) default NULL,
    `mapped_T` int(11) default NULL,
    `unmapped` int(11) default NULL,
    `map_freq_U` decimal(10,5) default NULL,
    `map_freq_M` decimal(10,5) default NULL,
    `map_freq_T` decimal(10,5) default NULL
    )";

    $dbh_sqlite_results->do($query_table);

    foreach my $sample (keys %$stat){

        my $prev_ref;
        foreach my $ref (keys %{$stat->{$sample}}) {
            #The previous ma
            if (uc($mapper) eq 'STAR' || uc($mapper) eq 'TOPHAT2' || uc($mapper) eq 'BOWTIE' || uc($mapper) eq 'BOWTIE2') {
                if (uc($readtype) eq 'RIBO') {
                    $prev_ref = ($ref eq "genomic") ? "rRNA" : "";
                    $total = ($ref eq "rRNA" || $ref eq "phix" || $ref eq "snRNA" || $ref eq "tRNA") ? $stat->{$sample}->{$ref}->{"fastq"} : $stat->{$sample}->{$prev_ref}->{"unhit"};

                } elsif (uc($readtype) =~ m/POLYA/) {
                    $prev_ref = "";
                    $total = ($ref eq "genomic") ? $stat->{$sample}->{$ref}->{"fastq"} : "";
                }
            }

            my $query;


            if ((uc($mapper) eq "STAR" || uc($mapper) eq "TOPHAT2" || uc($mapper) eq 'BOWTIE' || uc($mapper) eq 'BOWTIE2') && $ref eq "genomic" ) {

                my $freq_U =  $stat->{$sample}->{$ref}->{"hitU"} / $stat->{$sample}->{$ref}->{"fastq"};
                my $freq_M =  $stat->{$sample}->{$ref}->{"hitM"} / $stat->{$sample}->{$ref}->{"fastq"};
                my $freq_T = $freq_U + $freq_M;
                my $map_T  = $stat->{$sample}->{$ref}->{"hitM"} + $stat->{$sample}->{$ref}->{"hitU"};
                my $total  = $stat->{$sample}->{$ref}->{"fastq"};

                $query = "INSERT INTO statistics (sample,type,total,mapped_U,mapped_M,mapped_T,unmapped,map_freq_U,map_freq_M,map_freq_T) VALUES (\'".$sample."\',\'".$ref."\',\'".$total."\',\'".$stat->{$sample}->{$ref}->{"hitU"}."\',\'".$stat->{$sample}->{$ref}->{"hitM"}."\',\'".$map_T."\',\'".$stat->{$sample}->{$ref}->{"unhit"}."\',\'".$freq_U."\',\'".$freq_M."\',\'".$freq_T."\')";

            } else {

                my $freq =  ($prev_ref eq "")       ? $stat->{$sample}->{$ref}->{"hit"} / $stat->{$sample}->{$ref}->{"fastq"}
                :  ($prev_ref eq "rRNA")   ? $stat->{$sample}->{$ref}->{"hit"} / $stat->{$sample}->{$prev_ref}->{"unhit"}
                :  ($prev_ref eq "cDNA")   ? $stat->{$sample}->{$ref}->{"hit"} / $stat->{$sample}->{$prev_ref}->{"unhit"} : "";

                $query = "INSERT INTO statistics (sample,type,total,mapped_T,unmapped,map_freq_T) VALUES (\'".$sample."\',\'".$ref."\',\'".$total."\',\'".$stat->{$sample}->{$ref}->{"hit"}."\',\'".$stat->{$sample}->{$ref}->{"unhit"}."\',\'".$freq."\')";
            }

            $dbh_sqlite_results->do($query);

        }
    }

    system("rm ".$stat_file);

}

###  Store all input variables in an SQLite table
sub store_input_vars {

    # Catch
    my $dsn = $_[0];
    my $us = $_[1];
    my $pw =  $_[2];
    my $run_name = $_[3];
    my $ensembl_version = $_[4];
    my $species = $_[5];
    my $mapper = $_[6];
    my $unique = $_[7];
    my $adaptor= $_[8];
    my $readlength= $_[9];
    my $readtype= $_[10];
    my $IGENOMES_ROOT = $_[11];
    my $nr_of_cores = $_[12];

    my $dbh_sqlite_results = dbh($dsn,$us,$pw);

    my $query_table = "CREATE TABLE IF NOT EXISTS `arguments` (
    `variable` varchar(200) default NULL,
    `value` varchar(200) default NULL
    )";


    $dbh_sqlite_results->do($query_table);


    my $query = "INSERT INTO arguments (variable,value) VALUES (\'run_name\',\'".$run_name."\')";
    $dbh_sqlite_results->do($query);

    $query = "INSERT INTO arguments (variable,value) VALUES (\'ensembl_version\',\'".$ensembl_version."\')";
    $dbh_sqlite_results->do($query);

    $query = "INSERT INTO arguments (variable,value) VALUES (\'species\',\'".$species."\')";
    $dbh_sqlite_results->do($query);

    $query = "INSERT INTO arguments (variable,value) VALUES (\'mapper\',\'".$mapper."\')";
    $dbh_sqlite_results->do($query);

    $query = "INSERT INTO arguments (variable,value) VALUES (\'unique\',\'".$unique."\')";
    $dbh_sqlite_results->do($query);

    $query = "INSERT INTO arguments (variable,value) VALUES (\'adaptor\',\'".$adaptor."\')";
    $dbh_sqlite_results->do($query);

    $query = "INSERT INTO arguments (variable,value) VALUES (\'readlength\',\'".$readlength."\')";
    $dbh_sqlite_results->do($query);

    $query = "INSERT INTO arguments (variable,value) VALUES (\'readtype\',\'".$readtype."\')";
    $dbh_sqlite_results->do($query);

    $query = "INSERT INTO arguments (variable,value) VALUES (\'igenomes_root\',\'".$IGENOMES_ROOT."\')";
    $dbh_sqlite_results->do($query);

    $query = "INSERT INTO arguments (variable,value) VALUES (\'nr_of_cores\',\'".$nr_of_cores."\')";
    $dbh_sqlite_results->do($query);
}

### GET INDEX LOCATION ###
sub get_ref_loc {

    # Catch
    my $mapper  = $_[0];

    my $ref_loc = (uc($mapper) eq "BOWTIE") ? $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Sequence/BowtieIndex/"                 #Bowtie indexes
    : (uc($mapper) eq "BOWTIE2" || uc($mapper) eq "TOPHAT2") ? $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Sequence/Bowtie2Index/"   #Bowtie2 indexes
    : $STAR_ref_loc; #STAR indexes

    return($ref_loc);

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

### SYSTEMERROR ###
sub systemError {
    my ($command,$returnValue,$errorMessage) = @_;
    if ($returnValue == -1){
        die "$command failed!\n$errorMessage\n\n";
    }
}

sub parseLogSTAR  {

    # Catch
    my $dir  = $_[0];

    open (LOG,"<".$dir."Log.final.out") || die "path $dir";

    my ($line,@lineSplit,$inReads,$mappedReadsU,$mappedReadsM,$mappedReads,$unmappedReads);

    while($line = <LOG>) {
        @lineSplit = split(/\|/,$line);
        if ($lineSplit[0] =~ /Number of input reads/) {
            $lineSplit[1] =~s/^\s+//;
            $lineSplit[1] =~s/\s+$//;
            $inReads = $lineSplit[1];
        }
        if ($lineSplit[0] =~ /Uniquely mapped reads number/) {
            $lineSplit[1] =~s/^\s+//;
            $lineSplit[1] =~s/\s+$//;
            $mappedReadsU = $lineSplit[1];
        }
        if ($lineSplit[0] =~ /Number of reads mapped to multiple loci/) {
            $lineSplit[1] =~s/^\s+//;
            $lineSplit[1] =~s/\s+$//;
            $mappedReadsM = $lineSplit[1];
        }
    }
    $mappedReads = $mappedReadsU + $mappedReadsM;
    $unmappedReads = $inReads - $mappedReads;
    close(LOG);
    return ($inReads,$mappedReadsU,$mappedReadsM,$unmappedReads);
}


sub parseLogTopHat {

    # Catch
    my $seqFileName  = $_[0];
    my $directory = $work_dir."/".$mapper."/".$seqFileName."/";

    open (LOG1,"<".$directory."align_summary.txt") || die "ERROR";

    my ($line,$inReads,$mappedReadsU,$mappedReadsM,$mappedReads,$unmappedReads,$mappedReadsMcorr);
    # Quick Fix for TopHat2 parsing (check the parsing of the file)
    $mappedReads = 1; $inReads = 1;

    while($line = <LOG1>) {

        if ($line =~ /Input\s+:\s+(\d+)/) {
            $inReads= $1;
        }
        if ($line =~ /Mapped\s+:\s+(\d+)/) {
            $mappedReads= $1;
        }

        if ($line =~ /of these:\s+(\d+)\s+\(.*\((\d+)/) {
            $mappedReadsM= $1;
            $mappedReadsMcorr=$2;
        }
    }

    $mappedReadsU  = $mappedReads - $mappedReadsM;
    $mappedReadsM  = $mappedReadsM - $mappedReadsMcorr;
    $mappedReads   = $mappedReads - $mappedReadsMcorr;
    $unmappedReads = $inReads - $mappedReads;


    close(LOG1);

    return ($inReads,$mappedReadsU,$mappedReadsM,$unmappedReads);

}

sub parseLogBowtie {

  #  205829 reads; of these:
  #  205829 (100.00%) were unpaired; of these:
  #    205829 (100.00%) aligned 0 times
  #    0 (0.00%) aligned exactly 1 time
  #    0 (0.00%) aligned >1 times
  # 0.00% overall alignment rate


   # Catch
    my $dir  = $_[0];

    open (LOG1,"<".$dir."bowtie.log") || die "path $dir";
    my ($line,$inReads,$mappedReadsU,$mappedReadsM,$mappedReads,$unmappedReads,$mappedReadsMcorr);

    while($line = <LOG1>) {

        if ($line =~ /(\d+)\s+reads; of these:/) {
            $inReads= $1;
            print "inReads = $inReads\n";
        }
        if ($line =~ /\s+(\d+).*aligned 0 times/) {
            $unmappedReads= $1;
            print "unmappedReads = $unmappedReads\n";
        }
        if ($line =~ /\s+(\d+).*aligned exactly 1 time/) {
            $mappedReadsU= $1;
            print "mappedReadsU = $mappedReadsU\n";
        }
        if ($line =~ /\s+(\d+).*aligned >1 times/) {
            $mappedReadsM= $1;
            print "mappedReadsM = $mappedReadsM\n";
        }

    }

    $mappedReads   = $mappedReadsU + $mappedReadsM;
    close(LOG1);
    return ($inReads,$mappedReadsU,$mappedReadsM,$unmappedReads);

}


sub RIBO_parse_store {

    # Catch
    my $seqFile = $_[0];
    my $seqFileName = $_[1];
	my $uniq = $_[2];

    my $bedgr_s = ($seqFileName  eq 'fastq1') ? $out_bg_s_untr : $out_bg_s_tr;
    my $bedgr_as = ($seqFileName  eq 'fastq1') ? $out_bg_as_untr : $out_bg_as_tr;
    my $sam = ($seqFileName  eq 'fastq1') ? $out_sam_untr : $out_sam_tr;

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
    split_SAM_per_chr(\%chr_sizes,$work_dir,$seqFileName,$run_name,$sam,$uniq);

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

    # Init multi core
    my $pm = new Parallel::ForkManager($cores);
    print "   Using ".$cores." core(s)\n   ---------------\n";

    foreach my $chr (keys %chr_sizes){

        ### Start parallel process
        $pm->start and next;

        ### DBH per process
        my $dbh = dbh($dsn_sqlite_results,$us_sqlite_results,$pw_sqlite_results);

        ### RIBO parsing
        my ($hits,$hits_splitRPF) = RIBO_parsing_genomic_per_chr($work_dir,$seqFileName,$run_name,$sam,$chr);

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
    system("rm -rf ".$TMP."/genomic/");
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
    split_SAM_per_chr(\%chr_sizes, $work_dir, $seqFileName, $run_name, $sam, $uniq);

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
            next;
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
            next;
        } else {
            print BEDALLGRAS "chr".$_;
        }
    }
    close(OLDBEDGRAS);
    system("rm -rf ".$TMP."/genomic/bedgraph_old_antisense.bedgraph");
    close(BEDALLGRAS);
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

    # Get chrs from Chr_File
    open (Q,"<".$chr_file) || die "Cannot open chr sizes input\n";
    while ($line = <Q>){
        $line =~ /^(\S*)/;
        push (@chr,$1);
    }

    # Get correct coord_system_id
    my $query = "SELECT coord_system_id FROM coord_system where name = 'chromosome' and version = '".$assembly."'";
	my $sth = $dbh->prepare($query);
	$sth->execute();
    @coord_system = $sth->fetchrow_array;
    $coord_system_id = $coord_system[0];
    $sth->finish();

    # Get chrs with seq_region_id
    my $chr;
    foreach (@chr){
        if (uc($species) eq "FRUITFLY"){
            if($_ eq "M"){
                $chr = "dmel_mitochondrion_genome";
            }else{
                $chr = $_;
            }
        }else {
            $chr = $_;
        }

        my $query = "SELECT seq_region_id FROM seq_region where coord_system_id = ".$coord_system_id."  and name = '".$chr."' ";
        my $sth = $dbh->prepare($query);
        $sth->execute();
        @ids = $sth->fetchrow_array;
        $seq_region_id = $ids[0];
        $chrs->{$_}{'seq_region_id'} = $seq_region_id;
        $sth->finish();
    }

    #Disconnect DBH
    $dbh->disconnect();

	# Return
	return($chrs);

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
                    if ( $genmatchL >= 25 && $genmatchL <= 34) {
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
            ($offset,$genmatchL,$intron_total,$extra_for_min_strand) = parse_RIBO_CIGAR($CIGAR,$strand);
            $lendistribution->{$genmatchL}++;
            #Determine genomic position based on CIGAR string output and mapping position and direction
            $start = ($strand eq "+") ? $mapping_store[3] + $offset + $intron_total : ($strand eq "-") ? $mapping_store[3] - $offset - $intron_total + $extra_for_min_strand -1 : "";

            $hits_genomic_splitRPF->{$chr}->{$start}->{$genmatchL}->{$strand}++;

            if ( $genmatchL >= 26 && $genmatchL <= 34) {
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
    my $offset = get_offset($genmatchL);
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

    my $offset = ($len >= 34) ? 14 :
    ($len <= 30) ? 12 :
    13;

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


### PROGRESS BAR ###
sub progressbar {
    $| = 1;
    my $a=$_[0];
    if($a<0) {
        for(my $s=0;$s<-$a/1000000;$s++) {
            print " ";
        }
        print "                                \r";
    }
    my $seq="";
    if($a%50==0) {
        for(my $s=0;$s<$a/1000000;$s++) {
            $seq.="►";
        }
    }
    if(($a/1000000)%2==0) { print "  ✖ Crunching ".$seq."\r"; }
    if(($a/1000000)%2==1) { print "  ✚ Crunching ".$seq."\r"; }
}
