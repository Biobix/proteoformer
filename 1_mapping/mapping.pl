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
# ./mapping.pl --name mESC --species mouse --ensembl 72 --cores 20 --readtype ribo --unique N --inputfile1 file1 --inputfile2 file2 --igenomes_root IGENOMES_ROOT (--mapper STAR --adaptor CTGTAGGCACCATCAAT --readlength 36 --truseq Y --firstrankmultimap N --out_bg_s_untr bg_s_untr --out_bg_as_untr bg_as_untr --out_bg_s_tr bg_s_tr --out_bg_as_tr bg_as_tr --out_sam_untr sam_untr --out_sam_tr sam_tr --out_bam_untr bam_untr --out_bam_tr bam_tr --out_sqlite sqliteDBName --work_dir getcwd --tmpfolder $TMP --suite custom)

#For GALAXY
#mapping.pl --name "${experimentname}" --species "${organism}" --ensembl "${ensembl}" --cores "${cores}" --readtype $readtype.riboSinPair --unique "${unique}" --mapper "${mapper}" --readlength $readtype.readlength --adaptor $readtype.adaptor --inputfile1 $readtype.input_file1 --inputfile2 $readtype.input_file2 --out_bg_s_untr "${untreat_s_bg}"  --out_bg_as_untr "${untreat_as_bg}" --out_bg_s_tr "${treat_s_bg}" --out_bg_as_tr "${treat_as_bg}" --out_sam_untr "${untreat_sam}" --out_sam_tr "${treat_sam}" --out_sqlite "${out_sqlite}" --igenomes_root "${igenomes_root}"

# get the command line arguments
my ($work_dir,$run_name,$species,$ensemblversion,$cores,$mapper,$readlength,$readtype,$truseq,$tmpfolder,$adaptorSeq,$unique,$seqFileName1,$seqFileName2,$fastqName,$min_l_plastid,$max_l_plastid,$offset_img_untr,$offset_img_tr,$min_l_parsing,$max_l_parsing,$out_bg_s_untr,$out_bg_as_untr,$out_bg_s_tr,$out_bg_as_tr,$out_sam_untr,$out_sam_tr,$out_bam_untr,$out_bam_tr,$out_sqlite,$IGENOMES_ROOT,$ref_loc,$clipper,$phix,$rRNA,$snRNA,$tRNA,$tr_coord,$maxmultimap,$mismatch,$out_bam_tr_untr,$out_bam_tr_tr,$splicing,$FirstRankMultiMap,$rpf_split,$suite,$suite_tools_loc);

GetOptions(
"inputfile1=s"=>\$seqFileName1,         	# the fastq file of the untreated data for RIBO-seq (no,CHX,EMT) or the 1st fastq for single/paired-end RNA-seq                  mandatory argument
"inputfile2=s"=>\$seqFileName2,         	# the fastq file of the treated data for RIBO-seq (PUR,LTM,HARR) or the 2nd fastq for paired-end RNA-seq                         mandatory argument
"name=s"=>\$run_name,                   	# Name of the run,                                                  			mandatory argument
"species=s"=>\$species,                 	# Species, eg mouse/rat/human/fruitfly/arabidopsis/SL1344/zebrafish/yeast       mandatory argument
"ensembl=i"=>\$ensemblversion,          	# Ensembl annotation version, eg 66 (Feb2012),                      			mandatory argument
"cores=i"=>\$cores,                     	# Number of cores to use for Bowtie Mapping,                        			mandatory argument
"readtype=s"=>\$readtype,              		# The readtype (ribo, ribo_untr, PE_polyA, SE_polyA, PE_total, SE_total)       	mandatory argument (default = ribo)
"mapper:s"=>\$mapper,                   	# The mapper used for alignment (STAR,TopHat2)       			                optional  argument (default = STAR)
"readlength:i"=>\$readlength,           	# The readlength (if RiboSeq take 50 bases),                        			optional  argument (default = 50)
"adaptor:s"=>\$adaptorSeq,              	# The adaptor sequence that needs to be clipped with fastx_clipper, 			optional  argument (default = CTGTAGGCACCATCAAT) => Ingolia paper (for ArtSeq = AGATCGGAAGAGCACAC)
"unique=s" =>\$unique,                  	# Retain the uniquely (and multiple) mapping reads (Y or N),        			mandatory argument
"tmp:s" =>\$tmpfolder,                  	# Folder where temporary files are stored,                          			optional  argument (default = $TMP or $CWD/tmp env setting)
"work_dir:s" =>\$work_dir,              	# Working directory ,                                               			optional  argument (default = $CWD env setting)
"min_l_plastid:s" =>\$min_l_plastid,        # Minimum length for plastid                                                    optional  argument (default = 22)
"max_l_plastid:s" =>\$max_l_plastid,        # Maximum length for plastid                                                    optional  argument (default = 34)
"offset_img_untr:s" =>\$offset_img_untr,    # Path to save the offset image of plastid in (untreated)                       optional  argument (default = CWD/plastid/run_name_untreated_p_offsets.png)
"offset_img_tr:s" =>\$offset_img_tr,        # Path to save the offset image of plastid in (treated)                         optional  argument (default = CWD/plastid/run_name_treated_p_offsets.png)
"min_l_parsing:s" =>\$min_l_parsing,        # Minimum length for count table parsing                                        optional  argument (default = 26, 25 for fruitfly)
"max_l_parsing:s" =>\$max_l_parsing,        # Maximum length for count table parsing                                        optional  argument (default = 34)
"out_bg_s_untr:s" =>\$out_bg_s_untr,    	# Output file for sense untreated count data (bedgraph)             			optional  argument (default = untreat_sense.bedgraph)
"out_bg_as_untr:s" =>\$out_bg_as_untr,  	# Output file for antisense untreated count data (bedgraph)         			optional  argument (default = untreat_antisense.bedgraph)
"out_bg_s_tr:s" =>\$out_bg_s_tr,        	# Output file for sense treated count data (bedgraph)               			optional  argument (default = treat_sense.bedgraph)
"out_bg_as_tr:s" =>\$out_bg_as_tr,      	# Output file for antisense treated count data (bedgraph)           			optional  argument (default = treat_antisense.bedgraph)
"out_sam_untr:s" =>\$out_sam_untr,      	# Output file for alignments of untreated data (sam)                			optional  argument (default = untreat.sam)
"out_sam_tr:s" =>\$out_sam_tr,         	  # Output file for alignments of treated data (sam)                  			optional  argument (default = treat.sam)
"out_bam_untr:s" =>\$out_bam_untr,      	# Output file for alignments of untreated data (bam)                			optional  argument (default = untreat.bam)
"out_bam_tr:s" =>\$out_bam_tr,         	  # Output file for alignments of treated data (bam)                  			optional  argument (default = treat.bam)
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
"firstrankmultimap=s" =>\$FirstRankMultiMap,  # Only retain the first ranked alignment of multimapper (Y or N)           optional argument (default = N)
"rpf_split=s" =>\$rpf_split,                 #If the program needs to construct RPF specific bedgraph files (Y or N)     optional argument (default = N)
"suite=s" =>\$suite,                       #Option to execute different mapping modules all together for ribo data ('custom', 'standard', 'plastid')     optional argument (default = custom)
                                                        #Custom: only mapping, other modules manually afterwards
                                                        #Standard: mapping + parsing with standard offset
                                                        #Plastid: mapping + default p offset calculation with plastid + parsing based on these offsets
"suite_tools_loc=s" =>\$suite_tools_loc     #The folder with script of subsequent tools when using a suite                  optional argument (default = workdir)
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
} elsif (!defined($seqFileName2) && (uc($readtype) eq 'RIBO' || uc($readtype) eq 'PE_polyA' || uc($readtype) eq 'PE_total')) {
    die "\nDon't forget to pass the FastQ file for treated RIBO-seq (PUR,LTM,HARR) or 2nd fastq for paired-end RNA-seq using the --file or -f argument!\n\n";
} else {
    if ($readtype eq 'ribo_untr'){
        $seqFileName2 = "Ribo-seq only untreated";
    } else {
        $seqFileName2 = "Single End RNA sequencing";
    }
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
if ($min_l_plastid){
    print "The minimum length for plastid                           : $min_l_plastid\n";
} else {
    #Choose default value
    $min_l_plastid = 22;
    print "The minimum length for plastid                           : $min_l_plastid\n";
}
if ($max_l_plastid){
    print "The maximum length for plastid                           : $max_l_plastid\n";
} else {
    #Choose default value
    $max_l_plastid = 34;
    print "The maximum length for plastid                           : $max_l_plastid\n";
}
if ($offset_img_untr){
    print "The path for the untreated offset image                   : $offset_img_untr\n";
} else {
    $offset_img_untr = $work_dir."/plastid/".$run_name."_untreated_p_offsets.png";
    print "The path for the untreated offset image                    : $offset_img_untr\n";
}
if ($offset_img_tr){
    print "The path for the treated offset image                     : $offset_img_tr\n";
} else {
    $offset_img_tr = $work_dir."/plastid/".$run_name."_treated_p_offsets.png";
    print "The path for the treated offset image                      : $offset_img_tr\n";
}
if ($min_l_parsing){
    print "The minimum length for count table parsing               : $min_l_parsing\n";
} else {
    if(uc($species) eq "FRUITFLY"){
        $min_l_parsing = 25;
    } else {
        $min_l_parsing = 26;
    }
    print "The minimum length for count table parsing               : $min_l_parsing\n";
}
if ($max_l_parsing){
    print "The maximum length for count table parsing               : $max_l_parsing\n";
} else {
    $max_l_parsing = 34;
    print "The maximum length for count table parsing               : $max_l_parsing\n";
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
if ($rpf_split){
    print "RPF specific bedgraphs                                   : $rpf_split\n";
} else {
    $rpf_split = 'N';
    print "RPF specific bedgraphs                                   : $rpf_split\n";
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
if($readtype eq "ribo" || $readtype eq "ribo_untr"){
    if($suite eq "custom" || $suite eq "standard" || $suite eq "plastid"){
        print "Mapping suite                                               : $suite\n";
        if ($suite_tools_loc){
            print "Mapping suite tools folder                                       : $suite_tools_loc\n";
        } else {
            $suite_tools_loc = $work_dir;
            print "Mapping suite tools folder                                       : $suite_tools_loc\n";
        }
        if($suite eq "standard"){
            if(!-e $suite_tools_loc."/mapping_parsing.pl"){
                die "Parsing script not found!!!";
            }
        } elsif ($suite eq "plastid"){
            if(!-e $suite_tools_loc."/mapping_parsing.pl"){
                die "Parsing script not found!!!";
            } elsif(!-e $suite_tools_loc."/mapping_plastid.pl"){
                die "Plastid script not found!!!";
            }
        }
    } else {
        $suite = "custom";
        print "Mapping suite                                               : $suite\n";
    }
} else {
    if($suite eq "standard" || $suite eq "plastid"){
        die "Standard and plastid suite only for RIBO data!";
    } elsif (!defined($suite)){
        $suite = "custom";
    }
    print "Mapping suite                                               : $suite\n";
}

# Create output directory
system "mkdir -p ".$work_dir."/output/";
system "mkdir -p ".$work_dir."/fastq/";
if ($rpf_split eq 'Y'){
    system "mkdir -p ".$work_dir."/output/rpf_specific_bedgraph/";
}

if (!defined($out_bg_s_untr))  		{$out_bg_s_untr     = $work_dir."/output/untreat_sense.bedgraph";}
if (!defined($out_bg_as_untr)) 		{$out_bg_as_untr    = $work_dir."/output/untreat_antisense.bedgraph";}
if (!defined($out_bg_s_tr))    		{$out_bg_s_tr       = $work_dir."/output/treat_sense.bedgraph";}
if (!defined($out_bg_as_tr))   		{$out_bg_as_tr      = $work_dir."/output/treat_antisense.bedgraph";}
if (!defined($out_sam_untr))   		{$out_sam_untr      = $work_dir."/".$mapper."/fastq1/untreat.sam";}
if (!defined($out_sam_tr))     		{$out_sam_tr        = $work_dir."/".$mapper."/fastq2/treat.sam";}
if (!defined($out_bam_untr))   		{$out_bam_untr      = $work_dir."/".$mapper."/fastq1/untreat.bam";}
if (!defined($out_bam_tr))     		{$out_bam_tr        = $work_dir."/".$mapper."/fastq2/treat.bam";}
if (!defined($out_bam_tr_untr)) 	{$out_bam_tr_untr    = $work_dir."/".$mapper."/fastq1/untreat_tr.bam";}
if (!defined($out_bam_tr_tr))		{$out_bam_tr_tr     = $work_dir."/".$mapper."/fastq2/treat_tr.bam";}
if (!defined($out_sqlite))     		{$out_sqlite        = $work_dir."/SQLite/results.db";}
print "SQLite database                                          : $out_sqlite\n";

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
my $spec = (uc($species) eq "MOUSE") ? "Mus_musculus" : (uc($species) eq "RAT") ? "Rattus_norvegicus" : (uc($species) eq "SL1344") ? "SL1344" : (uc($species) eq "HUMAN") ? "Homo_sapiens" : (uc($species) eq "ARABIDOPSIS") ? "Arabidopsis_thaliana" : (uc($species) eq "FRUITFLY") ? "Drosophila_melanogaster" : (uc($species) eq "YEAST") ? "Saccharomyces_cerevisiae" : (uc($species) eq "ZEBRAFISH") ? "Danio_rerio" : "";
my $spec_short = (uc($species) eq "MOUSE") ? "mmu" : (uc($species) eq "RAT") ? "rnor" : (uc($species) eq "SL1344") ? "sl1344" : (uc($species) eq "HUMAN") ? "hsa" : (uc($species) eq "ARABIDOPSIS") ? "ath" : (uc($species) eq "FRUITFLY") ? "dme" : (uc($species) eq "YEAST") ? "sce" : (uc($species) eq "ZEBRAFISH") ? "dre" : "";
#Old mouse assembly = NCBIM37, new one is GRCm38. Old human assembly = GRCh37, the new one is GRCh38
my $assembly = (uc($species) eq "MOUSE" && $ensemblversion >= 70 ) ? "GRCm38"
: (uc($species) eq "MOUSE" && $ensemblversion < 70 ) ? "NCBIM37"
: (uc($species) eq "RAT" && $ensemblversion >=80 ) ? "Rnor_6.0"
: (uc($species) eq "RAT" && $ensemblversion < 80) ? "Rnor_5.0"
: (uc($species) eq "HUMAN" && $ensemblversion >= 76) ? "GRCh38"
: (uc($species) eq "HUMAN" && $ensemblversion < 76) ? "GRCh37"
: (uc($species) eq "ARABIDOPSIS") ? "TAIR10"
: (uc($species) eq "SL1344") ? "ASM21085v2"
: (uc($species) eq "ZEBRAFISH") ? "GRCz10"
: (uc($species) eq "YEAST") ? "R64-1-1"
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
store_input_vars($dsn_sqlite_results,$us_sqlite_results,$pw_sqlite_results,$run_name_short,$ensemblversion,$species,$mapper,$unique,$adaptorSeq,$readlength,$readtype,$IGENOMES_ROOT,$cores, $seqFileName1, $seqFileName2,$rpf_split,$FirstRankMultiMap,$truseq,$out_bg_s_untr,$out_bg_as_untr,$out_bg_s_tr,$out_bg_as_tr,$out_sam_untr,$out_sam_tr,$maxmultimap,$out_bam_untr,$out_bam_tr,$min_l_plastid,$max_l_plastid,$min_l_parsing,$max_l_parsing);


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
        }
        if ($readtype eq "ribo_untr"){
            map_topHat2($_,$fastqName,$clipper,$mismatch);
        }
        if (uc($readtype) eq "SE_POLYA") {
            map_topHat2($_,$fastqName,$clipper,$mismatch);
        }
        if (uc($readtype) =~ m/PE/) {
            map_topHat2($_,$fastqName,$clipper,$mismatch);
        }
        if (uc($readtype) eq 'SE_TOTAL') {
            map_STAR($_,$fastqName,$clipper,$mismatch);
        }
        my $end = time - $start;
        printf("runtime TopHat against genomic: %02d:%02d:%02d\n\n",int($end/3600), int(($end % 3600)/60), int($end % 60));
    }
    elsif (uc($mapper) eq "STAR") {

		print "Mapping with STAR\n";
        my $start = time;
        if (uc($readtype) eq "RIBO") {
            map_STAR($_,$fastqName,$clipper,$mismatch);
        }
        if ($readtype eq "ribo_untr"){
            map_STAR($_,$fastqName,$clipper,$mismatch);
        }
        if (uc($readtype) eq "SE_POLYA") {
            map_STAR($_,$fastqName,$clipper,$mismatch);
        }
        if (uc($readtype) =~ m/PE/) {
            map_STAR($_,$fastqName,$clipper,$mismatch);
        }
        if (uc($readtype) eq 'SE_TOTAL') {
            map_STAR($_,$fastqName,$clipper,$mismatch);
        }

        my $end = time - $start;
        printf("runtime STAR against genomic: %02d:%02d:%02d\n\n",int($end/3600), int(($end % 3600)/60), int($end % 60));
    }
    # Store statistics in DB
    store_statistics($stat_file,$dsn_sqlite_results,$us_sqlite_results,$pw_sqlite_results);

}



### Suite options ###
if($suite eq "standard"){
    print "\n\n\n\n\n\n\n\t\tS U I T E:     GO THROUGH WITH MAPPING PARSING\n\n\n";
    system("perl ".$suite_tools_loc."/mapping_parsing.pl --out_sqlite ".$out_sqlite." --offset ".$suite);
} elsif ($suite eq "plastid"){
    print "\n\n\n\n\n\n\n\t\tS U I T E:     GO THROUGH WITH PLASTID (UNTREATED)\n\n\n";
    system("perl ".$suite_tools_loc."/mapping_plastid.pl --out_sqlite ".$out_sqlite." --treated untreated  --offset_img ".$offset_img_untr);
    
    if ($readtype eq "ribo"){#Treated sample only for ribo experiment with two sample treatment types
        print "\n\n\n\n\n\n\n\t\tS U I T E:     GO THROUGH WITH PLASTID (TREATED)\n\n\n";
        system("perl ".$suite_tools_loc."/mapping_plastid.pl --out_sqlite ".$out_sqlite." --treated treated --offset_img ".$offset_img_tr);
    }
    
    print "\n\n\n\n\n\n\n\t\tS U I T E:     GO THROUGH WITH MAPPING PARSING\n\n\n";
    system("perl ".$suite_tools_loc."/mapping_parsing.pl --out_sqlite ".$out_sqlite." --offset ".$suite);
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
    system ($fastx_trim_loc." -Q33 -v -t 28 -l 20 -i ".$work_dir."/fastq/$seqFileName"."_clip.fastq -o ".$work_dir."/fastq/$seqFileName"."_trim.fastq");

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
        system ($fastx_trim_loc." -Q33 -v -t 28 -l 20 -i ".$work_dir."/fastq/$seqFileName"."_clip.fastq -o ".$work_dir."/fastq/$seqFileName"."_trim.fastq");
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
    if (uc($readtype) eq "RIBO" || $readtype eq "ribo_untr" || uc($readtype) =~ m/SE/) {
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
        my $clip_command = $fastx_clip_loc." -Q33 -a ".$adaptorSeq." -l 20 -n –v -i ".$fasta." -o ".$work_dir."/fastq/".$seqFileName."_clip.fastq";
        print "     ".$clip_command."\n";
        system ($clip_command);
        $fasta = $work_dir."/fastq/$seqFileName"."_clip.fastq";
        print "clipfasta= $fasta\n";

	    print "     Trimming $seqFileName"." using fastq_quality_trimmer tool\n";
		my $trim_command = $fastx_trim_loc." -Q33 -v -t 28 -l 20  -i ".$fasta." -o ".$work_dir."/fastq/".$seqFileName."_clip_trim.fastq";
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
		print "Phix mapping ...\n";
		print "$command\n";
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
				system ($fastx_trim_loc." -Q33 -v -t 28 -l 20 -i ".$work_dir."/fastq/$seqFileName"."_clip.fastq -o ".$work_dir."/fastq/$seqFileName"."_trim.fastq");

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

    if (uc($readtype) eq "RIBO" || $readtype eq "ribo_untr") {

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

    # # Output bamfile
    my $outbam = ($seqFileName eq 'fastq1') ? $out_bam_untr : $out_bam_tr;
    system("cp ".$directory."Aligned.sortedByCoord.out.bam ".$outbam);
    
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
		system($samtools_loc." sort -@ ". $cores. " -m 1000M ".$directory."Aligned.toTranscriptome.out.bam -o ".$directory."Aligned.toTranscriptome.out.sorted.bam 2>&1" );
		systemError("Samtools sort",$?,$!);

		#  convert BAM back to SAM file
		#print "converting BAM to SAM...\n";
		#system($samtools_loc." view -h -o ".$directory."Aligned.toTranscriptome.out.sorted.sam ".$directory."Aligned.toTranscriptome.out.sorted.bam > /dev/null 2>&1");
		#systemError("Samtools view",$?,$!);

		# rename SAM output file
		#print "renaming SAM output file...\n";

		# Bam file depends on what fastq file is processed (fastq1 = untreated, fastq2 = treated; that is for RIBO-seq experiments)
		my $bamf_tr = ($seqFileName  eq 'fastq1') ? $out_bam_tr_untr : $out_bam_tr_tr;
		system("mv ".$directory."Aligned.toTranscriptome.out.sorted.bam ".$bamf_tr);

		print "transcript coordinate $bamf_tr\n";
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
            my $PATH_TO_GTF   = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Annotation/Genes/genes_".$ensemblversion.".gtf";
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
    if (uc($readtype) eq "RIBO" || $readtype eq "ribo_untr" || uc($readtype) =~ m/SE/) {
        $fasta = $seqFile;
    } elsif (uc($readtype) =~ m/PE/) {
        $fasta1 = $seqFile;
        $fasta2 = $seqFile2;
    }

    # We need to first clip the adapter with fastx_clipper?
    if (uc($clipper) eq "FASTX") {

        print "     Clipping $seqFileName"." using fastx_clipper tool\n";
        # Without length cut-off and adaptor presence
        my $clip_command = $fastx_clip_loc." -Q33 -a ".$adaptorSeq." -l 20 -n –v -i ".$fasta." -o ".$work_dir."/fastq/$seqFileName"."_clip.fastq";
        print "     ".$clip_command."\n";
        system ($clip_command);
        $fasta = $work_dir."/fastq/$seqFileName"."_clip.fastq";
        print "clipfasta= $fasta\n";

        print "     Trimming $seqFileName"." using fastq_quality_trimmer tool\n";
        my $trim_command = $fastx_trim_loc." -Q33 -v -t 28 -l 20  -i ".$fasta." -o ".$work_dir."/fastq/$seqFileName"."_clip_trim.fastq";
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
    my $PATH_TO_GTF   = $IGENOMES_ROOT."/".$spec."/Ensembl/".$assembly."/Annotation/Genes/genes_".$ensemblversion.".gtf";
    my $PATH_TO_GENOME_INDEX = $ref_loc."".$IndexGenome;

    check_Bowtie_index($IndexGenome,$prev_mapper);

    # alignment dependent on read type
    # Tophat parameters:
    #   --max-multihits: only ouput that many hits if read maps to multiple locations, if more equally scoring options are available than set with this value, a limited number is randomly selected)
    #   --report-secondary-alignments: this parameter shouldn't be set, suboptimal (i.e. with lower score) are not retained
    if (uc($readtype) eq "RIBO" || $readtype eq "ribo_untr") {

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
    
    # # Output bamfile
    my $outbam = ($seqFileName eq 'fastq1') ? $out_bam_untr : $out_bam_tr;
    system("cp ".$directory."Aligned.sortedByCoord.out.bam ".$outbam);

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
                if (uc($readtype) eq 'RIBO' || $readtype eq "ribo_untr") {
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
    my $seqFileName1 = $_[13];
    my $seqFileName2 = $_[14];
    my $rpf_split = $_[15];
    my $FirstRankMultiMap = $_[16];
    my $truseq = $_[17];
    my $out_bg_s_untr = $_[18];
    my $out_bg_as_untr = $_[19];
    my $out_bg_s_tr = $_[20];
    my $out_bg_as_tr = $_[21];
    my $out_sam_untr = $_[22];
    my $out_sam_tr = $_[23];
    my $maxmultimap = $_[24];
    my $out_bam_untr = $_[25];
    my $out_bam_tr = $_[26];
    my $min_l_plastid = $_[27];
    my $max_l_plastid = $_[28];
    my $min_l_parsing = $_[29];
    my $max_l_parsing = $_[30];

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
    
    $query = "INSERT INTO arguments (variable,value) VALUES (\'seqFileName1\',\'".$seqFileName1."\')";
    $dbh_sqlite_results->do($query);
    
    $query = "INSERT INTO arguments (variable,value) VALUES (\'seqFileName2\',\'".$seqFileName2."\')";
    $dbh_sqlite_results->do($query);
    
    $query = "INSERT INTO arguments (variable,value) VALUES (\'rpf_split\',\'".$rpf_split."\')";
    $dbh_sqlite_results->do($query);
    
    $query = "INSERT INTO arguments (variable,value) VALUES (\'firstRankMultiMap\',\'".$FirstRankMultiMap."\')";
    $dbh_sqlite_results->do($query);
    
    $query = "INSERT INTO arguments (variable,value) VALUES (\'truseq\',\'".$truseq."\')";
    $dbh_sqlite_results->do($query);
    
    $query = "INSERT INTO arguments (variable,value) VALUES (\'out_bg_s_untr\',\'".$out_bg_s_untr."\')";
    $dbh_sqlite_results->do($query);
    
    $query = "INSERT INTO arguments (variable,value) VALUES (\'out_bg_as_untr\',\'".$out_bg_as_untr."\')";
    $dbh_sqlite_results->do($query);
    
    $query = "INSERT INTO arguments (variable,value) VALUES (\'out_bg_s_tr\',\'".$out_bg_s_tr."\')";
    $dbh_sqlite_results->do($query);
    
    $query = "INSERT INTO arguments (variable,value) VALUES (\'out_bg_as_tr\',\'".$out_bg_as_tr."\')";
    $dbh_sqlite_results->do($query);
    
    $query = "INSERT INTO arguments (variable,value) VALUES (\'out_sam_untr\',\'".$out_sam_untr."\')";
    $dbh_sqlite_results->do($query);
    
    $query = "INSERT INTO arguments (variable,value) VALUES (\'out_sam_tr\',\'".$out_sam_tr."\')";
    $dbh_sqlite_results->do($query);
    
    $query = "INSERT INTO arguments (variable,value) VALUES (\'maxmultimap\',\'".$maxmultimap."\')";
    $dbh_sqlite_results->do($query);
    
    $query = "INSERT INTO arguments (variable,value) VALUES (\'out_bam_untr\',\'".$out_bam_untr."\')";
    $dbh_sqlite_results->do($query);
    
    $query = "INSERT INTO arguments (variable,value) VALUES (\'out_bam_tr\',\'".$out_bam_tr."\')";
    $dbh_sqlite_results->do($query);
    
    $query = "INSERT INTO arguments (variable,value) VALUES (\'min_l_plastid\',\'".$min_l_plastid."\')";
    $dbh_sqlite_results->do($query);
    
    $query = "INSERT INTO arguments (variable,value) VALUES (\'max_l_plastid\',\'".$max_l_plastid."\')";
    $dbh_sqlite_results->do($query);
    
    $query = "INSERT INTO arguments (variable,value) VALUES (\'min_l_parsing\',\'".$min_l_parsing."\')";
    $dbh_sqlite_results->do($query);
    
    $query = "INSERT INTO arguments (variable,value) VALUES (\'max_l_parsing\',\'".$max_l_parsing."\')";
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
