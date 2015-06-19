#!/usr/bin/perl -w
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use v5.10;
use Parallel::ForkManager;

## Command
## ./MAFtoPhyloCSFmaf.pl --cores 23 --MAF_root /data2/MAF/ --species human

# get the command line arguments
my ($cores,$maf_root,$IGENOMES_ROOT,$species);

GetOptions(
"cores=i"=>\$cores,                         # Number of cores to use for Bowtie Mapping,            mandatory argument
"maf_root=s" =>\$maf_root,                  # MAF ROOT FOLDER                                       mandatory argument
"species=s"=>\$species                      # Species, eg mouse/human/fruitfly,                     mandatory argument
);

if ($maf_root){
    print "The maf_root folder used is                              : $maf_root\n";
} else {
    die "\nDon't forget to pass the maf_root folder --maf_root argument!\n\n";
}
if ($cores){
    print "Number of cores to use for assembling                    : $cores\n";
} else {
    die "\nDon't forget to pass number of cores to use for mapping using the --cores or -c argument!\n\n";
}
if ($species){
    print "Species                                                  : $species\n";
} else {
    die "\nDon't forget to pass the Species name using the --species or -sp argument!\n\n";
}

#Conversion for species terminology
my $spec = ($species eq "mouse") ? "Mus_musculus" : ($species eq "human") ? "Homo_sapiens" : ($species eq "arabidopsis") ? "Arabidopsis_thaliana" : ($species eq "fruitfly") ? "Drosophila_melanogaster" : "";

############################################
##                                        ##
## PhyloCSF specific MAF files            ##
##                                        ##
############################################

# Start time
my $start = time;

print "\nGet chromosomes... \n";
## Get chromosomes based on seq_region_id ##
my @chrs = get_chrs($species);

# Init multi core
my $pm = new Parallel::ForkManager($cores);

print "\nExtract PhyloCSF multiple alignments from MAF... \n";
foreach(@chrs){
    
    ### Start parallel process
    $pm->start and next;
    
    ## Copy MAF to PhyloCSF maf
    MAFtoPhyloCSFmaf($_,$maf_root,$spec,$species);

    print "\nFinished Chromosome ".$_."\n";
    ### Finish childs
    $pm->finish;
}

#Waiting for all childs to finish
$pm->wait_all_children();

# End time
print "   DONE! \n";
my $end = time - $start;
printf("runtime assembly: %02d:%02d:%02d\n\n",int($end/3600), int(($end % 3600)/60), int($end % 60));

############
# THE SUBS #
############

### Maf to MAf ###

sub MAFtoPhyloCSFmaf{

    #Catch
    my $chr         =   $_[0];
    my $maf_root    =   $_[1];
    my $spec        =   $_[2];
    my $species     =   $_[3];
    
    #Init
    open INPUT, "<".$maf_root."/".$spec."/chr".$chr.".maf";
	open OUTPUT, "+>".$maf_root."/".$spec."/tmp_chr".$chr.".maf";
    
    
    my @NOT_PhyloCSF_species;
    if ($species eq "fruitfly"){
        @NOT_PhyloCSF_species =  ('anoGam1','apiMel3','triCas2');
    }elsif ($species eq "human"){
        @NOT_PhyloCSF_species =  ('gorGor3','ponAbe2','nomLeu3','macFas5','papHam1','chlSab1','calJac3','saiBol1','jacJac1','micOch1','criGri1','mesAur1','hetGla2','chiLan1','octDeg1','susScr3','camFer1','orcOrc1','panHod1','oviAri3','capHir1','cerSim1','musFur1','ailMel1','odoRosDiv1','lepWed1','pteAle1','myoDav1','eptFus1','conCri1','eleEdw1','triMan1','chrAsi1','oryAfe1','monDom5','sarHar1','macEug2','ornAna1','falChe1','falPer1','ficAlb2','zonAlb1','geoFor1','taeGut2','pseHum1','melUnd1','amaVit1','araMac1','colLiv1','anaPla1','galCal4','melGal1','allMis1','cheMyd1','chrPic1','pelSin1','apaSpi1','anoCar2','xenTro7','latCha1','tetNig2','fr3','takFla1','oreNil2','neoBri1','hapBur1','mayZeb1','punNye1','oryLat2','xipMac1','gasAcu1','gadMor1','danRer7','astMex1','lepOcu1','petMar2');
    }elsif ($species eq "mouse"){
        @NOT_PhyloCSF_species =  ('hetGla2','calJac3','gorGor3','nomLeu2','papHam1','ponAbe2','saiBol1','ailMel1','oviAri1','susScr3','triMan1','anoCar2','chrPic1','danRer7','fr3','gadMor1','galGal4','gasAcu1','latCha1','macEug2','melGal1','melUnd1','monDom5','oreNil2','ornAna1','oryLat2','petMar1','sarHar1','taeGut1','tetNig2','xenTro3');
    }

    #Copy files
    while ( my $line = <INPUT> ) {
        PRINT: {    foreach( @NOT_PhyloCSF_species){
                        if ($line =~ m/$_/){
                            last PRINT;
                        }
                    }
            print OUTPUT $line;
        }
    }
    close(INPUT);
    close(OUTPUT);
    
    #remove old maf mv tmp maf
    #system("rm -rf ".$maf_root."/".$spec."/chr".$chr.".maf");
    #system("mv ".$maf_root."/".$spec."/tmp_chr".$chr.".maf ".$maf_root."/".$spec."/chr".$chr.".maf");

}

### GET CHRs ###

sub get_chrs {
    
    # Catch
    my $species    =   $_[0];
    
    # Init
    my @chrs;
    
    # Get chrs from Chr_File
    if ($species eq "fruitfly"){
        @chrs   =   ('X','M','4','3R','3L','2R','2L');
    }elsif ($species eq "human"){
        @chrs   =   ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','M','X','Y');
    }elsif ($species eq "mouse"){
        @chrs   =   ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','MT','X','Y');
    }
	# Return
	return(@chrs);
    
}
