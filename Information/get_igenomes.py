__author__ = 'SV'

###############################################################################################
## Construct custom igenomes information structure for PROTEOFORMER based on data of Ensembl ##
###############################################################################################

###########
## USAGE ##
###########

'''
    ARGUMENTS:
    
    -d | --dir                                  Directory wherein the igenomes tree structure will be installed
    -v | --version                              Ensembl annotation version to download
                                                (Ensembl plant (for arabidopsis) has seperate annotation versions!)
    -s | --species                              Specify the desired species for which gene annotation files should be downloaded
    -r | --remove                               If any, overwrite the existing igenomes structure for that species
    -c | --cores                                The amount of cores that will be used for downloading chromosomes files
                                                (Do not use more than 15 cores as the download server can only establish 15 connections at once)
    
    currently supported species:
    
    human                       |   Homo_sapiens
    mouse                       |   Mus_musculus
    fruitfly                    |   Drosophila_melanogaster
    yeast                       |   Saccharomyces_cerevisiae
    zebrafish                   |   Danio_rerio
    arabidopsis                 |   Arabidopsis_thaliana
    c.elegans                   |   Caenorhabditis_elegans
    
    EXAMPLE:
    
    python get_igenomes.py -v 82 -s human -d /path/to/dir -r -c 15
    
    python get_igenomes.py -v 29 -s arabidopsis -d /path/to/dir -r -c 15
    
    DEPENDENCIES:
    
    This program depends upon wget, rsync and gzip, which are normally pre-installed on any unix system
    
    '''

import os
import shutil
import getopt
import sys
import re
from multiprocessing import Pool
import datetime

try:
    myopts, args = getopt.getopt(sys.argv[1:],"d:s:v:rc:",["dir=","version=","species=","remove","cores="])
except getopt.GetoptError as err:
    print(err)
    sys.exit()

#################################
# o == option                   #
# a == argument passed to the o #
#################################

#
# Catch arguments
#

removeExisting=False
for o, a in myopts:
    if o in ('-d','--dir'):
        instalDir=a;
    if o in ('-s','--species'):
        species=a
    if o in ('-v','--version'):
        ens_v=a
    if o in ('-r','--remove'):
        removeExisting=True
    if o in ('-c', '--cores'):
        stringCores=a
        cores=int(a)

#
# Check for correct arguments
#

if(instalDir == ''):
    print("Error: do not forget to pass the directory where the igenomes structure should be installed!")
    sys.exit()
if(species == ''):
    print("Error: do not forget to pass the species argument!")
    sys.exit()
if(ens_v == ''):
    print("Error: do not forget to pass the ensembl version argument!")
    sys.exit()
if(cores == ''):
    print("Error: do no forget to pass the amount of cores!")
    sys.exit()
elif(cores>15):
    print("Error: the amount of cores cannot be larger than 15!")
    sys.exit()
#Remove last "/" from instal dir path
pattern=re.compile('^(\S+)/$')
m = pattern.match(instalDir)
if(m):
    instalDir=m.group(1)
#Report on input
print("The igenomes structure will be installed in : " + instalDir)
print("Ensembl version used                        : " + ens_v)
print("Selected species                            : " + species)
print("Amount of cores                             : " + stringCores)
print("")

#Convert species and construct additional info. New assemblies can be modified here.
if(species=='human'):
    speciesLong='Homo_sapiens'
    assembly='GRCh38'
    ucscCode='hg38'
elif(species=='mouse'):
    speciesLong='Mus_musculus'
    assembly='GRCm38'
    ucscCode='mm10'
elif(species=='fruitfly'):
    speciesLong='Drosophila_melanogaster'
    assembly='BDGP6'
    ucscCode='dm6'
elif(species=='yeast'):
    speciesLong='Saccharomyces_cerevisiae'
    assembly='R64-1-1'
    ucscCode='sacCer3'
elif(species=='zebrafish'):
    speciesLong='Danio_rerio'
    assembly='GRCz10'
    ucscCode='danRer10'
elif(species=='arabidopsis'):
    speciesLong='Arabidopsis_thaliana'
    assembly='TAIR10'
elif(species=='c.elegans'):
    speciesLong='Caenorhabditis_elegans'
    assembly='WBcel235'
    ucscCode='ce10'
else:
    print("Species has to be one of the following list: human, mouse, fruitfly, yeast, zebrafish, arabidopsis, c.elegans")
    sys.exit()


os.chdir(instalDir)

#Check if the igenomes folder already exists
if os.path.isdir("igenomes"):
    if os.path.isdir("igenomes/"+speciesLong):
        if(removeExisting==False):
            print("There is already folder called igenomes/"+speciesLong+" for "+species+" in "+instalDir)
            print("If you want to overwrite the existing structure for that species in "+instalDir+", please use the -r or --remove option.")
            sys.exit()
        else:
            print("Overwriting existing igenomes folder for "+species)
            shutil.rmtree("igenomes/"+speciesLong)
else:
    os.system("mkdir igenomes")



#construct the basic folder structure
os.system("mkdir igenomes/"+speciesLong)
os.system("mkdir igenomes/"+speciesLong+"/Ensembl")
os.system("mkdir igenomes/"+speciesLong+"/Ensembl/"+assembly)
os.system("mkdir igenomes/"+speciesLong+"/Ensembl/"+assembly+"/Annotation")
os.system("mkdir igenomes/"+speciesLong+"/Ensembl/"+assembly+"/Sequence")

os.chdir(instalDir+"/igenomes/"+speciesLong+"/Ensembl/"+assembly+"/Annotation")



## Chr sizes files
print("")
print("Get Chromosome lengths from UCSC")
chromList = {}
os.system("mkdir Genes")
os.chdir(instalDir+"/igenomes/"+speciesLong+"/Ensembl/"+assembly+"/Annotation/Genes")

#Fetch data from UCSC except for arabidopsis (not included in UCSC)
if(species=='arabidopsis'):
    #Manual from custom arabidopsis UCSC site: http://epigenomics.mcdb.ucla.edu/cgi-bin/hgTracks?hgsid=25975&chromInfoPage= (no txt file available)
    chromList['1']='30432563'
    chromList['2']='19705359'
    chromList['3']='23470805'
    chromList['4']='18585042'
    chromList['5']='26992728'
    chromList['Pt']='154478'
    chromList['M']='366923'
else:
    #Other species: download from UCSC
    os.system("wget -q ftp://hgdownload.cse.ucsc.edu/goldenPath/"+ucscCode+"/database/chromInfo.txt.gz")
    os.system("gzip -d chromInfo.txt.gz")
    os.system("mv chromInfo.txt tmpChromInfo.txt")
    #Retain only the standard chromosomes
    with open('tmpChromInfo.txt') as inFile:
        for line in inFile:
            pattern = re.compile('^chr(\w{1,4})\t(\d+)\t\S+\n$') #chrVIII (yeast) is the longest one to capture
            m = pattern.match(line)
            if m:
                chromList[m.group(1)]=m.group(2)
#Write standard chromosomes to chromosome sizes file
outFile = open('ChromInfo.txt','w')
for key in chromList:
    outFile.write(key+"\t"+chromList[key]+"\n")
outFile.close()
os.system("rm -rf tmpChromInfo.txt")




## Download chromosome files
print("\n")
print("Get chromosome fasta files")
os.chdir(instalDir+"/igenomes/"+speciesLong+"/Ensembl/"+assembly+"/Sequence")
os.system("mkdir Chromosomes")
os.chdir(instalDir+"/igenomes/"+speciesLong+"/Ensembl/"+assembly+"/Sequence/Chromosomes")

#Defenition of one download process
def downloadChromosomeFasta(chr):
    if(species=='arabidopsis'):#Arabidopsis is on the site of ensembl Plants instead of normal Ensembl. This site cannot use rsync yet.
        if(chr=='M'):
            os.system("wget -q ftp://ftp.ensemblgenomes.org/pub/plants/release-"+ens_v+"/fasta/"+speciesLong.lower()+"/dna//"+speciesLong+"."+assembly+"."+ens_v+".dna.chromosome.Mt.fa.gz")
            os.system("gunzip "+speciesLong+"."+assembly+"."+ens_v+".dna.chromosome.Mt.fa.gz")
            os.system("mv "+speciesLong+"."+assembly+"."+ens_v+".dna.chromosome.Mt.fa Mt.fa")
        else:
            os.system("wget -q ftp://ftp.ensemblgenomes.org/pub/release-"+ens_v+"/plants/fasta/"+speciesLong.lower()+"/dna//"+speciesLong+"."+assembly+"."+ens_v+".dna.chromosome."+chr+".fa.gz")
            os.system("gunzip "+speciesLong+"."+assembly+"."+ens_v+".dna.chromosome."+chr+".fa.gz")
            os.system("mv "+speciesLong+"."+assembly+"."+ens_v+".dna.chromosome."+chr+".fa "+chr+".fa")
        print("\t\t*) Chromosome "+chr+" finished")
    else:#use rsync for other species
        if(chr=='M'): #Ensembl uses MT instead of M for mitochondrial genome
            if(species=='fruitfly'):#Other name 'dmel_mitochondrion_genome' for fruitfly
                os.system("rsync -avq rsync://ftp.ensembl.org/ensembl/pub/release-"+ens_v+"/fasta/"+speciesLong.lower()+"/dna//"+speciesLong+"."+assembly+".dna.chromosome.dmel_mitochondrion_genome.fa.gz "+instalDir+"/igenomes/"+speciesLong+"/Ensembl/"+assembly+"/Sequence/Chromosomes/M.fa.gz") #Give name 'M' instead of 'MT' like in real igenomes bundles
                os.system("gunzip M.fa.gz")
            elif(species=='yeast'):#Other name 'Mito' for yeast
                os.system("rsync -avq rsync://ftp.ensembl.org/ensembl/pub/release-"+ens_v+"/fasta/"+speciesLong.lower()+"/dna//"+speciesLong+"."+assembly+".dna.chromosome.Mito.fa.gz "+instalDir+"/igenomes/"+speciesLong+"/Ensembl/"+assembly+"/Sequence/Chromosomes/MT.fa.gz")
                os.system("gunzip MT.fa.gz")
            elif(species=='c.elegans'):#Other name 'MtDNA' for c elegans
                os.system("rsync -avq rsync://ftp.ensembl.org/ensembl/pub/release-"+ens_v+"/fasta/"+speciesLong.lower()+"/dna//"+speciesLong+"."+assembly+".dna.chromosome.MtDNA.fa.gz "+instalDir+"/igenomes/"+speciesLong+"/Ensembl/"+assembly+"/Sequence/Chromosomes/MT.fa.gz")
                os.system("gunzip MT.fa.gz")
            else:#MT for other species
                os.system("rsync -avq rsync://ftp.ensembl.org/ensembl/pub/release-"+ens_v+"/fasta/"+speciesLong.lower()+"/dna//"+speciesLong+"."+assembly+".dna.chromosome.MT.fa.gz "+instalDir+"/igenomes/"+speciesLong+"/Ensembl/"+assembly+"/Sequence/Chromosomes/MT.fa.gz")
                os.system("gunzip MT.fa.gz")
        else:
            os.system("rsync -avq rsync://ftp.ensembl.org/ensembl/pub/release-"+ens_v+"/fasta/"+speciesLong.lower()+"/dna//"+speciesLong+"."+assembly+".dna.chromosome."+chr+".fa.gz "+instalDir+"/igenomes/"+speciesLong+"/Ensembl/"+assembly+"/Sequence/Chromosomes/"+chr+".fa.gz")
            os.system("gunzip "+chr+".fa.gz")
        print("\t\t*) Chromosome "+chr+" finished")

#Multithreading chromosome downloads
pool = Pool(processes=cores)
[pool.apply_async(downloadChromosomeFasta, args=(key,)) for key in chromList]
pool.close()
pool.join()





## Make whole genome fasta
print("\n")
print("Make whole genome fasta file as a concatenation of the seperate chromosome files")
os.chdir(instalDir+"/igenomes/"+speciesLong+"/Ensembl/"+assembly+"/Sequence")
os.system("mkdir WholeGenomeFasta")
os.chdir(instalDir+"/igenomes/"+speciesLong+"/Ensembl/"+assembly+"/Sequence/WholeGenomeFasta")

#Whole genome fasta file is a concatenation of all chromosome fasta files: first, construct command
command = ""
for chr in chromList:
    if(chr=='M'):#Different names for mitochondrial chromosome file
        if(species=='arabidopsis'):
            command=command+" "+instalDir+"/igenomes/"+speciesLong+"/Ensembl/"+assembly+"/Sequence/Chromosomes/Mt.fa"
        elif(species=='fruitfly'):
            command=command+" "+instalDir+"/igenomes/"+speciesLong+"/Ensembl/"+assembly+"/Sequence/Chromosomes/M.fa"
        else:
            command=command+" "+instalDir+"/igenomes/"+speciesLong+"/Ensembl/"+assembly+"/Sequence/Chromosomes/MT.fa"
    else:
        command=command+" "+instalDir+"/igenomes/"+speciesLong+"/Ensembl/"+assembly+"/Sequence/Chromosomes/"+chr+".fa"
#Do concatenation
os.system("cat"+command+"> genome.fa")






## Download genes.gtf file
print("\n")
print("Download genes.gtf file")
os.chdir(instalDir+"/igenomes/"+speciesLong+"/Ensembl/"+assembly+"/Annotation/Genes")
if(species=='arabidopsis'):#Arabidopsis from Ensembl Plants
    os.system("wget -q ftp://ftp.ensemblgenomes.org/pub/release-"+ens_v+"/plants/gtf/"+speciesLong.lower()+"//"+speciesLong+"."+assembly+"."+ens_v+".gtf.gz")
    os.system("mv "+speciesLong+"."+assembly+"."+ens_v+".gtf.gz genesTmp.gtf.gz")
else:
    os.system("rsync -avq rsync://ftp.ensembl.org/ensembl/pub/release-"+ens_v+"/gtf/"+speciesLong.lower()+"//"+speciesLong+"."+assembly+"."+ens_v+".gtf.gz genesTmp.gtf.gz")
os.system("gunzip genesTmp.gtf.gz")

#The first lines are comments and are unwanted: delete them
genesFile = open('genes.gtf','w')
with open('genesTmp.gtf') as genesTmp:
    for line in genesTmp:
        pattern = re.compile('^#!')
        m = pattern.match(line)
        if(m==None):
            genesFile.write(line)
genesFile.close()
os.system("rm -rf genesTmp.gtf")



##Download supplemental abundant sequences
print("\n")
print("PhiX fasta file downloading in "+instalDir+"igenomes/"+speciesLong+"/Ensembl/"+assembly+"/Sequence/AbundantSequences. Other abundant sequences (e.g. rRNA) can be added in this folder too.")
os.chdir(instalDir+"/igenomes/"+speciesLong+"/Ensembl/"+assembly+"/Sequence")
os.system("mkdir AbundantSequences")
os.chdir(instalDir+"/igenomes/"+speciesLong+"/Ensembl/"+assembly+"/Sequence/AbundantSequences")
os.system("wget -q http://bcb.dfci.harvard.edu/~vwang/phix.fasta")
os.system("mv phix.fasta phix.fa")





## Make README.txt file
os.chdir(instalDir+"/igenomes/"+speciesLong+"/Ensembl/"+assembly+"/Annotation")
readmeFile = open('README.txt','w')
#Get date
downloadDate = datetime.date.today()
year=downloadDate.year
monthinteger=downloadDate.month
month=datetime.date(1900, monthinteger, 1).strftime('%B')
day=downloadDate.day
readmeFile.write("The contents of the annotation directories were downloaded from Ensembl on: "+month+" "+str(day)+", "+str(year)+".\n")
if(species=='arabidopsis'):
    readmeFile.write("Gene annotation files were downloaded from Ensembl Plants release "+ens_v+".")
else:
    readmeFile.write("Gene annotation files were downloaded from Ensembl release "+ens_v+".")
readmeFile.close()






##Change permissions for free consultation
os.chdir(instalDir)
os.system("chmod 777 igenomes")
os.system("chgrp biobix igenomes")
os.chdir(instalDir+"/igenomes")
os.system("chmod -R 777 "+speciesLong)
os.system("chgrp -R biobix "+speciesLong)

print("\n")
print("   (***) igenomes folder download complete (***)")
print("\n")
