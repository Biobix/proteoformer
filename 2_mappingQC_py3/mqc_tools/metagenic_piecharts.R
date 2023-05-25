# This script makes a functional annotation pie chart of the ribo-seq reads that map in ensembl transcripts
###############

###########
## COMMAND LINE
## $ Rscript metagenic_piecharts.R in_coding in_noncoding out_coding out_noncoding
## **example: Rscript metagenic_piecharts.R mESC_GA_ens72_STAR_untreated_annotation_coding.txt mESC_GA_ens72_STAR_untreated_annotation_noncoding.txt mESC_GA_ens72_STAR_untreated_annotation_coding.pdf mESC_GA_ens72_STAR_untreated_annotation_noncoding.pdf
###########

# Parameters
args <- commandArgs(TRUE)
file1 <- as.character(args[1])
file2 <- as.character(args[2])
out_c <- as.character(args[3])
out_nc <- as.character(args[4])

# Set outputfiles
#out_c <- paste(substr(file1,1,nchar(file1)-4),".pdf",sep="")
#out_nc <- paste(substr(file2,1,nchar(file2)-4),".pdf",sep="")

# PIE CHART FOR READS IN PROTEIN (=coding) TRANSCRIPTS
#########
# Read in data
inputfile_c <- file1
data_c <- read.table(inputfile_c,header=TRUE,sep="\t");
#column1: chr
#column2: ribo_reads
#column3: exon_reads
#column4: 5'utr_reads
#column5 3'utr_reads
#column6: intron_reads
#column7: non protein_coding biotype transcripts
#column7: intergenic

# Calculate relative amount of reads in exons, 5'utr and 3'utr
ribo_reads_c = colSums(data_c["ribo"])
exon_reads_c = colSums(data_c["exon"])
utr5_reads_c = colSums(data_c["X5utr"])
utr3_reads_c = colSums(data_c["X3utr"])
intron_reads_c = colSums(data_c["intron"])
nPC_reads_c = colSums(data_c["non_protein_coding"])
intergenic_reads_c = colSums(data_c["intergenic"])

# Make functional region annotation pie chart
slices_c <- c(intergenic_reads_c, exon_reads_c, utr5_reads_c, utr3_reads_c, intron_reads_c,nPC_reads_c) 
pct_c <- round(slices_c/sum(slices_c)*100,2) 
lbls_c <- paste(pct_c,"%",sep="") # add % to labels 
png(file=out_c, width=2600, height=2600, res=384)
colors <- c("orangered2","olivedrab1","limegreen","blue","orchid4","yellow")
pie(slices_c,labels = lbls_c, col=colors,radius=0.7,main="Metagenic classification",cex=0.8)
legend("topleft", c("Intergenic", "Coding region", "5'UTR", "3'UTR", "Intron", "Other biotypes"), cex=0.7, fill=colors)
dev.off()

# PIE CHART FOR READS IN NON-PROTEIN (=non-coding) TRANSCRIPTS
#########
# Read in data
inputfile_nc <- file2
data_nc <- read.table(inputfile_nc,header=TRUE,sep="\t");
#column1: chr                             
#column2: ribo_reads nPC biotypes
#column3..N: nPC biotypes

# Calculate relative amount of reads in exons, 5'utr and 3'utr
#ribo_reads_nc = colSums(data_nc["non_protein_coding"])
biotypes <- rep(0,dim(data_nc)[2]-2)
for(i in 1:length(biotypes)){
  biotypes[i] <- colSums(data_nc[names(data_nc)[i+2]])
  names(biotypes)[i] <- names(data_nc)[i+2]
}
biotypes_covered <- biotypes[which(biotypes>0)]
minors <- biotypes_covered[which(round(biotypes_covered/sum(biotypes_covered)*100,2)<1)]
majors <- biotypes_covered[which(round(biotypes_covered/sum(biotypes_covered)*100,2)>=1)]

# Make functional region annotation pie chart
slices_nc <- c(sum(minors),majors)
names(slices_nc) <- c("Others",names(majors))
pct_nc <- round(slices_nc/sum(slices_nc)*100,2)
lbls_nc <- paste(pct_nc,"%",sep="") # add % to labels 
png(file=out_nc, width=2600, height=2600, res=384)
colors <- c("yellow",rainbow(length(slices_nc)-1))
pie(slices_nc,labels = lbls_nc, col=colors,radius=0.7,main="Overview other Ensembl biotypes",cex=0.8)
legend("topleft", names(slices_nc), cex=0.5, fill=colors)
dev.off()
