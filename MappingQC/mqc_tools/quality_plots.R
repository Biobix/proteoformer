# This script makes plots of the gene distribution of the ribo-seq reads (dynamic range etc.)
###############

###########
## COMMAND LINE
## $ Rscript quality_plots.R in_table out_rankedgenes out_cumulative out_density
## **example: Rscript quality_plots.R mESC_GA_ens72_STAR_untreated_genedistribution.txt mESC_GA_ens72_STAR_untreated_rankedgenes.pdf mESC_GA_ens72_STAR_untreated_cumulative.pdf mESC_GA_ens72_STAR_untreated_density.pdf 
###########

# Parameters
args <- commandArgs(TRUE)
file<- as.character(args[1])
plotA <- as.character(args[2])
plotB <- as.character(args[3])
plotC <- as.character(args[4])

# Set outputfiles
#out <- substr(file,1,nchar(file)-4)

# Read in data gene distribution (column 1: genes, column 2: #reads)
input <- file
data <- read.table(input,header=TRUE,sep="\t",colClasses=c("character",NA))
data[,2] <- as.integer(data[,2])
colnames(data)[1] <- "geneID"
colnames(data)[2] <- "counts"

# Rank genes (descending)
sort <- order(-data$counts)
data_sort <- data[sort,]
genes <- length(data$geneID)
data_sort$rank <- seq(1,genes,1) #rank of gene 

# Normalize read counts
all <- sum(data$counts)
data_sort$counts_norm <- (data_sort$counts/all) #Normalized counts
data_sort$log2 <- log2(data_sort$counts) #Log2(#reads)
cumul <- cumsum(data_sort$counts_norm) 
data_sort$cumul <- cumul/cumul[genes] #cumulative #reads

# Plot transcript abundance of ranked genes
#plotA <- paste(out,"_ranked_genes.pdf",sep="")
#plotB <- paste(out,"_cumulative.pdf",sep="")
#plotC <- paste(out,"_density.pdf",sep="")
d <- density(data_sort$log2) #density log2(#reads)

png(plotA)
plot(data_sort$rank,data_sort$log2,main="Ranked gene abundance",xlab="Ranked genes",ylab="log2(#reads)",col="dark red")
dev.off()
png(plotB)
plot(data_sort$rank,data_sort$cumul,main="Cumulative abundance",xlab="Ranked genes",ylab="Cumulative #reads",col="dark red")
dev.off()
png(plotC)
plot(d,main="Density log2(#reads)",xlab="log2(#reads)",col="red")
dev.off()
