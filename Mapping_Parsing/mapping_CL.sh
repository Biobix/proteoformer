#!/bin/bash
#$ -S /bin/bash

# Run as:  qsub -cwd -pe threaded 20-28 mapping.sh

module load bowtie
module load bowtie2
module load star 
module load tophat
module load igenomes
module load bedtools
module load samtools
module load picard-tools
module load fastqc

#Variables for SGE environment
echo $TMP
echo $NSLOTS

$HOME/HCT116/1_mapping_v2.pl --dir $HOME/HCT116/ --file IngHar --name mESC --species mouse --ensembl 66 --cores $NSLOTS --readtype ribo --mapper STAR --unique N --tmp $TMP
