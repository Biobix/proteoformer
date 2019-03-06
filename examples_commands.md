# PROTEOFORMER 2.0: Commands used in the examples

Author: Steven Verbruggen (March, 2019)

In this file, I want to give an overview of all commands used to generate the results for the HCT116 and Jurkat cell lines, described in the PROTEOFORMER 2.0 paper. 
These commands are intended to be used as an example for other proteogenomic studies and as a general illustration of the PROTEOFORMER 2.0 pipeline.
All used scripts can be found in this GitHub directory.

## HCT116

### Ribosome profiling data download

Data of SRR1333393 up until SRR1333394 can be downloaded with the download_sra_parallel tool:

```
source activate download_sra_parallel
./download_sra_parallel.sh -f 1333393 -l 1333394 –c 25
```

For most other steps, the proteoformer Conda environment is used:

```
source activate proteoformer
```

### FastQC on raw files

Check the quality of the raw FASTQ files.

```
fastqc HCT116_CHX.fastq -t 10
fastqc HCT116_LTM.fastq -t 10
```

### Mapping

Reference information includes reference sequences and Ensembl databases.
These can be installed as follows:
```
python get_igenomes.py -v 92 -s human -d /data/igenomes -r -c 15
python ENS_db.py -v 92 -s human
```
Extra reference sequences (rRNA, tRNA and sn(-o-)RNA) can be downloaded from [NCBI](https://www.ncbi.nlm.nih.gov/nuccore) , [BioMart](https://www.ensembl.org/biomart/martview) or [GtRNAdb](http://gtrnadb.ucsc.edu/).
These sequences were downloaded as FASTA and added to the `./igenomes/Homo_sapiens/Ensembl/GRCh38/Sequence/abundant_sequences` folder of the newly installed iGenomes directory.

Then, we do the mapping (with STAR). This includes different subsequent steps:
* clipping off adapters (AGATCGGAAGAGCACAC)
* filtering out PhiX bacteriophage sequences
* filtering out rRNA, tRNA and sn(-o-)RNA
* mapping against the human genome non-uniquely
* calculating Plastid offsets
* parsing the alignments into output files: count tables in SQLite, BedGraph files,...

```
perl ./mapping.pl --inputfile1 HCT116_CHX.fastq --inputfile2 HCT116_LTM.fastq --name HCT116_nonuniq --species human --ensembl 92 --cores 20 --readtype ribo --mapper STAR --readlength 36 --adaptor AGATCGGAAGAGCACAC --unique N --igenomes_root /data/igenomes/ --clipper fastx --phix Y --rRNA Y --snRNA Y --tRNA Y --rpf_split Y --price_files Y --suite plastid
```

### Quality control of the alignments

This includes performing both FastQC and mQC on both alignment files:

```
fastqc STAR/fastq1/untreat.sam
fastqc STAR/fastq2/treat.sam
./mappingQC.pl --samfile STAR/fastq1/untreat.sam --cores 20 --treated untreated --result_db SQLite/results.db --unique N --ens_db ENS_hsa_92.db --offset plastid --offset_img plastid/HCT116_nonuniq_untreated_p_offsets.png --output_folder mappingQC_untr_output --plotrpftool pyplot3D
./mappingQC.pl --samfile STAR/fastq2/treat.sam --cores 20 --treated treated --result_db SQLite/results.db --unique N --ens_db ENS_hsa_92.db --offset plastid --offset_img plastid/HCT116_nonuniq_treated_p_offsets.png --output_folder mappingQC_tr_output --plotrpftool pyplot3D
```
After visually checking all quality information in a HTML browser, we continued with the pipeline.

### Transcript calling

Actively translated transcripts were called (rule-based):

```
perl ribo_translation.pl --in_sqlite SQLite/results.db --out_sqlite SQLite/results.db --ens_db ENS_hsa_92.db 
```

### TIS ID 1: classic proteoform calling

#### TIS calling

Calling translation initiation sites based on default values (optimal values as determined in the first PROTEOFORMER manuscript (Crappé et al. 2014)):

```
perl TIScalling_categorised.pl --sqlite_db SQLite/results.db 
```

#### Assembly

Assembly of ORF products:

```
perl assembly.pl --sqliteRES SQLite/results.db --tis_ids 1
```

### TIS ID 2: PRICE

Run the PRICE wrapper:

```
python PRICE.py -r SQLite/results.db 
```

### TIS ID 3: SPECtre

Run the SPECtre wrapper:

```
source activate spectre
python SPECtre.py -r SQLite/results.db -b STAR/fastq1/untreat.bam --offsets 28:12,29:12,30:12 -c 60 -x 3
source deactivate
source activate proteoformer
```

SPECtre can only work with a certain number of offsets, so we used the offsets of the main RPF lengths.

### FASTA file generation

Generate non-redundant FASTA files for all TIS IDs:

```
perl generate_translation_db.pl --result_db SQLite/results.db --mflag 4 --tis_ids 1
perl generate_translation_db.pl --result_db SQLite/results.db --mflag 4 --tis_ids 2 
perl generate_translation_db.pl --result_db SQLite/results.db --mflag 4 --tis_ids 3 
```

Combine these non-redundant files:

```
python combine_dbs.py -i human_TIS_1_transcripts.fasta,human_TIS_2_transcripts.fasta,human_TIS_3_transcripts.fasta -o HCT116_nonuniq_combined_output.fa -t fasta_combination_overview.txt --venn_file venn_diagram_nonredundant.png
```

Also, generate the redundant FASTA files and combine these:

```
perl generate_translation_db.pl --result_db SQLite/results.db --tis_ids 1 --mflag 5 --translation_db human_TIS_1_transcripts_redundant.fasta --var_file human_TIS_1_transcripts_VAR_redundant.txt
perl generate_translation_db.pl --result_db SQLite/results.db --tis_ids 2 --mflag 5 --translation_db human_TIS_2_transcripts_redundant.fasta --var_file human_TIS_2_transcripts_VAR_redundant.txt
perl generate_translation_db.pl --result_db SQLite/results.db --tis_ids 3 --mflag 5 --translation_db human_TIS_3_transcripts_redundant.fasta --var_file human_TIS_3_transcripts_VAR_redundant.txt
python combine_dbs.py --in_files human_TIS_1_transcripts_redundant.fasta,human_TIS_2_transcripts_redundant.fasta,human_TIS_3_transcripts_redundant.fasta --output_file HCT116_nonuniq_combined_output_redundant.fa --overview_file fasta_combination_overview_redundant.txt --venn_file venn_diagram_redundant.png
```

Combined FASTA files can be converted to SQLite databases for easy querying:

```
python combfasta2sqlite.py --fasta HCT116_nonuniq_combined_output.fa
python combfasta2sqlite.py --fasta HCT116_nonuniq_combined_output_redundant.fa
```

### Combination with UniProt

Human UniProt database were downloaded from the [UniProt site](https://www.uniprot.org/) in FASTA format.
The combination was both performed with the canonical and the splice-included version of UniProt.

For the canonical version:

```
python combine_with_uniprot.py –f HCT116_nonuniq_combined_output_redundant.fa –u uniprotHumanProteomeCanonical.fasta
python combuniprot2sqlite.py --fasta comb_fasta_uniprot_proteoformer.fasta
```

For the splice-included version:

```
python combine_with_uniprot.py -f HCT116_nonuniq_combined_output_redundant.fa -u uniprotHumanProteomeIsoform.fasta
```

### MS validation with MaxQuant

For this validation we used MaxQuant version 1.6.1.0. MaxQuant is used as graphical user interface.
One of the different combined FASTA files was used as search space. Raw data can be found on the PRIDE repository under identifier PXD011353.
Most of the parameters were chosen as default, except the following:
* Methionine oxidation as fixed modification
* N-term acetylation as variable modification
* Match between runs
* MaxLFQ: True
* Only unique peptides for protein quantification, no razor peptides
* iBAQ: True

For each run, from the `combined/txt` MaxQuant folder, we used the following files:
* proteinGroups.txt
* peptides.txt
* msms.txt

### MS results parsing

For the combination files of the different proteoform calling techniques with removed redundancy:

```
python parse_maxquant.py -g proteinGroups.txt -p peptides.txt -x msms.txt -t fasta_combination_overview.txt -f HCT116_nonuniq_combined_output.fa -e ENS_hsa_92.db -r N -l protein
python parse_maxquant.py -g proteinGroups.txt -p peptides.txt -x msms.txt -t fasta_combination_overview.txt -f HCT116_nonuniq_combined_output.fa -e ENS_hsa_92.db -r N -l peptide
python parse_maxquant.py -g proteinGroups.txt -p peptides.txt -x msms.txt -t fasta_combination_overview.txt -f HCT116_nonuniq_combined_output.fa -e ENS_hsa_92.db -r N -l PSM
```
For the combination files of the different proteoform calling techniques with redundancy:

```
python parse_maxquant.py -g proteinGroups.txt -p peptides.txt -x msms.txt -t fasta_combination_overview.txt -f HCT116_nonuniq_combined_output_redundant.fa -e ENS_hsa_92.db -r Y -l protein
python parse_maxquant.py -g proteinGroups.txt -p peptides.txt -x msms.txt -t fasta_combination_overview.txt -f HCT116_nonuniq_combined_output_redundant.fa -e ENS_hsa_92.db -r Y -l peptide
python parse_maxquant.py -g proteinGroups.txt -p peptides.txt -x msms.txt -t fasta_combination_overview.txt -f HCT116_nonuniq_combined_output_redundant.fa -e ENS_hsa_92.db -r Y -l PSM
```
For the combination with UniProt canonical:

```
python parse_maxquant_uniprot.py -g proteinGroups.txt -p peptides.txt -x msms.txt -f comb_fasta_uniprot_proteoformer.fasta -l protein -e ENS_hsa_92.db
python parse_maxquant_uniprot.py -g proteinGroups.txt -p peptides.txt -x msms.txt -f comb_fasta_uniprot_proteoformer.fasta -l peptide -e ENS_hsa_92.db
python parse_maxquant_uniprot.py -g proteinGroups.txt -p peptides.txt -x msms.txt -f comb_fasta_uniprot_proteoformer.fasta -l PSM -e ENS_hsa_92.db
python analyse_proteoforms.py -m max_proteins_uniprot.db -e ENS_hsa_92.db -M HUMAN_9606_idmapping_selected.tab -f ../comb_fasta_uniprot_proteoformer.fasta -g proteinGroups.txt -p peptides.txt -c data/tools/clustalo/clustalo
```
(For this last tool, an ID mapped overview file is needed, available from the [UniProt FTP server](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz))

For the combination with UniProt spliced:

```
python parse_maxquant_uniprot.py -g proteinGroups.txt -p peptides.txt -x msms.txt -f comb_fasta_uniprotspliced_proteoformer.fasta -l protein -e ENS_hsa_92.db
python parse_maxquant_uniprot.py -g proteinGroups.txt -p peptides.txt -x msms.txt -f comb_fasta_uniprotspliced_proteoformer.fasta -l peptide -e ENS_hsa_92.db
python parse_maxquant_uniprot.py -g proteinGroups.txt -p peptides.txt -x msms.txt -f comb_fasta_uniprotspliced_proteoformer.fasta -l PSM -e ENS_hsa_92.db
python analyse_proteoforms.py -m max_proteins_uniprot.db -e ENS_hsa_92.db -M HUMAN_9606_idmapping_selected.tab -f ../comb_fasta_uniprotspliced_proteoformer.fasta -g proteinGroups.txt -p peptides.txt -c data/tools/clustalo/clustalo
```


## Jurkat

### Ribosome profiling data download

Data of NCBI project GSE74279 needs to be downloaded:

```
source activate download_sra_parallel
./download_sra_parallel.sh -f SRR2732970 -l SRR2732970 –c 25
./download_sra_parallel.sh -f SRR2733100 -l SRR2733100 –c 25
```
### FastQC on raw files

Check the quality of the raw FASTQ files.

```
fastqc -o fastqc_chx/ -t 20 raw/jurkat_chx.fastq 
fastqc -o fastqc_ltm/ -t 20 raw/jurkat_ltm.fastq
```
### Mapping

Reference information includes reference sequences and Ensembl databases.
These can be installed as follows:
```
python get_igenomes.py -v 92 -s human -d /data/igenomes -r -c 15
python ENS_db.py -v 92 -s human
```
Extra reference sequences (rRNA, tRNA and sn(-o)RNA) can be downloaded from [NCBI](https://www.ncbi.nlm.nih.gov/nuccore) , [BioMart](https://www.ensembl.org/biomart/martview) or [GtRNAdb](http://gtrnadb.ucsc.edu/).
These sequences were downloaded as FASTA and added to the `./igenomes/Homo_sapiens/Ensembl/GRCh38/Sequence/abundant_sequences` folder of the newly installed iGenomes directory.

Then, we do the mapping (STAR). This includes different subsequent steps:
* clipping off adapters (CTGTAGGCACCATCAATAGATCGGAAGAGCACAC)
* filtering out PhiX bacteriophage sequences
* filtering out rRNA, tRNA and sn(-o)RNA
* mapping against the human genome uniquely
* calculating Plastid offsets
* parsing the alignments into output files: count tables in SQLite, BedGraph files,...

```
perl mapping.pl --inputfile1 ../raw/jurkat_chx.fastq --inputfile2 ../raw/jurkat_ltm.fastq --name jurkat_uniq --species human --ensembl 92 --cores 20 --readtype ribo --readlength 36 --mapper STAR --adaptor CTGTAGGCACCATCAATAGATCGGAAGAGCACAC --unique Y --igenomes_root /data/igenomes/ --clipper fastx --phix Y --rRNA Y --snRNA Y --tRNA Y --rpf_split Y --suite plastid --price_files Y
```
 
### Quality control of the alignments

This includes performing FastQC on both alignment files:

```
fastqc -o fastqc_untreat/ STAR/fastq1/untreat.sam 
fastqc -o fastqc_treat/ STAR/fastq2/treat.sam 
```
Visual inspection of the FastQC alignment HTML reports, resulted in a small rest fraction of low quality reads >=37bp.
These reads were not present in high quantities but we decided to filter them out anyway:

```
#Filter for read lengths < 37bp
nohup awk 'length($10) < 37 || $1 ~ /^@/' untreat.sam > untreat_filt.sam &
nohup awk 'length($10) < 37 || $1 ~ /^@/' treat.sam > treat_filt.sam &
#Rename filtered files as untreat.sam and treat.sam for the rest of the pipeline
mv untreat_filt.sam untreat.sam
mv treat_filt.sam treat.sam
#samtools view to convert them to resp bam files
samtools view -b untreat.bam untreat.sam
samtools view -b treat.bam treat.sam
```

Redo Plastid and mapping parsing:

```
perl mapping_plastid.pl --out_sqlite SQLite/results.db --treated untreated
perl mapping_plastid.pl --out_sqlite SQLite/results.db --treated treated
perl mapping_parsing.pl --out_sqlite SQLite/results.db --offset plastid
```

Redo quality control on filtered alignment files:

```
fastqc -o fastqc_untreat/ STAR/fastq1/untreat.sam 
fastqc -o fastqc_treat/ STAR/fastq2/treat.sam
perl mappingQC.pl --samfile STAR/fastq1/untreat.sam --treated untreated --cores 25 --result_db SQLite/results.db --unique Y --ens_db ENS_hsa_92.db --offset plastid --offset_img plastid/jurkat_uniq_untreated_p_offsets.png --tool_dir mqc_tools/ --plotrpftool pyplot3D
perl mappingQC.pl --samfile STAR/fastq2/treat.sam --treated treated --cores 25 --result_db SQLite/results.db --unique Y --ens_db ENS_hsa_92.db --offset plastid --offset_img plastid/jurkat_uniq_treated_p_offsets.png --tool_dir mqc_tools/ --plotrpftool pyplot3D 
```

After visually checking all quality information in a HTML browser, we continued with the pipeline.

### Transcript calling

Actively translated transcripts were called (rule-based):

```
perl ribo_translation.pl --in_sqlite SQLite/results.db --out_sqlite SQLite/results.db --ens_db ENS_hsa_92.db 
```

### TIS ID 1: classic proteoform calling

#### TIS calling

Calling translation initiation sites based on default values (optimal values as determined in the first PROTEOFORMER manuscript (Crappé et al. 2014)):

```
perl TIScalling_categorised.pl --sqlite_db SQLite/results.db 
```

#### SNP calling

Get a SNP VCF file from the NCBI FTP server. More information about where to find this can be found in the readme in the SNP calling folder of this GitHub Repo.

```
ftp> get common_all_20180418.vcf.gz
gunzip common_all_20180418.vcf.gz
mv common_all_20180418.vcf human_snps.vcf
```

Execute SNP calling:

```
bash snp_calling --sqlitein SQLite/results.db --sqliteout SQLite/results.db --removeduplicates y -e ENS_hsa_92.db --picardpath /home/steven/tools/picard-tools-1.119/ --snpdbselected y --snpdb human_snp.vcf --toolsdir snp_tools/ --reads STAR/fastq1/untreat.sam
```

#### Assembly

Assembly of ORF products:

```
perl assembly.pl --sqliteRES SQLite/results.db --snp samtools_dbSNP --tis_ids 1
```

### TIS ID 2: PRICE

Run the PRICE wrapper:

```
python PRICE.py --result_db SQLite/results.db –m 16 
```

### TIS ID 3: SPECtre

Run the SPECtre wrapper:

```
source activate spectre
python SPECtre.py –-result_db SQLite/results.db -–offsets 29:12,30:12,31:12 –-untr_bam STAR/fastq1/untreat.bam --cores 60 –-threads_per_chrom 3
source deactivate
source activate proteoformer
```

SPECtre can only work with a certain amount of offsets, so we used the offsets of the main RPF lengths.

### FASTA file generation

Generate the redundant FASTA files and combine these:

```
perl generate_translation_db.pl --result_db SQLite/results.db --tis_ids 1 --mflag 5 --translation_db human_TIS_1_samtools_dbSNP_transcripts_redundant.fasta --var_file human_TIS_1_samtools_dbSNP_transcripts_VAR_redundant.txt
perl generate_translation_db.pl --result_db SQLite/results.db --tis_ids 2 --mflag 5 --translation_db human_TIS_2_transcripts_redundant.fasta --var_file human_TIS_2_transcripts_VAR_redundant.txt
perl generate_translation_db.pl --result_db SQLite/results.db --tis_ids 3 --mflag 5 --translation_db human_TIS_3_transcripts_redundant.fasta --var_file human_TIS_3_transcripts_VAR_redundant.txt
python combine_dbs.py --in_files human_TIS_1_samtools_dbSNP_transcripts_redundant.fasta,human_TIS_2_transcripts_redundant.fasta,human_TIS_3_transcripts_redundant.fasta --output_file jurkat_uniq_combined_output_redundant.fa --overview_file fasta_combination_overview_redundant.txt --venn_diagram venn_diagram_redundant.png
```

Combined FASTA files can be converted to SQLite databases for easy querying:

```
python combfasta2sqlite.py --fasta jurkat_uniq_combined_output_redundant.fa
```

### Combination with UniProt

Human UniProt database were downloaded from the [UniProt site](https://www.uniprot.org/) in FASTA format.
The combination was performed with the splice-included version of UniProt.

For the splice-included version:

```
python combine_with_uniprot.py -f ../compare_methods/jurkat_uniq_combined_output_redundant.fa -u ../../../HCT116/validation_uniprot_spliced/uniprotHumanProteomeIsoform.fasta
```

### MS validation with MaxQuant

For this validation we used MaxQuant version 1.6.1.0. MaxQuant is used as graphical user interface.
One of the different combined FASTA files was used as search space. Raw data can be found on the PRIDE repository under identifier PXD011353.
Most of the parameters were chosen as default, except the following:
* Methionine oxidation as fixed modification
* N-term acetylation as variable modification
* Match between runs
* MaxLFQ: True
* Only unique peptides for protein quantification, no razor peptides
* iBAQ: True

For each run, from the `combined/txt` MaxQuant folder, we used the following files:
* proteinGroups.txt
* peptides.txt
* msms.txt

### MS results parsing

For the combination files of the different proteoform calling techniques with redundancy:

```
python parse_maxquant.py -g proteinGroups.txt -p peptides.txt -x msms.txt -t ../fasta_combination_overview_redundant.txt -f ../jurkat_uniq_combined_output_redundant.fa -e ../ENS_hsa_92.db -r N -l protein
python parse_maxquant.py -g proteinGroups.txt -p peptides.txt -x msms.txt -t ../fasta_combination_overview_redundant.txt -f ../jurkat_uniq_combined_output_redundant.fa -e ../ENS_hsa_92.db -r N -l peptide
python parse_maxquant.py -g proteinGroups.txt -p peptides.txt -x msms.txt -t ../fasta_combination_overview_redundant.txt -f ../jurkat_uniq_combined_output_redundant.fa -e ../ENS_hsa_92.db -r N -l PSM
```

For the combination with UniProt spliced:

```
python parse_maxquant_uniprot.py -g proteinGroups.txt -p peptides.txt -x msms.txt -f ../comb_fasta_uniprotspliced.fasta -e ../../compare_methods/ENS_hsa_92.db  -l protein
python parse_maxquant_uniprot.py -g proteinGroups.txt -p peptides.txt -x msms.txt -f ../comb_fasta_uniprotspliced.fasta -e ../../compare_methods/ENS_hsa_92.db  -l peptide
python parse_maxquant_uniprot.py -g proteinGroups.txt -p peptides.txt -x msms.txt -f ../comb_fasta_uniprotspliced.fasta -e ../../compare_methods/ENS_hsa_92.db  -l PSM
python analyse_proteoforms.py -m max_proteins_uniprot.db -e ../../compare_methods/ENS_hsa_92.db -M ../../../../HCT116/validation_uniprot_canonical/proteoformeruniprot_validation1/HUMAN_9606_idmapping_selected.tab -f ../comb_fasta_uniprotspliced.fasta -g proteinGroups.txt -p peptides.txt -c ../../../../../../tools/clustalo/clustalo
```