Proteoformer
============

A proteogenomic pipeline that delineates true *in vivo* proteoforms and generates a protein sequence
 search space for peptide to MS/MS matching.

## Table of contents
1. [Introduction](#introduction)
2. [Dependencies](#dependencies)
3. [Prepations](#preparations)
    1. [iGenomes reference information download](#igenomes)
    2. [Ensembl download](#ensembl)
    3. [UTR simulation for Prokaryotes](#prokaryotutr)
    4. [SRA parallel download](#sra)
4. [Main pipeline](#main)
    1. [General quality check: fastQC](#fastqc1)
    2. [Mapping](#mapping)
    3. [General quality check: fastQC](#fastqc2)
    4. [Specific ribosome profiling quality check: mappingQC](#mappingqc)
    5. [Transcript calling](#transcriptcalling)
        1. [Rule-based transcript calling](#rulebased)
        2. [RiboZINB](#ribozinb_trcal)
    6. [ORF calling](#orf_calling)
        1. [PROTEOFORMER classic ORF calling](#classic_orf_calling)
            1. [TIS calling](#tis_calling)
            2. [SNP calling](#snp_calling)
            3. [ORF assembly](#assembly)
        2. [PRICE](#price)
        3. [SPECtre](#spectre_tool)
        4. [Analysis ID overview table](#tis_overview)
    7. [Fasta file generation](#fasta_generation)
5. [Optional steps](#optional)
    1. ORF-based counts
    2. FLOSS calculation
    3. Feature summarization
6. MS validation
    1. SearchGUI and PeptideShaker
    2. MaxQuant
7. [Copyright](#copyright)
8. [More information](#moreinformation)


## Introduction <a name="introduction"></a>

PROTEOFORMER is a proteogenomic pipeline that delineates true *in vivo* proteoforms and generates a protein sequence
 search space for peptide to MS/MS matching. It can be combined with canonical protein databases or used independently
 for identification of novel translation products. The pipeline makes use of the recently developed next generation 
 sequencing strategy termed ribosome profiling (RIBO-seq) that provides genome-wide information on protein synthesis
 *in vivo*. RIBO-seq is based on the deep sequencing of ribosome protected mRNA fragments. RIBO-seq allows for the mapping
 of the location of translating ribosomes on mRNA with sub codon precision, it can indicate which portion of the genome 
 is actually being translated at the time of the experiment as well as account for sequence variations such as single 
 nucleotide polymorphism and RNA splicing.

The pipeline 
* aligns your ribosome profiling data to a reference genome
* checks the quality and general features of this alignments
* searches for translated transcripts
* searches for all possible proteoforms in these transcripts
* constructs counts for different feature levels and calculates FLOSS scores
* constructs fasta files which allow mass spectrometry validation

Most modules of this pipeline are provided with a built-in help message. Execute the script of choice with the `-h` or 
`--help` to get the full help message printed in the command line.

PROTEOFORMER is also available in galaxy: [http://galaxy.ugent.be](http://galaxy.ugent.be)

## Dependencies <a name="dependencies"></a>

Proteoformer is built in Perl 5, Python 2.7 and Bash. All necessary scripts are included in this GitHub repository.

To prevent problems with missing dependencies, we included all necessary dependencies in a [Conda](https://conda.io/docs/) environment.
For more information about Conda installation, click [here](https://conda.io/docs/user-guide/install/index.html).

Once conda is installed, make sure to have the right channel order by executing following commands in the same order as listed here:

```
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```

Then you can install all dependencies based on the yml file in the `dependency_envs` folder of this GitHub repository with following command:

```
conda env create -f Dependency_envs/proteoformer.yml
```

This installs a new Conda environment in which all needed Conda dependencies are installed and available, including Perl
and Python. If not mentioned otherwise, all tools of the PROTEOFORMER pipeline should be executed in this environment.
 To activate this new Conda environment:

```
source activate proteoformer
```

Some Perl packages are not included in Conda, so after installation and first activation of the new environment,
 execute following script:

```
perl install_add_perl_tools.pl
```

If you want to exit the proteoformer Conda environment:

```
source deactivate
```

### Additional environments for RiboZINB, SPECtre and SRA download <a name="add_envs"></a>

For some tools, we needed to construct separate environments with different versions of the underlying tools. For all 
the other tools, the proteoformer environment is used.

#### RiboZINB

```
conda env create -f Dependency_envs/ribozinb.yml
source activate ribozinb
```

#### SPECtre

```
conda env create -f Dependency_envs/spectre.yml
source activate spectre
```

#### SRA download

```
conda env create -f Dependency_envs/download_sra_parallel.yml
source activate download_sra_parallel
```

## Preparations <a name="preparations"></a>

### iGenomes reference information download <a name="igenomes"></a>

Mapping is done based on reference information in the form of iGenomes directories. These directories can easily
downloaded and constructed with the `get_igenomes.py` script in the `Additional_tools` folder. For example:

```
python get_igenomes.py -v 92 -s human -d /path/to/dir -r -c 15
```

Input arguments:

| Argument       | Default   | Description                                                                                                                                                              |
|----------------|-----------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| -d / --dir     | Mandatory | Directory wherein the igenomes tree structure will be installed                                                                                                          |
| -v / --version | Mandatory | Ensembl annotation version to download (Ensembl plant (for arabidopsis) has seperate annotation versions!)                                                               |
| -s / --species | Mandatory | Specify the desired species for which gene annotation files should be downloaded                                                                                         |
| -r / --remove  |           | If any, overwrite the existing igenomes structure for that species                                                                                                       |
| -c / --cores   | Mandatory | The amount of cores that will be used for downloading chromosomes files (Do not use more than 15 cores as the download server can only establish 15 connections at once) |
| -h / --help    |           | This useful help message                                                                                                                                                 |

The tool currently supports following species:

| Species                                                             | Input value species argument |
|---------------------------------------------------------------------|------------------------------|
| Homo sapiens                                                        | human                        |
| Mus musculus                                                        | mouse                        |
| Rattus norvegicus                                                   | rat                          |
| Drosophila melanogaster                                             | fruitfly                     |
| Saccharomyces cerevisiae                                            | yeast                        |
| Danio rerio                                                         | zebrafish                    |
| Arabidopsis thaliana                                                | arabidopsis                  |
| Caenorhabditis elegans                                              | c.elegans                    |
| Salmonella enterica subsp. enterica serovar Typhimurium str. SL1344 | SL1344                       |
| Mycobacterium abscessus atcc 19977                                  | MYC_ABS_ATCC_19977           |


### Ensembl download <a name="ensembl"></a>

After mapping, mostly reference annotation is used from Ensembl (exons, splicing, canonical translation initiation,...). 
This information is available as an SQLite database and is downloadable by using the `ENS_db.py` script of the `Additional_tools`
 folder. For example:
 
```
python ENS_db.py -v 92 -s human
```

Input arguments:

| Argument       | Default   | Description                                                                      |
|----------------|-----------|----------------------------------------------------------------------------------|
| -v / --version | Mandatory | Ensembl annotation version to download (supported versions: from 74)             |
| -s / --species | Mandatory | Specify the desired species for which gene annotation files should be downloaded |
| -h / --help    |           | Print this useful help message                                                   |

Currently supported species:

| Species                                                             | Input value species argument |
|---------------------------------------------------------------------|------------------------------|
| Homo sapiens                                                        | human                        |
| Mus musculus                                                        | mouse                        |
| Drosophila melanogaster                                             | fruitfly                     |
| Saccharomyces cerevisiae                                            | yeast                        |
| Caenorhabditis elegans                                              | c.elegans                    |
| Salmonella enterica subsp. enterica serovar Typhimurium str. SL1344 | SL1344                       |
| Mycobacterium abscessus atcc 19977                                  | MYC_ABS_ATCC_19977           |

### UTR simulation for Prokaryotes <a name="prokaryotutr"></a>

For Prokaryotes, no untranslated upstream regions (UTRs) exist. Although, offset callers, used during mapping, need these
regions in order to calculate P-site offsets. Therefore, for Prokaryotes, these UTRs need to be simulated with the 
`simulate_utr_for_prokaryotes.py` script in the `Additional_tools` folder. For example:

```
python simulate_utr_for_prokaryotes.py igenomes/SL1344/Ensembl/ASM21085v2/Annotation/Genes/genes_32.gtf > igenomes/SL1344/Ensembl/ASM21085v2/Annotation/Genes/genes_32_with_utr.gtf
# Move and copy GTFs
mv igenomes/SL1344/Ensembl/ASM21085v2/Annotation/Genes/genes_32.gtf igenomes/SL1344/Ensembl/ASM21085v2/Annotation/Genes/genes_32_without_utr.gtf
cp igenomes/SL1344/Ensembl/ASM21085v2/Annotation/Genes/genes_32_with_utr.gtf igenomes/SL1344/Ensembl/ASM21085v2/Annotation/Genes/genes_32.gtf
```

This outputs a new GTF file. Best to rename the old GTF file and copy the new one under the name of the original GTF file 
as shown in the example.

Additional documentation can be found in the help message of the module.

### SRA parallel download <a name="sra"></a>

If you download the raw data (FASTQ) from SRA, you can use the [SRA toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=std).
However, we made a module to speed up this downloading process by using multiple cores of your system for multiple files
at once. Use the specific environment for this tool, if using this module. For example:

```
source activate download_sra_parallel
./download_sra_parallel.sh -c 20 -f 3034567 -l 3034572 #This downloads all fastq data  from SRR3034572 up until SRR3034572 on 20 cores
```

## Main pipeline <a name="main"></a>

### General quality check: fastQC <a name="fastqc1"></a>

Before starting off the analysis, we believe it is useful to check the general features and quality of the raw data.
This can be done by running [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). FastQC is included in 
the proteoformer conda environment, so you do not need to download it separately.

For example, to run it on 20 cores:

```
fastqc yourdata.fastq -t 20
```

This generates an HTML output report with figures for assessing the quality and general features and a ZIP file 
(for exporting the results to another system). More information about the tool can be found on the [help page](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)
 of FastQC. A tutorial video of what to expect of the figures can be found [here](https://www.youtube.com/watch?v=bz93ReOv87Y).

### Mapping <a name="mapping"></a>

The first big task of the pipeline is mapping the raw data on the reference genome. All reference data should be 
downloaded as an [iGenomes folder](#igenomes).

First, some prefiltering of bad and low-quality reads is done by using the [FastX toolkit](http://hannonlab.cshl.edu/fastx_toolkit/).
 Also, the adapters are eventually clipped off, using the [FastQ Clipper](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastx_clipper_usage)
 (recommended) or the clipper included in [STAR](https://github.com/alexdobin/STAR) (faster).
 
Mapping can be done by using [STAR](https://github.com/alexdobin/STAR), [TopHat](https://ccb.jhu.edu/software/tophat/index.shtml)
 or [BowTie](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml). BowTie is less preferable as this not
includes splicing. Before mapping against the genome, it is possible to filter out contaminants of PhiX, rRNA, sn(o)RNA 
and tRNA with the same mapping tool you chose to map against the genome.

After mapping, SAM and BAM files with aligned data are obtained. [Plastid](https://plastid.readthedocs.io/en/latest/examples/p_site.html)
 can be used to determine the P site offsets per RPF length. These offsets allow to pinpoint all reads to an exact base position 
 as explained [here](https://plastid.readthedocs.io/en/latest/examples/p_site.html). These offsets are very important 
 further down the pipeline to assign reads to the correct base position.
 
Next, these offsets are used to parse the alignment files into count tables. These count tables will be placed inside a 
 results SQLite database in which also the mapping statistics and the arguments will be put. During all following steps 
 of the pipeline, all results will be stored in this database and are available for consultation. We recommend to use the
 [sqlite3](https://www.sqlite.org/cli.html) command line shell for easy consultation of this database. More information 
 can be found on their [website](https://www.sqlite.org/cli.html). Use the following command for starting up sqlite3:
 
```sqlite3 SQLite/results.db```

For visualization, BedGraph files are generated. These can be used on different genome browsers like [UCSC](https://genome.ucsc.edu).

The mapping can be run as in one of following examples:

```
#Combination of untreated and treated data
perl mapping.pl --inputfile1 link/to/your/untr_data.fq --inputfile2 link/to/your/tr_data.fq --name my_experiment --species human --ensembl 92 --cores 20 --unique Y --igenomes_root /user/igenomes --readtype ribo

#Only untreated data with fastx clipping
perl mapping.pl --inputfile1 link/to/your/untr_data.fq --name my_experiment --species human --ensembl 92 --cores 20 --unique Y --igenomes_root /user/igenomes --readtype ribo_untr --clipper fastx --adaptor AGATCGGAAGAGCACAC

#To automatically determine P site offsets with Plastid and parse the results into count tables
perl mapping.pl --inputfile1 link/to/your/untr_data.fq --inputfile2 link/to/your/tr_data.fq --name my_experiment --species human --ensembl 92 --cores 20 --unique Y --igenomes_root /user/igenomes --readtype ribo --suite plastid

#With extra prefiltering, clipping, P site offset calling, count table parsing, extra file generation
perl mapping.pl --inputfile1 link/to/your/untr_data.fq --inputfile2 link/to/your/tr_data.fq --name my_experiment --species human --ensembl 92 --cores 20 --unique Y --igenomes_root /user/igenomes --readtype ribo --clipper fastx --adaptor CTGTAGGCACCATCAAT --phix Y --rRNA Y --snRNA Y --tRNA Y --rpf_split Y --pricefiles Y --suite plastid
```

Input arguments:

| Argument            | Default                                           | Description                                                                                                                                                                                                                                                                                                                                                      |
|---------------------|---------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| --inputfile1        | Mandatory                                         | The fastq file of the 'untreated' data for RIBO-seq (no drug,CHX,EMT) or the 1st fastq for single/paired-end RNA-seq                                                                                                                                                                                                                                                  |
| --inputfile2        | Mandatory for ribo and PE readtypes               | The fastq file of the treated data for RIBO-seq (PUR,LTM,HARR) or the 2nd fastq for paired-end RNA-seq,mandatory argument                                                                                                                                                                                                                                        |
| --name              | Mandatory                                         | Name of the run                                                                                                                                                                                                                                                                                                                                                  |
| --species           | Mandatory                                         | Species, possibilities: Mus musculus (mouse), Rattus norvegicus (rat), Homo sapiens (human), Drosophila melanogaster (fruitfly), Arabidopsis thaliana (arabidopsis), Salmonella (SL1344), Danio rerio (zebrafish), Saccharomyces cerevisiae (yeast)                                                                                                              |
| --ensembl           | Mandatory                                         | Ensembl annotation version                                                                                                                                                                                                                                                                                                                                       |
| --cores             | Mandatory                                         | Available cores for the pipeline to run on                                                                                                                                                                                                                                                                                                                       |
| --unique            | Mandatory                                         | Y/N. Retain only unique mapped reads (Y) or also multimapped reads (N)                                                                                                                                                                                                                                                                                           |
| --igenomes_root     | Mandatory                                         | Root folder of igenomes folder structure                                                                                                                                                                                                                                                                                                                         |
| --readtype          | ribo                                              | The readtype, possbilities: RIBOseq untreated and treated sample (ribo), RIBOseq only untreated sample (ribo_untr), RNAseq Paired end PolyA enriched (PE_polyA), RNAseq Single end PolyA enriched (SE_polyA), RNAseq Paired end total (PE_total) or RNAseq Single end total (SE_total)                                                                           |
| --adaptor           | CTGTAGGCACCATCAAT                                 | The adaptor sequence that needs to be clipped of by the clipper                                                                                                                                                                                                                                                                                                  |
| --work_dir          | $CWD env setting                                  | Working directory                                                                                                                                                                                                                                                                                                                                                |
| --tmp               | work_dir/tmp                                      | Temporary files folder                                                                                                                                                                                                                                                                                                                                           |
| --min_l_plastid     | 22                                                | Minimum RPF length for Plastid                                                                                                                                                                                                                                                                                                                                   |
| --max_l_plastid     | 34                                                | Maximum RPF length for Plastid                                                                                                                                                                                                                                                                                                                                   |
| --offset_img_untr   | work_dir/plastid/run_name_untreated_p_offsets.png | Path to save the offset image of plastid in for untreated sample                                                                                                                                                                                                                                                                                                 |
| --offset_img_tr     | work_dir/plastid/run_name_treated_p_offsets.png   | Path to save the offset image of plastid in for treated sample                                                                                                                                                                                                                                                                                                   |
| --min_l_parsing     | 26 (25 for fruitfly)                              | Minimum length for count table parsing                                                                                                                                                                                                                                                                                                                           |
| --max_l_parsing     | 34                                                | Maximum length for count table parsing                                                                                                                                                                                                                                                                                                                           |
| --out_bg_s_untr     | untreat_sense.bedgraph                            | Output visual file for sense untreated count data                                                                                                                                                                                                                                                                                                                |
| --out_bg_as_untr    | untreat_antisense.bedgraph                        | Output visual file for antisense untreated count data                                                                                                                                                                                                                                                                                                            |
| --out_bg_s_tr       | treat_sense.bedgraph                              | Output visual file for sense treated count data                                                                                                                                                                                                                                                                                                                  |
| --out_bg_as_tr      | treat_antisense.bedgraph                          | Output visual file for antisense treated count data                                                                                                                                                                                                                                                                                                              |
| --out_sam_untr      | untreat.sam                                       | Output alignment file of untreated data                                                                                                                                                                                                                                                                                                                          |
| --out_sam_tr        | treat.sam                                         | Output alignment file of treated data                                                                                                                                                                                                                                                                                                                            |
| --out_bam_untr      | untreat.bam                                       | Binary output alignment file of untreated data                                                                                                                                                                                                                                                                                                                   |
| --out_bam_tr        | treat.bam                                         | Binary output alignment file of treated data                                                                                                                                                                                                                                                                                                                     |
| --out_bam_tr_untr   | untreat_tr.bam                                    | Binary output alignment file of untreated data in transcript coordinates (if tr_coord option selected)                                                                                                                                                                                                                                                           |
| --out_bam_tr_tr     | treat_tr.bam                                      | Binary output alignment file of treated data in transcript coordinates (if tr_coord option selected)                                                                                                                                                                                                                                                             |
| --out_sqlite        | work_dir/SQLite/results.db                        | Output results SQLite database                                                                                                                                                                                                                                                                                                                                   |
| --clipper           | none                                              | If and which clipper has to be used, possibilities: none, STAR, fastx                                                                                                                                                                                                                                                                                            |
| --phix              | N                                                 | Map to Phix data as prefilter before genomic mapping: Yes (Y) or No (N)                                                                                                                                                                                                                                                                                          |
| --rRNA              | Y                                                 | Map to rRNA data as prefilter before genomic mapping: Yes (Y) or No(N)                                                                                                                                                                                                                                                                                           |
| --snRNA             | N                                                 | Map to sn(-o-)RNA data as prefilter before genomic mapping: Yes (Y) or No (N)                                                                                                                                                                                                                                                                                    |
| --tRNA              | N                                                 | Map to tRNA data as prefilter before genomic mapping: Yes (Y) or No (N)                                                                                                                                                                                                                                                                                          |
| --tr_coord          | N                                                 | Generate alignment files based on transcript coordinates: Yes (Y) or No (N)                                                                                                                                                                                                                                                                                      |
| --truseq            | Y                                                 | Whether strands (+ and -) are assigned as in TruSeq for RNAseq: Yes (Y) or No (N)                                                                                                                                                                                                                                                                                |
| --mismatch          | 2                                                 | Alignment will be outputted only if it has fewer mismatches than this value                                                                                                                                                                                                                                                                                      |
| --maxmultimap       | 16                                                | Alignments will be outputted only if the read maps fewer times than this value                                                                                                                                                                                                                                                                                   |
| --splicing          | Y                                                 | Allow splicing for genome alignment: Yes (Y) or No (N)                                                                                                                                                                                                                                                                                                           |
| --firstrankmultimap | N                                                 | Only retain the first ranked alignment when multimapping is allowed: Yes (Y) or No (N)                                                                                                                                                                                                                                                                           |
| --rpf_split         | N                                                 | Whether the program generates RPF specific bedgraph files: Yes (Y) or No (N)                                                                                                                                                                                                                                                                                     |
| --pricefiles        | N                                                 | If the program needs to generate SAM files specifically for PRICE: Yes (Y) or No (N), choose Y if you plan to execute PRICE later on                                                                                                                                                                                                                             |
| --suite             | custom                                            | Options of how to run the different mapping modules (mapping, mapping_plastid, mapping_parsing) all together for ribo data: only execute mapping and startup plastid and parsing manually later (custom), mapping and parsing with standard offsets (standard), mapping with afterwards plastid offset calculation and parsing based on these offsets (plastid)  |
| --suite_tools_loc   | work_dir                                          | The foder with scripts of the subsequent mapping tools when using plastid or standard suite                                                                                                                                                                                                                                                                      |
| --help              |                                                   | Generate help message                                                                                                                                                                                                                                                                                                                                            |

If you chose the custom suite, you can execute the Plastid and the parsing part of the mapping yourself after `mapping.pl`.
 This gives you the opportunity to test different offset sets. As you can see, also a tab-separated offset list can be used
 in the parsing by inputting them as a txt file. In that way, you can let the program use your offsets of choice.
 
Examples of how to run these separated programs:

```
#Plastid for untreated sample
perl mapping_plastid.pl --out_sqlite SQLite/results.db --treated untreated
#Plastid for treated sample
perl mapping_plastid.pl --out_sqlite SQLite/results.db --treated treated
#Parsing of all results into counts (mapping itself, Plastid untreated, Plastid treated)
perl mapping_parsing.pl --out_sqlite SQLite/results.db --offset plastid
```

Input arguments of mapping_plastid.pl:

| Argument     | Default                                             | Description                                                        |
|--------------|-----------------------------------------------------|--------------------------------------------------------------------|
| --work_dir   | CWD env setting                                     | The working directory                                              |
| --tmp        | work_dir/tmp                                        | Folder to store the temporary files in                             |
| --out_sqlite | work_dir/SQLite/results.db                          | SQLite results database                                            |
| --treated    | untreated                                           | Which sample to calculate offsets for (options: untreated/treated) |
| --offset_img | work_dir/plastid/run_name_(un)treated_p_offsets.png | P-site offsets output image path                                   |
| --help       |                                                     | Print help message                                                 |

Input arguments of mapping_parsing.pl:

| Argument           | Default                       | Description                                                                                                                                                                                            |
|--------------------|-------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| --work_dir         | CWD env setting               | The working directory                                                                                                                                                                                  |
| --tmp              | work_dir/tmp                  | Folder to store the temporary files in                                                                                                                                                                 |
| --out_sqlite       | work_dir/SQLite/results.db    | SQLite results database                                                                                                                                                                                |
| --offset           | standard                      | Offset option: you can use the standard offsets from Ingolia 2009 (standard), you can add offsets as a tab-separated txt file (from_file) or you can use the offsets determined with Plastid (plastid) |
| --offset_file_untr | mandatory if offset=from_file | Offset input tab-separated txt file for untreated sample                                                                                                                                               |
| --offset_file_tr   | mandatory if offset=from_file | Offset input tab-separated txt file for treated sample                                                                                                                                                 |
| --help             |                               | Print help message                                                                                                                                                                                     |

**Output:**

After running the full mapping step, different tables will be generated in the SQLite results database.

A table of all inputted arguments:

| Variable | Value       |
|----------|-------------|
| run_name | experiment1 |
| ...      | ...         |

A table with mapping statistics:

| Sample             | type    | total    | mapped_U | mapped_M | mapped_T | unmapped | map_freq_U        | map_freq_M        | map_freq_T      |
|--------------------|---------|----------|----------|----------|----------|----------|-------------------|-------------------|-----------------|
| experiment1.fastq1 | genomic | 31133026 | 17607183 | 1240756  | 18847939 | 12285087 | 0.565546792656775 | 0.039853369863554 | 0.3945998375403 |
| ...                | ...     | ...      | ...      | ...      | ...      | ...      | ...               | ...               | ...             |

A table with counts per base position for each sample:

| chr | strand | start  | count |
|-----|--------|--------|-------|
| 8   | -1     | 116252 | 2.0   |
| ... | ...    | ...    | ...   |

A table with counts per base position and per RPF length for each sample:

| chr | strand | start  | RPF | count |
|-----|--------|--------|-----|-------|
| 8   | -1     | 116252 | 28  | 2.0   |
| ... | ...    | ...    | ... | ...   |

Furthermore, BED and BedGraph files are generated, which allow visualizing these counts in a genome browser.

**Extra:** To obtain normalized BedGraph files based on the different library sizes, there is an extra tool for that in 
the `Additional_tools` folder. As input, it takes the different original BedGraph files and library sizes of both samples 
(i.e. the total mapped reads against the genomic reference, which can be found in the statistics table of the results database).

```
bash normBedgraph --untrs output/untreat_sense.bedgraph --untras output/untreat_antisense.bedgraph --nttrs output/treat_sense.bedgraph --nttras output/treat_antisense.bedgraph --libuntr 37873493 --libtr 45427218
```

### General quality check: fastQC <a name="fastqc2"></a>

This step is quite similar to the [first quality check step with FastQC](#fastqc1), although this time FastQC will be 
performed on the alignment files (sam). This generates again HTML reports but due to the adapter clipping and pre-filtering
during the mapping step, quality will normally have remarkably improved.

```
fastqc output/untreat.sam -t 20
fastqc output/treat.sam -t 20
```

### Specific ribosome profiling quality check: mappingQC <a name="mappingqc"></a>

MappingQC generates some figures which give a nice overview of the quality and the general feature of the aligned
 ribosome profiling data. More specific, it gives an overview of the P site offset calculation, the gene distribution
 and the metagenic classification. Furthermore, MappingQC does a thorough analysis of the triplet periodicity and the
 linked triplet phase (typical for ribosome profiling) in the canonical transcript of your data. Especially, the link 
 between the phase distribution and the RPF length, the relative sequence position and the triplet identity are taken 
 into account. MappingQC is also available as a stand-alone tool, separated from PROTEOFORMER and its SQLite results 
 structure, so that you can apply it on every samfile you like, independent of the mapping source.

For PROTEOFORMER, mappingQC needs a SAM file (from one of both samples), an [Ensembl database](#ensembl) and the results database to run. It generates an
HTML file which gives an overview of all generated figures. For exporting this report, a ZIP file is generated with the 
HTML and all figures included. The tool needs a tool directory with different background scripts. The directory (`mqc_tools`)
is available in our GitHub repository under `2_mappingQC`. Input the path of where you put this directory in the `--tool_dir`
argument.

Example:

```
perl mappingQC.pl --samfile output/untreat.sam --treated untreated --cores 20 --result_db SQLite/results.db --unique Y --ens_db ENS_hsa_92.db --offset plastid --plastid plastid/experiment1_untreated_p_offsets.png --tool_dir mqc_tools --plotrpftool pyplot3D
```

Input arguments:

| Argument        | Default                            | Description                                                                                                                                              |
|-----------------|------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------|
| --work_dir      | CWD env setting                    | The working directory                                                                                                                                    |
| --tmp           | work_dir/tmp                       | The temporary files folder                                                                                                                               |
| --samfile       | Mandatory                          | The SAM file to do the analysis for                                                                                                                      |
| --treated       | untreated                          | Whether the SAM file is from the treated or untreated sample (possiblities: untreated/treated)                                                           |
| --cores         | 5                                  | The amount of cores that can be used by the program                                                                                                      |
| --result_db     | Mandatory                          | The result database with all mapping results                                                                                                             |
| --unique        | Y                                  | Whether to use only the unique alignments (Y/N) (has to be Y if the mapping was done uniquely)                                                           |
| --ens_db        | Mandatory                          | The Ensembl database with annotation info                                                                                                                |
| --offset        | standard                           | The source of offsets: calculated with Plastid during mapping (plastid), standard offsets from Ingolia 2019 (standard) or from an input file (from_file) |
| --offset_file   | Mandatory if offset=from_file      | The offsets input file                                                                                                                                   |
| --offset_img    | Mandatory if offset=plastid        | The offsets image Plastid generated during mapping                                                                                                       |
| --output_folder | work_dir/mappingQC_output          | The output folder for storing output files                                                                                                               |
| --tool_dir      | work_dir/mqc_tools                 | The directory with all necessary underlying tools                                                                                                        |
| --plotrpftool   | grouped2D                          | Module used for plotting the RPF-phase figure: Seaborn grouped 2D chart (grouped2D), mplot3d 3D bar chart (pyplot3D) or mayavi 3D bar chart (mayavi)     |
| --html          | work_dir/mappingqc_out.html        | The output HTML file                                                                                                                                     |
| --zip           | work_dir/mappingQC_(un)treated.zip | The output ZIP file                                                                                                                                      |
| --help          |                                    | This helpful screen                                                                                                                                      |

**Caution:** The mplot3d package in Pyplot is vulnerable for so-called Escher effects. Sometimes, certain boxes are plotted
with a wrong z-order, but overall, figures are clear. On the other hand, the Mayavi package requires the usage of a graphical
card, so on servers, this is mostly not an option.

### Transcript calling <a name="transcriptcalling"></a>

After checking the aligned data for quality and general features, you can search for the translated transcripts.

#### Rule-based transcript calling <a name="rulebased"></a>

A first way to determine these translated transcript, is based on general rules. Transcript without RIBO-seq counts are 
ignored from the start. Then, for each exon of the transcript, the counts of ribosome reads are calculated and normalized 
by exon length. If the normalized count of an exon is 5 times lower than the mean normalized exon count of the given transcript,
the exon is classified as a noise exon. Transcripts with less than 15% noise exons, are called as actively-translated 
transcripts (`exon_coverage = Yes` in the output table).

An example of how to run this tool:

```
perl ribo_translation.pl --in_sqlite SQLite/results.db --out_sqlite SQLite/results.db --ens_db ENS_hsa_92.db
```

Input arguments:

| Argument     | Default                         | Description                                         |
|--------------|---------------------------------|-----------------------------------------------------|
| --work_dir   | CWD env setting                 | The working directory                               |
| --tmp        | work_dir/tmp                    | The temporary files folder                          |
| --in_sqlite  | SQLite/results.db               | The SQLite results database from previous steps     |
| --out_sqlite | The same as in_sqlite atrgument | The SQLite results database to store all results in |
| --ens_db     | Mandatory                       | The Ensembl database with annotation info           |
| --help       |                                 | Generate help message                               |

Output table structure in the SQLite results database:

| transcript_id | stable_id       | chr | seq_region_id | seq_region_strand | seq_region_start | seq_region_end | read_counts | normalized_counts | biotype        | exon_coverage | canonical | ccds      | gene_stable_id  | FPKM             |
|---------------|-----------------|-----|---------------|-------------------|------------------|----------------|-------------|-------------------|----------------|---------------|-----------|-----------|-----------------|------------------|
| 196519        | ENST00000371471 | 20  | 131538        | -1                | 53567065         | 53593839       | 585         | 0.187680461982676 | protein_coding | Yes           | Yes       | CCDS13443 | ENSG00000171940 | 10.6593122808274 |
| ...           | ...             | ... | ...           | ...               | ...              | ...            | ...         | ...               | ...            | ...           | ...       | ...       | ...             | ...              |

#### RiboZINB <a name="ribozinb_trcal"></a>

Another way to call transcripts, is by using the [RiboZINB](https://github.com/Biobix/RiboZINB) tool. This tool makes use 
of the zero-inflated binomial model to determine actively translated transcript isoforms. RiboZINB was directly included 
in the PROTEOFORMER pipeline. **The RiboZINB tool itself is still in beta and statistics are still in development and validation.**

The RiboZINB tool requires its own [earlier installed](#add_envs) Conda environment to run:

```
source activate ribozinb
```

An example of how to run this tool:

```
python RiboZINB.py -p SQLite/results.db
```

Input arguments:

| Argument           | Default                | Description                                                                                        |
|--------------------|------------------------|----------------------------------------------------------------------------------------------------|
| -w/--work_dir      | CWD env setting        | The working directory                                                                              |
| -x/--tmpfolder     | work_dir/tmp           | The temporary files folder                                                                         |
| -i/--in_sqlite     | Mandatory              | The input SQLite results database                                                                  |
| -o/--out_sqlite    | Same path as in_sqlite | The output SQLite results database                                                                 |
| -e/--ens_db        | Mandatory              | The Ensembl annotation database                                                                    |
| -m/--mincount      | 5                      | The minimum reads for a transcript to be called                                                    |
| -n/--no_of_samples | 30                     | The number of iterations when generating a negative set                                            |
| -f/--fdr           | 0.05                   | The false discovery rate                                                                           |
| -s/--default_score | d                      | Use the default score threshold (d) or estimate the threshold by performing a permutation test (p) |
| -v/--cutoff        | 0.1                    | The default score threshold                                                                        |
| -a/--alpha         | 1                      | Proportion of noise when generating the negative set                                               |

Output table structure in the SQLite results database:

| transcript_id | stable_id       | chr | seq_region_id | seq_region_strand | seq_region_start | seq_region_end | read_counts | normalized_counts | biotype        | exon_coverage | canonical | ccds     | gene_stable_id  | FPKM          |
|---------------|-----------------|-----|---------------|-------------------|------------------|----------------|-------------|-------------------|----------------|---------------|-----------|----------|-----------------|---------------|
| 173461        | ENST00000314167 | 4   | 131552        | -1                | 849278           | 932298         | 1513        | 0.340612336785    | protein_coding | NA            | Yes       | CCDS3340 | ENSG00000178950 | 19.3450784708 |
| ...           | ...             | ... | ...           | ...               | ...              | ...            | ...         | ...               | ...            | ...           | ...       | ...      | ...             | ...           |

### ORF calling <a name="orf_calling"></a>

Based on the translated transcripts, you can search for the translatiion product candidates. Three methods are currently 
implemented: [PROTEOFORMER classic ORF calling](#classic_orf_calling), [PRICE](#price) and [SPECtre](#spectre_tool). For an 
overview of all analyses that you performed, an [overview table](#tis_overview) for all the analysis ID's can be generated.

#### PROTEOFORMER classic ORF calling <a name="classic_orf_calling"></a>

A first method to determine the candidate translation products based on ribosome profiling is PROTEOFORMER classic method. 
This method is rule-based and therefore does not work based on a score. However, this method is developped with MS validation 
afterwards in mind. Therefore, defaultly it works quite unstringent, allowing to filter stringently later on during MS. As this 
method is also based on TIS calling, it needs a initiation ribosome profile to work (i.e. a second sample treated with 
LTM or HARR). The classic ORF calling can be separated in three subsequent steps: [TIS calling](#tis_calling), 
[SNP calling](#snp_calling) and [ORF assembly](#assembly). The SNP calling step is optional.

##### TIS calling <a name="tis_calling"></a>

The first step of this method determines the translation initiation sites based on the initiation profile (HARR/LTM).
Only (near-)cognate start codons are considered as a possible translation initiation site (TIS). This means that a codon 
can only differ in one base from the canonical 'ATG' start codon to be considered as a possible TIS. Furthermore, a local 
max in the initiation profile counts needs to be observable in the inputted `local_max` range. The TIS peak needs also to 
contain more counts than the `min_count` parameter for its annotation class to be called and the R<sub>LTM</sub>-R<sub>CHX</sub>
value needs to be above the threshold set for its annotation class.

An example of how to run this tool:

```
perl TIScalling_categorised --sqlite_db SQLite/results.db
```

Input arguments:

| Argument           | Default         | Description                                                                                   |
|--------------------|-----------------|-----------------------------------------------------------------------------------------------|
| --dir              | CWD env setting | The working directory                                                                         |
| --tmp              | work_dir/tmp    | The temporary files folder                                                                    |
| --sqlite_db        | Mandatory       | The SQLite results database with all mapping and translated transcript calling results        |
| --local_max        | 1               | The range wherein the local maximum can be located (e.g. 1 means +/- one triplet)             |
| --min_count_aTIS   | 5               | The minimum count of ribosome profiling reads that need to map to the aTIS                    |
| --R_aTIS           | 0.01            | The Rltm - Rchx value threshold for aTIS ORFs                                                 |
| --min_count_5      | 10              | The minimum count of ribosome profiling reads that need to map to a 5'UTR ORF TIS             |
| --R_5              | 0.05            | The Rltm - Rchx value threshold for 5'UTR ORFs                                                |
| --min_count_CDS    | 15              | The minimum count of ribosome profiling reads that need to map to a CDS ORF TIS               |
| --R_CDS            | 0.15            | The Rltm - Rchx value threshold for CDS ORFs                                                  |
| --min_count_3      | 10              | The minimum count of ribosome profiling reads that need to map to a 3'UTR ORF TIS             |
| --R_3              | 0.05            | The Rltm - Rchx value threshold for 3'UTR ORFs                                                |
| --min_count_ntr    | 10              | The minimum count of ribosome profiling reads that need to map to a non-codign region ORF TIS |
| --R_ntr            | 0.05            | The Rltm - Rchx value threshold for non-coding region ORFs                                    |
| --out_sqlite       | /               | Parameter is only used in Galaxy, not in CLI                                                  |
| --transcriptfilter | none            | Use certain filters at transcript level (options: none, canonical, ccds)                      |

Called TIS's are saved in the SQLite results database in following table format:

| transcript_id | stable_id       | biotype        | chr | strand | start    | dist_to_transcript_start | dist_to_aTIS | annotation | exon | aTIS_call | start_codon | peak_shift | count | Rltm_min_Rchx     |
|---------------|-----------------|----------------|-----|--------|----------|--------------------------|--------------|------------|------|-----------|-------------|------------|-------|-------------------|
| 669998        | ENST00000337132 | protein_coding | 1   | 1      | 16440762 | 91                       | 0            | aTIS       | 1    | True      | ATG         | +1 0 -1    | 130.0 | 0.961529860001382 |
| ...           | ...             | ...            | ... | ...    | ...      | ...                      | ...          | ...        | ...  | ...       | ...         | ...        | ...   | ...               |

##### SNP calling <a name="snp_calling"></a>

Additional to the information about called TIS positions, it is also possible to search for SNP positions based on the 
ribosome profilng data.
More information about the SNP calling can be found in its separate [README file](https://github.com/Biobix/proteoformer/blob/master/4_ORF_calling/using_treated_TIS_calling/4b_variation_calling/SNPanalysis/READMEsnpcalling.txt).

With [Picard tools](https://broadinstitute.github.io/picard/)(for removing duplicates) and [SAM tools](http://www.htslib.org), SNP positions 
will be called based on the aligned ribosome profiling data. Data from [SNPdb](https://www.ncbi.nlm.nih.gov/snp) can 
also be included in the analysis.

An example of how to run this tool:

```
bash snp_calling --sqlitein path/to/results/database.db --sqliteout path/to/output/database.db --removeduplicates [y|n] --picardpath /path/to/picardmap --snpdbselected [y|n] --snpdb path/to/snpdb --toolsdir /path/to/tooldir --reads path/to/mapped/reads.sam --mincoverage 3 --maxcoverage 100 --high_af 0.95 --lower_af 0.3 --upper_af 0.7
```

Input arguments:

| Argument           | Default                             | Description                                                                                                                                                                              |
|--------------------|-------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| -s / --sqlitein    | Mandatory                           | Path to the database with results from previous steps                                                                                                                                    |
| --sqliteout        | Mandatory                           | Path to store the database with all results after this step                                                                                                                              |
| --removeduplicates | Mandatory                           | 'y' (yes) or 'n' (no), whether to remove duplicate reads with Picard                                                                                                                     |
| --picardpath       | Mandatory if previous option is 'y' | Path to the Picard tools jar file                                                                                                                                                        |
| --snpdbselected    | Mandatory                           | 'y' (yes) or 'n' (no), whether he mapped reads will be searched for known SNPs data                                                                                                      |
| --snpdb            | Mandatory if previous option is 'y' | Path to the SNPdb                                                                                                                                                                        |
| --toolsdir         | Mandatory                           | Path to the folder where filterSAMfile.pl, snpIndexBuilder.pl and splitVCFaltRecords.pl are                                                                                              |
| -r / --reads       | Mandatory                           | Path to the SAM file containing alignments ('untreated')                                                                                                                                 |
| --mincoverage      | 3                                   | SAMtools parameter, the minimal number of reads that need to map at a location so that a SNP can be called there                                                                         |
| --maxcoverage      | 100                                 | SAMtools parameter, the maximal number of reads that need to map at a location so that a SNP can be called there                                                                         |
| --high_af          | 0.95                                | High, lower and upper allelic frequency, input variables for mergeVCFfiles.pl, select SNPs & INDELS when their allelic frequency is between lower_af and upper_af or higher than high_af |
| --lower_af         | 0.3                                 | Lower allelic frequency                                                                                                                                                                  |
| --upper_af         | 0.7                                 | Upper allelic frequency                                                                                                                                                                  |

The reasoning behind the high, lower and upper allelic frequency is providing the ability to select both hetero- and
monozygotic SNP variations.

Called SNPs are stored in the SQLite results database in following table format:

| id           | chr | pos    | ref | alt | dp  | af  | new | seq_region_id |
|--------------|-----|--------|-----|-----|-----|-----|-----|---------------|
| 100001506993 | 10  | 150699 | A   | C   | 1   | 0.5 | m   | 131544        |
| ...          | ... | ...    | ... | ... | ... | ... | ... | ...           |

The column "new" is added when the results are compared to SNPdb. A "y" (yes) indicates that a variant was new, ie NOT
found in SNPdb. A "n" (no) means not new (ie found in SNPdb) and an "m" (mismatch) is used for those mismatches in the 
mapped reads that were found in SNPdb. There is no allelic frequency information for the "m" variants, so this has been
set to 0.0.

##### ORF assembly <a name="assembly"></a>

Information about called [TIS's](#tis_calling) and [SNPs](#snp_calling) can be used to construct the candidate translation products *in silico*. In this 
process, the program keeps track of the exonic structure of the transcript and known selenocysteines (encoded by an 'UGA'
codon, which can also function as a STOP signal). Both the translation product with the selenocysteine and the earlier 
terminated product will be kept. For near-cognate start sites, the first codon is replace for a cognate
methionine. Furthermore, translation products will only be outputted if the translation stops at a valid stop signal.
The coverage and FPKM value per translation product will be calculated as well.

An example of how to run this tool:

```
perl assembly.pl --sqliteRES SQLite/results.db --snp samtools_dbSNP --tis_ids 1
```

Input arguments:

| Argument         | Default                     | Description                                                                                   |
|------------------|-----------------------------|-----------------------------------------------------------------------------------------------|
| --dir            | CWD env setting             | The working directory                                                                         |
| --tmp            | work_dir/tmp                | The temporary files folder                                                                    |
| --sqliteRES      | Mandatory                   | The SQLite results database with all previous results                                         |
| --local_max      | 1                           | The range wherein the local maximum can be located (e.g. 1 means +/- one triplet)             |
| --min_count_aTIS | 5                           | The minimum count of ribosome profiling reads that need to map to the aTIS                    |
| --R_aTIS         | 0.05                        | The Rltm - Rchx value threshold for aTIS ORFs                                                 |
| --min_count_5    | 10                          | The minimum count of ribosome profiling reads that need to map to a 5'UTR ORF TIS             |
| --R_5            | 0.05                        | The Rltm - Rchx value threshold for 5'UTR ORFs                                                |
| --min_count_CDS  | 15                          | The minimum count of ribosome profiling reads that need to map to a CDS ORF TIS               |
| --R_CDS          | 0.15                        | The Rltm - Rchx value threshold for CDS ORFs                                                  |
| --min_count_3    | 10                          | The minimum count of ribosome profiling reads that need to map to a 3'UTR ORF TIS             |
| --R_3            | 0.05                        | The Rltm - Rchx value threshold for 3'UTR ORFs                                                |
| --min_count_ntr  | 10                          | The minimum count of ribosome profiling reads that need to map to a non-codign region ORF TIS |
| --R_ntr          | 0.05                        | The Rltm - Rchx value threshold for non-coding region ORFs                                    |
| --out_sqlite     | same as sqliteRES parameter | The SQLite database holding all output                                                        |
| --snp            | NO                          | The applied SNP algorithm. Options: 'NO', 'samtools', 'samtools_dbSNP'                        |
| --tis_ids        | Mandatory                   | List of the analysis ids (one integer or a comma-separated list of integers)                  |

Candidate translation products are stored in the SQLite results database in following table format:

| tr_stable_id    | chr | strand | start     | start_codon | stop      | starts_list                                                                                         | ends_list                                                                                           | dist_to_transcript_start | dist_to_aTIS | annotation | aTIS_call | peak_shift | count  | Rltm_min_Rchx | coverage          | FPKM             | SNP         | INDEL | secs | tr_seq                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                | aa_seq                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
|-----------------|-----|--------|-----------|-------------|-----------|-----------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------|--------------------------|--------------|------------|-----------|------------|--------|---------------|-------------------|------------------|-------------|-------|------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| ENST00000368236 | 1   | 1      | 156581365 | ATG         | 156586550 | 156581365_156582231_156583042_156583338_156583797_156584877_156585125_156585713_156585949_156586469 | 156582070_156582434_156583170_156583450_156583895_156584974_156585231_156585826_156586045_156586550 | 37                       | 0            | aTIS       | no_data   | NA         | 2164.0 | NA            | 0.183533447684391 | 9.89861772252107 | 692_A_G_0.5 |       |      | ATGTCTTCCCCCAACCCTGAGGATGTGCCCCGGAGGCCAGAACCTGAGCCCTCAAGCTCCAATAAGAAAAAGAAGAAAAGAAAGTGGCTGCGGCAAGAAGCCAGCATCCAAGCCCTCACCAGGGCTGGCCATGGGGCCCTTCAGGCTGGCCAGAACCATGAAGCCTTGAACAACTTCCAGAGGGCCTTCCTTCTGGCCTCCAAGGCCCCACAAACCAGGGATACCCCTGTGCTCCAGGCCTGCGCCTTCAACCTGGGGGCTGCCTATGTGGAGACTGGGGACCCAGCCAGAGGCCTTGAGCTACTCCTGCGAGCCCACCCTGAAGAGAAGGCACAGGGCAGGCGACACGGCGACCAATGTTTCAATGTGGCTTTGGCCTACCATGCCCTCGGCGAGCTGCCTCAAGCTTTGGCCTGGTACCACAGGGCCCTGGGCCACTACCAGCCACAGGGTGACCAGGGAGAAGCCTGGGCAAAAATGGGAGCCTGCTACCAGGCTCTGGGACAGCCTGAGCTAGCAGCCCACTGCCTGCAGGAAGCAAGCCAGGCCTATGCTCAAGAGAGACAGCTGCGGGCCGCAGCCCTGGCACTGGGGGCTGCGGCAGGATGTATGCTGAAGAGTGGGCGGCATCGGGTGGGGGAAGTTGTGCAGGTGCTGGAGAAAAGCCGGAGGCTTGCCGAGAGGAGCACTGGGAGGCGACTGCTGGGGCACCTCTATAACGATCTAGGCCTGGGCTACTCCCAGCTCCAGCTGTTCCCGCTGGCCGTGGAGGCCTTCCTGCAGGCCCTGCCCCTGTGCTGGGTGCCAGGAGAGCAGGCCACAGTGCTAAGAAACCTCGGGATGGCCCACAATGCCCTCGGCAACTATCAGGAAGCTCGGGAGTTTCACCAGAAGGCTGCTGACCTACACGGCTCTGTGGGGCAGCGGTGGGAGCAGGGCCGGAGCTTTGGCAGCCTGGCCTTTGCATTGAGCCAGCTGGGGGACCACAAGGCTGCCAGAGACAACTACCTGCATGCCCTGCAGGCTGCCCGGGACTCTGGGGACATGAAGGGACAGTGGCAGGCCTGTGAGGGTCTGGGGGCTGCTGCAGCCAGGCTGGGGCAGTATGACCAGGCCTTGAAGTACTATAAGGAAGCACTGGCCCAGTGTCAGAAGGAGCCAGATTCTGTGCGAGAACGGCTGGTGGCCAAGCTGGCAGACACCGTGAGGACGCGCTTGGCCCAGGTGGGGCTGGTCCAGACTCACACCCTGACTTCGGCTCCGGGAAGACTCCAGGCTCCAGGTGGGGCCAGCCAGGCGGAGGGGACCCCAGCAAAGGCAGGAAGCAGCACAGCAGGTGTCCAGCACAGATCTTCCAGTGGGTGGGAAGATGAAGAGTTTGAGGAGGGCCACCAGAAGAAAAAAGAGGAGAGGTCGGCAAACGTTCCGGTGAGGGCTGGGCCGGGAAGACCAGAGCTGTGTTTCCTTCCAGGCACAGTGAATCATTCGCACCATCTAGCTTCTAGTTGCCCCACGTTTACCAAGCACACGCCCTGCAGAGGGACAGTCCTCGGCAAAGCCTCCATCTATAGTCCAGGACCCAGGGCCCATCTTCCATTTGTAGGTCCAGGCCCTCCCAGAGCGGAGTACCCTAGCATCTTGGTACCCAATGGCCCTCAAGCCAATAGGTCATCCAGGTGGCCCAGGGAAAGCCTCAGCAGGAGCCGCCAGAGGAGACCCATGGAGTCGGGCATCTGCACTATTGTGTGA | MSSPNPEDVPRRPEPEPSSSNKKKKKRKWLRQEASIQALTRAGHGALQAGQNHEALNNFQRAFLLASKAPQTRDTPVLQACAFNLGAAYVETGDPARGLELLLRAHPEEKAQGRRHGDQCFNVALAYHALGELPQALAWYHRALGHYQPQGDQGEAWAKMGACYQALGQPELAAHCLQEASQAYAQERQLRAAALALGAAAGCMLKSGRHRVGEVVQVLEKSRRLAERSTGRRLLGHLYNDLGLGYSQLQLFPLAVEAFLQALPLCWVPGEQATVLRNLGMAHNALGNYQEAREFHQKAADLHGSVGQRWEQGRSFGSLAFALSQLGDHKAARDNYLHALQAARDSGDMKGQWQACEGLGAAAARLGQYDQALKYYKEALAQCQKEPDSVRERLVAKLADTVRTRLAQVGLVQTHTLTSAPGRLQAPGGASQAEGTPAKAGSSTAGVQHRSSSGWEDEEFEEGHQKKKEERSANVPVRAGPGRPELCFLPGTVNHSHHLASSCPTFTKHTPCRGTVLGKASIYSPGPRAHLPFVGPGPPRAEYPSILVPNGPQANRSSRWPRESLSRSRQRRPMESGICTIV* |
| ...             | ... | ...    | ...       | ...         | ...       | ...                                                                                                 | ...                                                                                                 | ...                      | ...          | ...        | ...       | ...        | ...    | ...           | ...               | ...              | ...         | ...   | ...  | ...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |

#### PRICE <a name="price"></a>

Another method uses the [PRICE](https://github.com/erhard-lab/gedi/wiki/Price) algorithm to pick up translation from 
ribosome profiling data. PRICE (Probabilistic inference of codon activities by an EM algorithm) is a method to identify
translated ORFs. Make sure you have the right Java distribution as stated on the 
[wiki](https://github.com/erhard-lab/gedi/wiki/Price) of PRICE.

PROTEOFORMER contains a wrapper that fully executes PRICE and uses its results in order to obtain a candidate translation 
product table similar to the assembly table. PRICE does not account for selenocysteines, so these cannot be included 
with this methodology. Execution of PRICE comprises genome reference preparation, the actual PRICE run and conversion of 
the output CIT file. Afterwards the output of PRICE will be combined with information from Ensembl in order to obtain all
needed features of the translation product table. All these steps are included in the PROTEOFORMER wrapper program. The 
wrapper let PRICE use both treatment samples if they are available but it is also possible to run the PRICE wrapper with 
only an untreated/CHX sample.

An example of how to run this wrapper:

```
python PRICE.py -r SQLite/results.db
```

Input arguments:

| Argument       | Default         | Description                                                                   |
|----------------|-----------------|-------------------------------------------------------------------------------|
| -w/--workdir   | CWD env setting | Working directory                                                             |
| -t/--tmp       | workdir/tmp     | Temporary files folder                                                        |
| -r/--result_db | Mandatory       | The SQLite database with results of mapping and translated transcript calling |
| -f/--fdr       | 0.1             | FDR value threshold used in PRICE                                             |
| -m/--max_ram   | PRICE default   | Maximum amount of gigabytes RAM available for running PRICE                   |
| -h/--help      | /               | Show help message                                                             |

The output table has the following format:

| tr_stable_id    | chr | strand | start     | start_codon | stop      | starts_list                                                           | ends_list                                                             | dist_to_transcript_start | dist_to_aTIS | annotation | biotype        | aTIS_call | peak_shift | count  | Rltm_min_Rchx | coverage       | FPKM          | SNP | INDEL | secs | tr_seq                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           | aa_seq                                                                                                                                                                                                                                                                                                                                                                                                                                           | price_id                  | price_annotation | price_pval |
|-----------------|-----|--------|-----------|-------------|-----------|-----------------------------------------------------------------------|-----------------------------------------------------------------------|--------------------------|--------------|------------|----------------|-----------|------------|--------|---------------|----------------|---------------|-----|-------|------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------|------------------|------------|
| ENST00000218516 | X   | -1     | 101407909 | GTG         | 101397809 | 101397809_101398370_101398785_101400666_101401632_101403811_101407710 | 101398099_101398567_101398946_101400757_101401809_101403985_101407909 | 17                       | -6           | 5UTR       | protein_coding | NA        | NA         | 9878.0 | NA            | 0.545524691358 | 60.9777121395 |     |       |      | GTGACAATGCAGCTGAGGAACCCAGAACTACATCTGGGCTGCGCGCTTGCGCTTCGCTTCCTGGCCCTCGTTTCCTGGGACATCCCTGGGGCTAGAGCACTGGACAATGGATTGGCAAGGACGCCTACCATGGGCTGGCTGCACTGGGAGCGCTTCATGTGCAACCTTGACTGCCAGGAAGAGCCAGATTCCTGCATCAGTGAGAAGCTCTTCATGGAGATGGCAGAGCTCATGGTCTCAGAAGGCTGGAAGGATGCAGGTTATGAGTACCTCTGCATTGATGACTGTTGGATGGCTCCCCAAAGAGATTCAGAAGGCAGACTTCAGGCAGACCCTCAGCGCTTTCCTCATGGGATTCGCCAGCTAGCTAATTATGTTCACAGCAAAGGACTGAAGCTAGGGATTTATGCAGATGTTGGAAATAAAACCTGCGCAGGCTTCCCTGGGAGTTTTGGATACTACGACATTGATGCCCAGACCTTTGCTGACTGGGGAGTAGATCTGCTAAAATTTGATGGTTGTTACTGTGACAGTTTGGAAAATTTGGCAGATGGTTATAAGCACATGTCCTTGGCCCTGAATAGGACTGGCAGAAGCATTGTGTACTCCTGTGAGTGGCCTCTTTATATGTGGCCCTTTCAAAAGCCCAATTATACAGAAATCCGACAGTACTGCAATCACTGGCGAAATTTTGCTGACATTGATGATTCCTGGAAAAGTATAAAGAGTATCTTGGACTGGACATCTTTTAACCAGGAGAGAATTGTTGATGTTGCTGGACCAGGGGGTTGGAATGACCCAGATATGTTAGTGATTGGCAACTTTGGCCTCAGCTGGAATCAGCAAGTAACTCAGATGGCCCTCTGGGCTATCATGGCTGCTCCTTTATTCATGTCTAATGACCTCCGACACATCAGCCCTCAAGCCAAAGCTCTCCTTCAGGATAAGGACGTAATTGCCATCAATCAGGACCCCTTGGGCAAGCAAGGGTACCAGCTTAGACAGGGAGACAACTTTGAAGTGTGGGAACGACCTCTCTCAGGCTTAGCCTGGGCTGTAGCTATGATAAACCGGCAGGAGATTGGTGGACCTCGCTCTTATACCATCGCAGTTGCTTCCCTGGGTAAAGGAGTGGCCTGTAATCCTGCCTGCTTCATCACACAGCTCCTCCCTGTGAAAAGGAAGCTAGGGTTCTATGAATGGACTTCAAGGTTAAGAAGTCACATAAATCCCACAGGCACTGTTTTGCTTCAGCTAGAAAATACAATGCAGATGTCATTAAAAGACTTACTTTAA | MTMQLRNPELHLGCALALRFLALVSWDIPGARALDNGLARTPTMGWLHWERFMCNLDCQEEPDSCISEKLFMEMAELMVSEGWKDAGYEYLCIDDCWMAPQRDSEGRLQADPQRFPHGIRQLANYVHSKGLKLGIYADVGNKTCAGFPGSFGYYDIDAQTFADWGVDLLKFDGCYCDSLENLADGYKHMSLALNRTGRSIVYSCEWPLYMWPFQKPNYTEIRQYCNHWRNFADIDDSWKSIKSILDWTSFNQERIVDVAGPGGWNDPDMLVIGNFGLSWNQQVTQMALWAIMAAPLFMSNDLRHISPQAKALLQDKDVIAINQDPLGKQGYQLRQGDNFEVWERPLSGLAWAVAMINRQEIGGPRSYTIAVASLGKGVACNPACFITQLLPVKRKLGFYEWTSRLRSHINPTGTVLLQLENTMQMSLKDLL* | ENST00000218516_Variant_0 | Variant          | 1.6931e-11 |
| ...             | ... | ...    | ...       | ...         | ...       | ...                                                                   | ...                                                                   | ...                      | ...          | ...        | ...            | ...       | ...        | ...    | ...           | ...            | ...           | ... | ...   | ...  | ...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              | ...                                                                                                                                                                                                                                                                                                                                                                                                                                              | ...                       | ...              | ...        |

#### SPECtre <a name="spectre_tool"></a>

A third method is based on [SPECtre](https://github.com/mills-lab/spectre). SPECtre is a tool for searching actively 
translated regions in ribosome profiling using a spectral coherence classifier. SPECtre runs solely based on an 
untreated/CHX sample, so no initiation profile is necessary. However, SPECtre starts from the reference information in the
GTF file and is therefore less suited to pick up unknown new ORFs.

SPECtre needs its own [earlier installed](#add_envs) Conda environment to run:

```
source activate spectre
```

In this environment, a PROTEOFORMER wrapper for SPECtre can run. The wrapper runs SPECtre in the background, multithreaded
over all chromosomes and merges all results afterwards. Next, these results are parsed and supplied with information from 
Ensembl to obtain all features needed for the PROTEOFORMER translation product table. Furthermore, SPECtre takes 
selenocysteines into account.

An example of how to run the SPECtre wrapper:

```
python SPECtre.py -r SQLite/results.db -o 28:12,29:12,30:12 -c 20 -x 1
```

Input arguments:

| Argument               | Default                   | Description                                                                                                       |
|------------------------|---------------------------|-------------------------------------------------------------------------------------------------------------------|
| -w/--workdir           | CWD env setting           | The working directory                                                                                             |
| -t/--tmp               | workdir/tmp               | The temporary files folder                                                                                        |
| -r/--result_db         | Mandatory                 | SQLite database with results from mapping and translated transcript calling                                       |
| -b/--untr_bam          | Get from results database | BAM file of the untreated/CHX sample                                                                              |
| -o/--offsets           | Get from results database | Custom list of offsets, but defaultly the program takes these from the results database                           |
| -c/--cores             | Get from results database | Defaultly taken from the results database, but you can suggest a different amount of cores especially for SPECtre |
| -x/--threads_per_chrom | 1                         | The number of threads used per chromosome to run SPECtre                                                          |
| -f/--fdr               | 0.05                      | FDR threshold used in SPECtre                                                                                     |
| -h/--help              | /                         | Show help message                                                                                                 |

**Caution:** SPECtre is currently not suited to work with large offset list like outputted by Plastid. Therefore, it is 
advisable to only use a list of the P-site offsets of the most abundant RPF lengths.

The output table has the following format:

| tr_stable_id    | chr | strand | start   | start_codon | stop    | starts_list                                             | ends_list                                               | dist_to_transcript_start | dist_to_aTIS | annotation | biotype        | aTIS_call | peak_shift | count   | Rltm_min_Rchx | coverage       | FPKM          | SNP | INDEL | secs  | tr_seq                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | aa_seq                                                                                                                                                                                                                                              | spectre_id |
|-----------------|-----|--------|---------|-------------|---------|---------------------------------------------------------|---------------------------------------------------------|--------------------------|--------------|------------|----------------|-----------|------------|---------|---------------|----------------|---------------|-----|-------|-------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------|
| ENST00000588919 | 19  | 1      | 1104125 | ATG         | 1106766 | 1104125_1105186_1105366_1105658_1106242_1106378_1106540 | 1104127_1105280_1105510_1105809_1106266_1106459_1106766 | 53                       | 0            | aTIS       | protein_coding | NA        | NA         | 27921.0 | NA            | 0.552812071331 | 277.962315247 |     |       | 97_96 | ATGTGCGCGTCCCGGGACGACTGGCGCTGTGCGCGCTCCATGCACGAGTTTTCCGCCAAGGACATCGACGGGCACATGGTTAACCTGGACAAGTACCGGGGCTTCGTGTGCATCGTCACCAACGTGGCCTCCCAGTGAGGCAAGACCGAAGTAAACTACACTCAGCTCGTCGACCTGCACGCCCGATACGCTGAGTGTGGTTTGCGGATCCTGGCCTTCCCGTGTAACCAGTTCGGGAAGCAGGAGCCAGGGAGTAACGAAGAGATCAAAGAGTTCGCCGCGGGCTACAACGTCAAATTCGATATGTTCAGCAAGATCTGCGTGAACGGGGACGACGCCCACCCGCTGTGGAAGTGGATGAAGATCCAACCCAAGGGCAAGGGCATCCTGGGAAATGCCATCAAGTGGAACTTCACCAAGTTTGGACACCGTCTCTCCACAGTTCCTCATCGACAAGAACGGCTGCGTGGTGAAGCGCTACGGACCCATGGAGGAGCCCCTGGTGATAGAGAAGGACCTGCCCCACTATTTCTAGCTCCACAAGTGTGTGGCCCCGCCCGAGCCCCTGCCCACGCCCTTGGAGCCTTCCACCGGCACTCATGACGGCCTGCCTGCAAACCTGCTGGTGGGGCAGACCCGAAAATCCAGCGTGCACCCCGCCGGAGGAAGGTCCCATGGCCTGCTGGGCTTGGCTCGGCGCCCCCACCCCTGGCTACCTTGTGGGAATAA | MCASRDDWRCARSMHEFSAKDIDGHMVNLDKYRGFVCIVTNVASQUGKTEVNYTQLVDLHARYAECGLRILAFPCNQFGKQEPGSNEEIKEFAAGYNVKFDMFSKICVNGDDAHPLWKWMKIQPKGKGILGNAIKWNFTKFGHRLSTVPHRQERLRGEALRTHGGAPGDREGPAPLFLAPQVCGPARAPAHALGAFHRHSURPACKPAGGADPKIQRAPRRRKVPWPAGLGSAPPPLATLWE* | 55220      |
| ...             | ... | ...    | ...     | ...         | ...     | ...                                                     | ...                                                     | ...                      | ...          | ...        | ...            | ...       | ...        | ...     | ...           | ...            | ...           | ... | ...   | ...   | ...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       | ...                                                                                                                                                                                                                                                 | ...        |

#### Analysis ID overview table <a name="tis_overview"></a>

Over the different analysis strategies ([translated transcript calling](#transcriptcalling) and [ORF calling](#orf_calling)),
the applied parameters are kept in a TIS overview table. For each analysis executed and stored in a certain SQLite results 
database, PROTEOFORMER keeps track of this analysis ID and its corresponding parameters in a separate table. At any time,
the content of this analysis ID table can be printed in a tab file using following program:

```
perl TIS_overview.pl --sqlite_db SQLite/results.db
```
Input arguments:

| Argument           | Default                      | Description                      |
|--------------------|------------------------------|----------------------------------|
| --work_dir         | CWD env setting              | Working directory                |
| --sqlite_db        | Mandatory                    | SQLite results database          |
| --out_tis_overview | work_dir/SQLite/overview.tis | Analysis ID overview output file |

The overview table has following format:

| ID  | local_max | min_count_aTIS | R_aTis | min_count_5UTR | R_5UTR | min_count_CDS | R_CDS | min_count_3UTR | R_3UTR | min_count_ntr | R_ntr | PRICE_FDR | SPECTRE_FDR | SNP            | indel | filter | tr_calling | TIS_calling |
|-----|-----------|----------------|--------|----------------|--------|---------------|-------|----------------|--------|---------------|-------|-----------|-------------|----------------|-------|--------|------------|-------------|
| 1   | 1         | 5              | 0.01   | 10             | 0.05   | 15            | 0.15  | 10             | 0.05   | 10            | 0.05  |           |             | samtools_dbSNP | NO    | none   | rule_based | Yes         |
| ... | ...       | ...            | ...    | ...            | ...    | ...           | ...   | ...            | ...    | ...           | ...   | ...       | ...         | ...            | ...   | ...    | ...        | ...         |

This table will be outputted as a tab separated text file.

###FASTA file generation <a name="fasta_generation"></a>

Translation products estimated in a particular analysis ID can be outputted in FASTA file format. This allows to submit 
them in a proteomics analysis as search space. During this output generation, additional processes can be executed. The 
candidate translation products can be mapped to reference information and redundancy on subsequence level in between 
translation products can be removed.

An example of how to run this program for analysis ID 1:

```
perl generate_translation_db.pl --result_db SQLite/results.db --tis_ids 1 --mflag 4 
```

Input arguments:

| Argument            | Default                                       | Description                                                                                                                                                                                                                                  |
|---------------------|-----------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| --work_dir          | CWD env setting                               | Working directory                                                                                                                                                                                                                            |
| --result_db         | Mandatory                                     | SQLite database holding all results                                                                                                                                                                                                          |
| --tis_ids           | Mandatory                                     | Analysis IDs for which an output file needs to be generated. Either an integer with one ID, a comma-separated list of IDs or 'all', leading to a FASTA/PEFF for all of the analysis IDs available                                            |
| --mslength          | 6                                             | Minimum sequence length allowed in the output FASTA/PEFF file                                                                                                                                                                                |
| --translation_db    | species_TIS_analysisID_transcripts.fasta/peff | Name of the output translation database file (FASTA/PEFF)                                                                                                                                                                                    |
| --var_file          | species_TIS_analysisID_transcripts_VAR.txt    | File to store SNP and indel information of sequences in the translation database                                                                                                                                                             |
| --tis_call          | Y                                             | Whether annotated TISs that did not pass the TIS calling algorithm are allowed in the output file. Possible options: Yes ('Y') or No ('N')                                                                                                   |
| --peff              | N                                             | Generate PEFF instead of FASTA format. Options: Yes ('Y') or No ('N')                                                                                                                                                                        |
| --mflag             | 4                                             | Flag describing the additional tasks that will be performed: 1 = remote BioMart mapping, 2 = local file BioMart mapping, 3 = sequence based blast mapping, 4 = no mapping, only redundancy removal, 5 = no mapping and no redundancy removal |
| --mapping           | Mandatory if mflag=1 or 2                     | Ensembl database to download BioMart mapped transcripts or local file to map data                                                                                                                                                            |
| --db_config_version | Mandatory if mflag=1                          | Ensembl database configuration version                                                                                                                                                                                                       |
| --external_ref      | Mandatory if mflag=1                          | External reference in biomart to map transcripts to                                                                                                                                                                                          |
| --blast_pgm         | ublast                                        | BLAST program to map based on sequence. Possible options: 'ublast' and 'blastp'                                                                                                                                                              |
| --blastdb           | Mandatory if mflag=3                          | The UBLAST or BLASTP formatted database                                                                                                                                                                                                      |
| --evalue            | 1e-10                                         | BLAST e-value                                                                                                                                                                                                                                |
| --min_blast_length  | 32                                            | Minimum sequence length to perform BLAST                                                                                                                                                                                                     |
| --identity          | 75                                            | Minimum alignment score for BLAST                                                                                                                                                                                                            |
| --coverage          | 30                                            | Minimum number of identical positions for BLAST                                                                                                                                                                                              |
| --word_size         | 3                                             | Minimum word size for BLAST                                                                                                                                                                                                                  |
| --gapopen           | 11                                            | Cost of gap open in BLAST                                                                                                                                                                                                                    |
| --gapextend         | 1                                             | Gap extension penalty for BLAST                                                                                                                                                                                                              |
| --matrix            | BLOSUM62                                      | Matrix for the BLAST search                                                                                                                                                                                                                  |

Redundancy will be default removed unless mflag 5 is selected or if a PEFF file is generated instead of a FASTA file.
Mapping can be done based on [BioMart](https://www.ensembl.org/biomart) (remote or with a local file) or sequence-based
using [UBLAST](https://drive5.com/usearch/manual/ublast_algo.html) or 
[BLASTP](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins). After mapping, inherent redundancy (also on sub-string
level) will be removed if desired. A VAR file will also be generated, next to the FASTA file in order to better describe 
which SNPs are linked to which SNP IDs.

Hereafter, a snippet of an example FASTA output file is given:

```
>generic|ENST00000647097_4_102760894_aTIS_1|ENSG00000109323 ATG True -1 0 +1 [ENST00000489071_20_35548863_ntr#ENST00000226578_4_102760894_aTIS_1]
MRLHLLLLLALCGAGTTAAELSYSLRGNWSICNGNGSLELPGAVPGCVHSALFQQGLIQD
SYYRFNDLNYRWVSLDNWTYSKEFKIPFEISKWQKVNLILEGVDTVSKILFNEVTIGETD
NMFNRYSFDITNVVRDVNSIELRFQSAVLYAAQQSKAHTRYQVPPDCPPLVQKGECHVNF
VRKEQCSFSWDWGPSFPTQGIWKDVRIEAYNICHLNYFTFSPIYDKSAQEWNLEIESTFD
VVSSKPVGGQVIIAIPKLQTQQTYSIELQPGKRIVELFVNISKNITVETWWPHGHGNQTG
YNMTVLFELDGGLNIEKSAKVYFRTVELIEEPIKGSPGLSFYFKINGFPIFLKGSNWIPA
DSFQDRVTSELLRLLLQSVVDANMNTLRVWGGGIYEQDEFYELCDELGIMVWQDFMFACA
LYPTDQGFLDSVTAEVAYQIKRLKSHPSIIIWSGNNENEEALMMNWYHISFTDRPIYIKD
YVTLYVKNIRELVLAGDKSRPFITSSPTNGAETVAEAWVSQNPNSNYFGDVHFYDYISDC
WNWKVFPKARFASEYGYQSWPSFSTLEKVSSTEDWSFNSKFSLHRQHHEGGNKQMLYQAG
LHFKLPQSTDPLRTFKDTIYLTQVMQAQCVKTETEFYRRSRSEIVDQQGHTMGALYWQLN
DIWQAPSWASLEYGGKWKMLHYFAQNFFAPLLPVGFENENTFYIYGVSDLHSDYSMTLSV
RVHTWSSLEPVCSRVTERFVMKGGEAVCLYEEPVSELLRRCGNCTRESCVVSFYLSADHE
LLSPTNYHFLSSPKEAVGLCKAQITAIISQQGDIFVFDLETSAVAPFVWLDVGSIPGRFS
DNGFLMTEKTRTILFYPWEPTSKNELEQSFHVTSLTDIY
>generic|ENST00000647097_4_102760894_aTIS_2|ENSG00000109323 ATG True -1 0 +1 [ENST00000226578_4_102760894_aTIS_2]
MRLHLLLLLALCGAGTTAAELSYSLRGNWSICNGNGSLELPGAVPGCVHSALFQQGLIQD
SYYRFNDLNYRWVSLDNWTYSKEFKIPFEISKWQKVNLILEGVDTVSKILFNEVTIGETD
NMFNRYSFDITNVVRDVNSIELRFQSAVLYAAQQSKAHTRYQVPPDCPPLVQKGECHVNF
VRKEQCSFSWDWGPSFPTQGIWKDVRIEAYNICHLNYFTFSPIYDKSAQEWNLEIESTFD
VVSSKPVGGQVIVAIPKLQTQQTYSIELQPGKRIVELFVNISKNITVETWWPHGHGNQTG
YNMTVLFELDGGLNIEKSAKVYFRTVELIEEPIKGSPGLSFYFKINGFPIFLKGSNWIPA
DSFQDRVTSELLRLLLQSVVDANMNTLRVWGGGIYEQDEFYELCDELGIMVWQDFMFACA
LYPTDQGFLDSVTAEVAYQIKRLKSHPSIIIWSGNNENEEALMMNWYHISFTDRPIYIKD
YVTLYVKNIRELVLAGDKSRPFITSSPTNGAETVAEAWVSQNPNSNYFGDVHFYDYISDC
WNWKVFPKARFASEYGYQSWPSFSTLEKVSSTEDWSFNSKFSLHRQHHEGGNKQMLYQAG
LHFKLPQSTDPLRTFKDTIYLTQVMQAQCVKTETEFYRRSRSEIVDQQGHTMGALYWQLN
DIWQAPSWASLEYGGKWKMLHYFAQNFFAPLLPVGFENENMFYIYGVSDLHSDYSMTLSV
RVHTWSSLEPVCSRVTERFVMKGGEAVCLYEEPVSELLRRCGNCTRESCVVSFYLSADHE
LLSPTNYHFLSSPKEAVGLCKAQITAIISQQGDIFVFDLETSAVAPFVWLDVGSIPGRFS
DNGFLMTEKTRTILFYPWEPTSKNELEQSFHVTSLTDIY
>generic|ENST00000647097_4_102760938_5UTR|ENSG00000109323 CTG NA -1 0 +1 [ENST00000645348_4_102760938_5UTR#ENST00000646727_4_102760938_5UTR#ENST00000644545_4_102760938_5UTR#ENST00000226578_4_102760938_5UTR#ENST00000642252_4_102760938_5UTR#ENST00000644159_4_102760938_5UTR#ENST00000644965_4_102760797_CDS#ENST00000646311_4_102760938_5UTR]
MPFDLSTSRWRGISRCASTCSCCSRCAVQAPPPRSSVTACVATGASAMGTARWSCPGRSL
AACTAPCSSRA
...
```

Another option is to generated a [PEFF file](http://www.psidev.info/peff), an extended FASTA format for proteogenomics 
recently devised by the HUPO-PSI. When the PEFF option is turned on, redundancy will not be removed but schematically 
included in the PEFF tags.

Hereafter, a snippet of an example FASTA output file is given:

```
# PEFF 1.0
# //
# DbName=genericDB
# DbSource=http://www.biobix.be/proteoform/
# DbVersion=2018-03-22
# HasAnnotationIdentifiers=true
# ProteoformDb=true
# Prefix=gen
# NumberOfEntries=4
# SequenceType=AA
# GeneralComment=Proteogenomics application to generate protein sequences from RIBO-seq
# //
>gen:ENST00000000412-1 \PName=mannose-6-phosphate receptor, cation dependent  \GName=mannose-6-phosphate receptor, cation dependent  \NcbiTaxId=9606 \TaxName=Homo sapiens \Length=277 \VariantSimple=(1:122|H) \VariantComplex=(2:1|120|M)(3:1|191|M)(4:1|235|M)  \Proteoform=(ENST00000000412_12_8946404_aTIS_pf1.0|1-277||canonical form)(ENST00000000412_12_8946404_aTIS_1_pf1.4|1-277|1|proteoform with SAV)(ENST00000000412_12_8943896_CDS_pf1.5|1-277|2|N-terminal truncation)(ENST00000000412_12_8943418_CDS_pf1.2|1-277|3|N-terminal truncation)(ENST00000000412_12_8942424_CDS_pf1.3|1-277|4|N-terminal truncation)(ENST00000000412_12_8943896_CDS_1_pf1.1|1-277|1,2|N-terminal truncation, proteoform with SAV)
MFPFYSCWRTGLLLLLLAVAVRESWQTEEKTCDLVGEKGKESEKELALVKRLKPLFNKSFESTVGQGSDTYIYIFRVCREAGNHTSGAGLVQINKSNGKETVVGRLNETHIFNGSNWIMLIYKGGDEYDNHCGKEQRRAVVMISCNRHTLADNFNPVSEERGKVQDCFYLFEMDSSLACSPEISHLSVGSILLVTFASLVAVYVVGGFLYQRLVVGAKGMEQFPHLAFWQDLGNLVADGCDFVCRSKPRNVPAAYRGVGDDQLGEESEERDDHLLPM
>gen:ENST00000000412-2 \PName=uORF of mannose-6-phosphate receptor, cation dependent  \GName=mannose-6-phosphate receptor, cation dependent  \NcbiTaxId=9606 \TaxName=Homo sapiens \Length=41 \VariantComplex=(1:1|11|M)  \Proteoform=(ENST00000000412_12_8949642_5UTR_pf2.0|1-41||canonical form)(ENST00000000412_12_8949612_5UTR_pf2.1|1-41|1|N-terminal truncation)
MGHSEALGGTLASETSSGTGVCGRPLAGPVSHPQKGIHWGF
>gen:ENST00000000412-3 \PName=Alternative reading frame product of mannose-6-phosphate receptor, cation dependent  \GName=mannose-6-phosphate receptor, cation dependent  \NcbiTaxId=9606 \TaxName=Homo sapiens \Length=16  \Proteoform=(ENST00000000412_12_8946336_CDS_pf3.0|1-16||canonical form)
MLADRRKNLRLGRRKG
>gen:ENST00000000412-4 \PName=Alternative reading frame product of mannose-6-phosphate receptor, cation dependent  \GName=mannose-6-phosphate receptor, cation dependent  \NcbiTaxId=9606 \TaxName=Homo sapiens \Length=56  \Proteoform=(ENST00000000412_12_8943516_CDS_pf4.0|1-56||canonical form)
MRSVAKSKIVSTSLRWIAAWPVHQRSPTSVWVPSYLSRLHHWLLFMLLGGSYTSDW
```
## Optional steps <a name="optional"></a>



## Copyright <a name="copyright"></a>

Copyright (C) 2014 G. Menschaert, J.Crappé, E. Ndah, A. Koch & S. Steyaert

Later updates: 2019 (in submission) S. Verbruggen, G. Menschaert, E. Ndah

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

## More information <a name="moreinformation"></a>

For more (contact) information visit [http://www.biobix.be/PROTEOFORMER](http://www.biobix.be/PROTEOFORMER)


