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
        2. RiboZINB
    6. ORF calling
        1. PROTEOFORMER
        2. PRICE
        3. SPECtre
    7. Fasta file generation
5. Optional steps
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

PROTEOFORMER is also available in galaxy: http://galaxy.ugent.be

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

```conda env create -f Dependency_envs/proteoformer.yml```

This installs a new Conda environment in which all needed Conda dependencies are installed and available, including Perl
and Python. If not mentioned otherwise, all tools of the PROTEOFORMER pipeline should be executed in this environment.
 To activate this new Conda environment:

```source activate proteoformer```

Some Perl packages are not included in Conda, so after installation and first activation of the new environment,
 execute following script:

```perl install_add_perl_tools.pl```

If you want to exit the proteoformer Conda environment:

```source deactivate```

### Additional environments for RiboZINB, SPECtre and SRA download

For some tools, we needed to construct separate environments with different versions of the underlying tools. For all 
the other tools, the proteoformer environment is used.

#### RiboZINB

```
cconda env create -f Dependency_envs/ribozinb.yml
ssource activate ribozinb
```

#### SPECtre

```
cconda env create -f Dependency_envs/spectre.yml
ssource activate spectre
```

#### SRA download

```
cconda env create -f Dependency_envs/download_sra_parallel.yml
ssource activate download_sra_parallel
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

| Species                  | Input value species argument |
|--------------------------|------------------------------|
| Homo sapiens             | human                        |
| Mus musculus             | mouse                        |
| Drosophila melanogaster  | fruitfly                     |
| Saccharomyces cerevisiae | yeast                        |
| Caenorhabditis elegans   | c.elegans                    |

The Ensembl database for SL1344 (Salmonella) is available under request.

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

```fastqc yourdata.fastq -t 20```

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

### Rule-based transcript calling <a name="rulebased"></a>

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

| Argument     | Default                        | Description                                         |
|--------------|--------------------------------|-----------------------------------------------------|
| --work_dir   | CWD env setting                | The working directory                               |
| --tmp        | work_dir/tmp                   | The temporary files folder                          |
| --in_sqlite  | SQLite/results.db              | The SQLite results database from previous steps     |
| --out_sqlite | The same as in_sqlite argument | The SQLite results database to store all results in |
| --ens_db     | Mandatory                      | The Ensembl database with annotation info           |
| --help       |                                | Generate help message                               |

Output table structure:

| transcript_id | stable_id       | chr | seq_region_id | seq_region_strand | seq_region_start | seq_region_end | read_counts | normalized_counts | biotype        | exon_coverage | canonical | ccds      | gene_stable_id  | FPKM             |
|---------------|-----------------|-----|---------------|-------------------|------------------|----------------|-------------|-------------------|----------------|---------------|-----------|-----------|-----------------|------------------|
| 196519        | ENST00000371471 | 20  | 131538        | -1                | 53567065         | 53593839       | 585         | 0.187680461982676 | protein_coding | Yes           | Yes       | CCDS13443 | ENSG00000171940 | 10.6593122808274 |
| ...           | ...             | ... | ...           | ...               | ...              | ...            | ...         | ...               | ...            | ...           | ...       | ...       | ...             | ...              |


## Copyright <a name="copyright"></a>

Copyright (C) 2014 G. Menschaert, J.Crappé, E. Ndah, A. Koch & S. Steyaert

Later updates: S. Verbruggen, G. Menschaert, E. Ndah

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

For more (contact) information visit http://www.biobix.be/PROTEOFORMER


