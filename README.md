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
    3. UTR simulation for Prokaryotes
    4. SRA parallel download
4. Main pipeline
    1. General quality check: fastQC
    2. Mapping
        norm bedgraphs
    3. General quality check: fastQC
    4. Specific ribosome profiling check: mappingQC
    5. Transcript calling
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

PROTEOFORMER is also available in galaxy: http://galaxy.ugent.be

## Dependencies <a name="dependencies"></a>

Proteoformer is built in Perl 5 and Python 2.7. All necessary scripts are included in this GitHub repository.

To prevent problems with missing dependencies, we included all necessary dependencies in a [Conda](https://conda.io/docs/) environment.
For more information about Conda installation, click [here](https://conda.io/docs/user-guide/install/index.html).

Once conda is installed, make sure to have the right channel order by executing following commands in the same order as listed here:

```
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```

Then you can install all dependencies based on the yml file in the dependency_envs folder of this GitHub repository with following command:

```conda env create -f Dependency_envs/proteoformer.yml```

To activate this new Conda environment:

```source activate proteoformer```

Some Perl packages are not included in Conda, so afterwards execute following script:

```perl install_add_perl_tools.pl```

### Additional environments for RiboZINB, SPECtre and SRA download

For some tools, we needed to construct separate environments with different versions of the underlying tools.

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
downloaded and constructed with the get_igenomes.py script in the 'Additional_tools' folder. For example:

```
python get_igenomes.py -v 82 -s human -d /path/to/dir -r -c 15
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


