Command line Script version of the PROTEOFORMER pipeline -  http://www.biobix.be/PROTEOFORMER


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

#################################################################
### Using the Ribosome Profiling Pipeline on the command line ###
#################################################################

1 Installing prerequisites: Perl packages, tool binaries and Pipeline scripts
-----------------------------------------------------------

The ribosome profiling pipeline is primarily written in Perl code, using different Perl specific packages. 
It is also dependant on a set of tool binaries which should all be installed on your system before the pipeline can execute all of its commands.

	Perl packages:
	--------------
		- Perl + v.5.10
		- DBI
		- DBD::SQLite
		- Parallel::ForkManager
		- Getopt::Long
		- Storable 'dclone'
		- Cwd	
		- BioPerl
		- LWP::UserAgent
		- XML::Smart

	Tool binaries:
	--------------
		- STAR, v2.4.2a or higher (https://code.google.com/p/rna-star/)
		- TOPHAT2, v.2.0.13 or higher  (http://tophat.cbcb.umd.edu/)
		- BLASTP ( ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) or USEARCH[*] (http://www.drive5.com/usearch/download.html)
		- R (http://www.r-project.org/)
		- SAMTOOLS, v. 1.19 or higher (http://sourceforge.net/projects/samtools/files/samtools/)
		- GATK (http://www.broadinstitute.org/gatk/download)
		- PICARD (http://sourceforge.net/projects/picard/files/picard-tools/)
		- SQLITE3 (http://www.sqlite.org/download.html)
		- FASTX toolkit, v.0.0.13 or higher  (http://hannonlab.cshl.edu/fastx_toolkit/)
		
		[*] The usearchX.Y.Z executatble should be renamed to "usearch"
		
	The tool binary paths should be included in the $PATH variable. 
	The path to the picard tool JAR files should be added to the $CLASSPATH variable. 
	This can be done globally by altering/adding these variables to the /etc/profile file in most linux distributions or by altering the $PATH for the galaxy user.  

	Ribosome Profiling pipeline scripts:
	------------------------------------

	1_mapping.pl
	2_assembly.pl
	filterSAMfile.pl
	gene_distribution.pl
	generate_transclation_db.pl
	metagenic_classification.pl
	metagenic_piecharts.R
	quality_plots.R
	ribo_translation.pl
	snp_calling
	snpIndexBuilder.pl
	splitVCFaltRecords.pl
	TIScalling_categorised.pl
	FlossProteoformer.pl
	run.sh

	All necessary scripts should be added to your working directory.

2 Installing prerequisites: Data dependencies
---------------------------------------------

	Igenomes:
	---------
	Install the igenomes in a specific igenomes_root folder and make sure it is accessible to the user running the pipeline.
	(http://support.illumina.com/sequencing/sequencing_software/igenome.ilmn)
	
		- gtf file: 					${IGENOMES_ROOT}/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf
		- reference whole genome sequence: 		${IGENOMES_ROOT}/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/
		- reference chromosome sequences: 		${IGENOMES_ROOT}/Mus_musculus/Ensembl/GRCm38/Sequence/Chromosomes/
		- PHIX-control sequences: 			${IGENOMES_ROOT}/Mus_musculus/Ensembl/GRCm38/Sequence/AbundantSequences/phix.fa
		- Chr size file: 				${IGENOMES_ROOT}/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/ChromInfo.txt
		- TopHat2 (Bowtie2) and STAR indexes: 		${IGENOMES_ROOT}/Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index
								${IGENOMES_ROOT}/Mus_musculus/Ensembl/GRCm38/Sequence/STARIndex (this folder and the STAR indexes are created 
								automatically when ther first STAR job is launched)
                - rRNA seq file	(custom made)                   ${IGENOMES_ROOT}/Mus_musculus/Ensembl/GRCm38/Sequence/AbundantSequences/species_rRNA.fa[***]
		- rRNA seq file	(custom made)			${IGENOMES_ROOT}/Mus_musculus/Ensembl/GRCm38/Sequence/AbundantSequences/species_tRNA.fa[***]
		- rRNA seq file	(custom made)			${IGENOMES_ROOT}/Mus_musculus/Ensembl/GRCm38/Sequence/AbundantSequences/species_sn-o-RNA.fa[***]
		
        [***] use 3-letter abbreviation for species: e.g. dme for Drosophila melanogaster, hsa for Homo sapiens, mmu for Mus musculus
              rRNA, tRNA, sn-o-RNA sequence fasta files have to be compiled. BioMart Ensembl allows you to download specific types of RNA sequences.

											
	Ensembl SQLite databases:
	-------------------------
	SQLite Ensembl DB with tables gene, coord_system, exon, exon_transcript, transcript, translation and seq_region. 
 		
 		- You can create these custom sqlite3 databases by using the mysqldump2sqlite3 perl script (mysql2sqlite3/mysqldump2sqlite3.pl).
 		  See script for a detailed how-to. Make sure your sqlite3 database follows the correct naming convention.
 		  		
 		  		Database naming convention:		ENS_[species short name]_[Ensembl annotation version].db
 		  		Species short name: 			Human:			hsa	
 		  							Mouse:			mmu
 		  							Fruitfly:		dme
 		  							Arabidopsis:	ath
 		  		
 		  		Examples: 				ENS_hsa_75.db
 		  							ENS_mmu_75.db
 		 							ENS_dme_75.db
 									ENS_ath_16.db
 		  					
 		- You can also download pre-formatted sqlite databases from the website (www.biobix.be/PROTEOFORMER) 
 		  or use the provided python script ENS_db.py (usage example: python ENS_db.py -v 78 -s human)
														
	Blast search databases:
	-----------------------
	The blast formated databases can be generated with the makeblastdb or the usearch -makeudb_ublast commands
	
		- blastp
			makeblastdb -in <protein sequence file in fasta> -out <database name> -dbtype prot
			
		- ublast
			usearch -makeudb_ublast <protein sequence file in fasta> -output <database name>
			
		- You can also download pre-formatted Blast search databases from the website (www.biobix.be/PROTEOFORMER)
		
3 Overview of the Ribosome Profiling Pipeline:
----------------------------------------------

Step 1: Mapping
Step 2: Transcript Calling
Step 3:	TIS Calling
Step 4:	SNP Calling
Step 5: Translation Assembly
Step 6: Translation Database
Step 7: Floss Calculation
Quality control 1: Metagenic Classification
Quality control 2: Gene Distribution

4 Step 1: Mapping
-----------------

	a) What it does
	
	This tool uses transcriptome mappers (STAR or TopHat2) to map RIBO-seq or RNA-seq next-generation sequencing
	reads against the reference genome using Ensembl annotation bundles (from the corresponding IGenomes:
	http://support.illumina.com/sequencing/sequencing_software/igenome.ilmn).
	The footprints are mapped both the NGS-reads of the untreated, cycloheximide (CHX), or emetine (EMT) sample
	and the NGS-reads of the treated: puromycin (PUR), harringtonine (HAR), or lactimidomycin (LTM)
	The footprint alignments are assigned to specific A site nucleotides by using the position and total 
	length of each alignment
	
	b) Input
	
	A FastQ file holding the Untreated, CHX, or EMT treated sample
	A FastQ file holding the PUR, HARR, or LTM treated sample
	A species and annotation-version specific sqlite Ensembl database (See 2 Ensembl SQLite databases).

	c) Command*
	
	./1_mapping.pl --name mESC --species mouse --ensembl 72 --cores 20 --readtype ribo --unique N --inputfile1 file1
	 --inputfile2 file2 --igenomes_root ${IGENOMES_ROOT} (--mapper STAR --adaptor CTGTAGGCACCATCAAT --clipper STAR 
	 --readlength 36 --phix N -- rRNA Y --snRNA N --tRNA N --splicing Y --mismatch 2 --maxmultimap 16
     --out_bg_s_untr bg_s_untr --out_bg_as_untr bg_as_untr --out_bg_s_tr bg_s_tr --out_bg_as_tr bg_as_tr 
     --out_sam_untr sam_untr --out_sam_tr sam_tr --out_sqlite sqliteDBName  --work_dir getcwd --tmpfolder $TMP)
	
	* See the specific Perl script for a more detailed explanation of all optional/mandatory arguments
	
	d) Output
	
	- The mapping output, placed under $work_dir/$mapper/$file 
	- BED and BEDGRAPH files holding the genome-wide ribosome-profile densities (only the first base of the A-site is used),
	placed under $work_dir/output
	- An SQLite database holding this genome-wide ribosome-profile information, next to some mapping statistics, 
	placed under $work_dir/SQLite
	
5 Step 2: Transcript Calling
----------------------------

	a) What it does
	
	Based on the sqlite database holding the experimental data from RNA-mapping (RIBO-seq) and using the species specific
	Ensembl annotation bundles (from the corresponding IGenomes) this tool determines the translation-level of transcripts
	(exon_level) for a given untreated/CHX/EMT-treated ribosome profiling data set (=DATA1 from RNA-mapping in 'Step 1: Mapping').
        
	b) Input
	
	- An sqlite database holding experimental data from the RIBO-seq-mapping (4 Step 1: Mapping).
    - A species and annotation-version specific sqlite Ensembl database (See 2 Ensembl SQLite databases).
	
	c) Command*
	
	./ribo_translation.pl (--work_dir getcwd --in_sqlite SQLite/results.db --ens_db SQLite/ens.db --tmp $TMP --out_sqlite SQLite/results.db)
	* See the specific Perl script for a more detailed explanation of all optional/manda
	
	d) Output
	
	An sqlite database holding all experimental data from the RNA-mapping (RIBO-seq) and transcript calling/translation tool. 
	After the transcript calling an extra table will be added (tr_translation):
	(NOTE: annotations as CCDS identifier and canonical transcript  are also inserted for further filtering of the transcripts)
        
    transcript_id  stable_id        chr         seq_region_id  seq_region_strand  seq_region_start  seq_region_end  read_counts  normalized_counts    biotype         exon_coverage  canonical   ccds        gene_stable_id 
    -------------  ---------------  ----------  -------------  -----------------  ----------------  --------------  -----------  -------------------  --------------  -------------  ----------  ----------  ---------------
    2427712        ENST00000445884  1           27511          1                  10002981.0        10010032.0      3.0          0.00681818181818182  sense_intronic  Yes            Yes         No          ENSG00000228150	
    ...	

6 Step 3: TIS Calling
---------------------

	a) What it does
	
	This script searches for all possible TISes within known Ensembl transcripts. Only TISes passing a number of specified arguments are being withheld.
	The minimal profile coverage and Rltm - Rchx parameter setting can be specified for each annotation class (5'UTR, annotated TIS (aTIS), 
	coding sequence (CDS), 3'UTR, no translation(no_trans))

    Local maximum: The newly identified TIS should have the maximal number of reads within a down- and upstream window of x basepairs.
    Minimal profile coverage: The minimum number of ribosome profiles on a TIS-site (after combining the reads, because of subcodon specificity, on ATG or near cognate start positions).
    Rltm - Rchx: Value calculated according to the function mentioned below. The TIS should have a value equal or higher to the parameter setting.

    Rk = (Xk/Nk) x 10 (k = LTM, CHX), Xk number of reads on that position in data k, Nk total number of reads for transcript.

	b) Input
	
	- An sqlite database holding experimental data from the RIBO-seq-mapping and transcript calling (5 Step 2: Transcript Calling).
	
	c) Command*
	
	./TIScalling_categorised.pl --sqlite SQLite/results.db --cores 22 --local_max 1 --R_aTIS 0.01 --min_count_aTIS 5 --R_5 0.05 --min_count_5 10
	--R_CDS 0.15 --min_count_CDS 15 --R_3 0.05 --min_count_3 10 --R_ntr 0.05 --min_count_ntr 10
	
	* See the specific Perl script for a more detailed explanation of all optional/mandatory arguments
	
	d) Output
	
	An sqlite database holding all experimental data from the RNA-mapping (RIBO-seq), transcript translation and TIScalling tool. After TIScalling 2 extra tables will be added and/or updated.
	An extra TIS_id table with all TISses identified using the specified parameters.

    transcript_id   |   stable_id               |   biotype                 |   chr |   strand  |   start       |   dist_to_transcript_start    |   dist_to_aTIS    |   annotation      |   aTIS_call   |start_codon     	|   peak_shift  |   count   |   Rltm_min_Rchx
    ---------------- --------------------------- --------------------------- ------- ----------- --------------- ------------------------------- ------------------- ------------------- ------------------- --------------- --------------- ----------- --------------------
    356625          |   ENSMUST00000161973      |   processed_pseudogene    |   1   |   1       |   143683266   |   37                          |   NA              |   no_translation  |   NA          |   AGG             |   +1 0 -1     |   151.0   |   0.822453056367707
    ...

	A TIS_overview table with an overview of the parameters used in the different TIScallings within the results databse. 
	
    ID  |   local_max   |   min_count_aTIS  |   R_aTis  |   min_count_5UTR  |   R_5UTR  |   min_count_CDS   |   R_CDS   |   min_count_3UTR  |   R_3UTR  |   min_count_no_trans  |   R_no_trans
    ---- --------------- ------------------- ----------- ------------------- ----------- ------------------- ----------- ------------------- ----------- ----------------------- -------------
    1   |   1           |   5               |   0.01    |   10              |   0.05    |   15              |   0.15    |   10              |   0.05    |   10                  |   0.05
    ...

	
7 Step 4: SNP Calling
---------------------

	a) What it does
	
	This script uses samtools to call variants in mapped next-generation sequencing reads.
	The user has the option to compare the mapping and snp calling results to SNP data from the Sanger Institute.
	When you use the appropriate input arguments, duplicate reads in the input file will be remove with picard. 
	Next, the variants are called by chromosome and the results are stored in a separate VCF file for each chromosome. 
	The benefit of calling variants per chromosome is that the tool can easily be multithreaded and can be run on up 
	to x different cores (where x = number of chromosomes). When the SNP calling is finished, all the VCF files are 
	combined in one single VCF file. This file is filtered based on the parameters high_af, lower_af and upper_af
	(allelic frequency) and the results are written to an SQLite database (table snp_samtools).
	
	b) Input
	
	A sam file with the mapped next-generation sequencing reads and an SQLite database where the results will be added to.
	
	c) Command*
	
	snp_calling -e experimentName -r pathToMappedReads -c pathToChromosomeSequencesFolder -o organism -t 1 
	--mincoverage 3 --maxcoverage 100 --high_af 0.95 --lower_af 0.3 --upper_af 0.7
	
	* See the specific bash script for a more detailed explanation of all optional/mandatory arguments
	
	d) Output
	
	This tool adds a table (snp_samtools), which contains the SNP calling result, to the SQLite database.
	The column "new" is added when the results are compared to SNPdb. A "y" (yes) indicates that a variant was new, 
	ie NOT found in SNPdb. A "n" (no) means not new (ie found in SNPdb) and an "m" (mismatch) is used for those 
	mismatches in the mapped reads that were found in SNPdb. There is no allelic frequency information for the
	"m" variants, so this has been set to 0.0.

8 Step 5: Translation assembly
------------------------------

	a) What it does
	
	This script assembles all translation products based on all information derived from the ribosome profiling.
	
	b) Input
	
	An sqlite database holding experimental data from previous steps: i.e. 1.mapping, 2.transcript calling, 
	3.TIS calling, 4. SNP calling (optional)
	
	c) Command*
	
	./2_assembly.pl  --sqliteRES sqlite_results_DB --sqliteENS sqlite_ensembl_DB --tis_ids list_of_analysis_ids 
	(--snp NO --dir /path/ --tmp /path/to/tmp/ --localmax 1 --mincount_aTIS 10 --R_aTIS .05 and the other mincount/R values
	 --out_sqlite $out_sqlite)
	
	* See the specific Perl script for a more detailed explanation of all optional/mandatory arguments
	
	d) Output
	
	An sqlite database holding all experimental data resulting from the 1.mapping, 2.transcript calling, 
	3.TIS calling, 4. SNP calling (optional) and 5.translation product assembly. After assembly step an
	extra table will be added:
	An extra TIS_id_transcripts table with all translation products.

	tr_stable_id        chr         strand      start       start_codon  dist_to_transcript_start  dist_to_aTIS  annotation      peak_shift  count       Rltm_min_Rchx      coverage    FPKM        SNP_        tr_seq                                                                                                                                                                                                                                                               aa_seq
	------------------  ----------  ----------  ----------  -----------  ------------------------  ------------  --------------  ----------  ----------  -----------------  ----------  ----------  ----------  ---------------------------------------------------------------------------------------------                                                                                                                                                                        ----------
	ENSMUST00000134291  1           -1          178335900   GTG          414                       NA            no_translation  -1 0        14.0        0.148350964210667                                      GTGGTTTGTCTTGATACTTGTAAGTGAGGCGGACTTTCCCGCTCTTTGCTAGTTTTCAAGTGCAGTGTTAACATTGGTTTTAAGTGTTAGGAAAAATGTGTTGATCTTAAAGTGGAAAAAACCAGTTAAGAAACTTATTTTAAAGCAGGAAGGTAGTAGCTTAATGCTCAGGCTACAGTAGTTTAATATTAATACCACTTCAGTTGTATTGATTTAACAAGCTCTAGTAATTGGTAATATTTACGAGTTGAAACTTCTC  MVCLDTCK*
	...
	
	
9 Step 6: Translation database
------------------------------

	a) What it does
	
	This tool generates a non redundant custom protein sequence database from transclation products generated
	based on known Ensemble transcripts. It removes all duplicate and sub-sequences of the generated transclation 
	products taking into considerations the highrarchy in the table below. It performs a Local BLAST search
	(using blastp or ublast) against a canonical protein database to map the protein sequences to already 
	annotated proteins.
	
	The annotations are ranked according to the table below with 1 the most important.

    Annotation                                     		 |  Rank
	--------------------------------------------------------------
	aTIS (annotated TIS)                                 |    1
	5UTR (5'untransalted region)                         |    2
	CDS (transclation starting within a coding sequence) |    3
	ntr (translation evidence at pseudo-genes)           |    4
	3UTR (3'untransalted region)                         |    5
	
	b) Input
	
	- A SQLite database with a table containing the tranacript ID, annotation, chromosomes, start positions and start codon.
	- A Blast or USearch formated canonical protein database (e.g. SWISSPROT). 
	- An Ensembl database name to download biomart mapped Ensembl to Swissprot IDs or a comma seperated text file containing
	 the mapped IDS in the following structure.

	Ensemble_transcript_ID,Swissprot_Accession,Gene_name
	ENSMUST00000000033,P09535,IGF2_MOUSE
	
	c) Command*
	
	perl generate_db.pl -blast_db /path/to/udb -sqlite_db /path/to/db -tis_ids 1 -blast_pgm ublast -mflag 0 -external_ref uniprot_swissprot_accession -mapping_db mmusculus_gene_ensembl -num_threads 3 -work_dir work_directory
	
	* See the specific Perl script for a more detailed explanation of all optional/mandatory arguments
	
	d) Output
	
	A Fasta file of non redundant translation products database.
	
10 Quality Control 1: Metagenic Classification
----------------------------------------------

	a) What it does
	
	This tool performs a metagenic classification (annotation) of the mapped RIBO-read positions firstly on i) transcripts defined in Ensembl
	as 'protein_coding' (biotype) (i.e. '5'UTR','3'UTR','Exon','Intron'), next ii) all other transcripts (i.e. 'Other biotypes') and finally 
	iii) intergenic regions ('Intergenic')
	Next, in a second step a classification is made of the RIBO-seq read positions mapping on transcripts with biotypes other than 
	'protein_coding' (the 'Other biotypes' of step 1).
	
	b) Input
	
	An sqlite database holding experimental data from the RNA-mapping (RIBO-seq).
	
	c) Command*
	
	./metagenic_classification.pl --cores 3 --treated untreated (or 'treated') (--work_dir getcwd --in_sqlite SQLite/results.db 
	--output_file run_name_annotation)
	
	* See the specific Perl script for a more detailed explanation of all optional/mandatory arguments
	
	d) Output
	
	Two tables (TAB-separated files; one for all RIBO-seq read positions and one for the mapped positions that fall into non-protein 
	coding transcripts) are outputted with the obtained counts per functional region (per chromosome):

	TABLE 1:

	chr |   ribo    |   exon    |   5utr    |   3utr    |   intron  |   non_protein_coding      |       intergenic
	--------------------------------------------------------------------------------------------------------
	10  |   251511  |   195295  |   8517        |   2412    |   8496    |       10136                   |       26655
	1   |   317695  |   237148  |   10252   |   2366    |   12974   |           20554                   |       34401
	...
	TABLE 2:

	chr |   non_protein_coding  |       lincRNA |       miRNA   |       misc_RNA        |       nonsense_mediated_decay | processed_pseudogene  | . .
	-----------------------------------------------------------------------------------------------------------------------------
	10  |               10136                   |       284             |       34              |               2               |                       8                               |               7560                    | . .
	1   |               20554                   |       773             |       54              |               16              |                       31                              |               15265                   | . .
	...

	To summarize, for each table a pie chart (pdf) is generated.
	
11 Quality Control 2: Gene Distribution
---------------------------------------

	a) What it does
	
	This tool determines in which genes the RIBO-seq reads fall, and determines the total read count for these genes.
	
	b) Input
	
	An sqlite database holding experimental data from the RNA-mapping (RIBO-seq).
	
	c) Command*
	
	./gene_distribution.pl --treated untreated (or 'treated') (--work_dir getcwd --in_sqlite SQLite/results.db --output_file run_name_genedistribution)
	
	* See the specific Perl script for a more detailed explanation of all optional/mandatory arguments
	
	d) Output
	
	The results are outputted in the form of a table (TAB-separated file with 2 columns: GeneID and read_count):

        GeneID          |   read_count
	------------------------------------
	ENSMUSG00000064367      |       4216
	ENSMUSG00000065947      |       2051
	ENSMUSG00000064342      |       75
	ENSMUSG00000064368      |       4279
	ENSMUSG00000064365      |       64
	ENSMUSG00000064352      |       71
	ENSMUSG00000064349      |       38
	ENSMUSG00000064347      |       301
	ENSMUSG00000064355      |       96
	ENSMUSG00000064336      |       73643
	...

	Additionally, to visualize this table, 3 overall gene abundance plots (pdf) are generated:

	-cumulative gene distribution
	-gene density
	-ranked gene abundance


12 Floss Calculation:
---------------------

	a) What it does

    This tool calculates:
    
        (1) fragment length reference fractions based on protein-coding transcripts
        (2) cutoff values as function of the amount of reads
        (3) a length distribution for each possible translation product
         
    Based on these results, a FLOSS score for each putative transcript can be calculated. With this FLOSS score and the cutoff values, the coding potential of each possible product can be assessed.

	b) Input
        
    An sqlite database holding experimental data from at least following previous steps: mapping, transcript calling, TIS calling, SNP calling (optional) and translation product assembly.
        
    The TIS ids for which a FLOSS score needs to be calculated. It can be a sequence of ids seperated by commas or it can be "all".

	c) Command*
	
	./FlossProteoformer.pl --sqlite SQLite/results.db  --tis_ids 1
	
	* See the specific Perl script for a more detailed explanation of all optional/mandatory arguments


	d) Output
        
    An sqlite database holding all experimental data resulting from the mapping, transcript calling, TIS calling, SNP calling (optional), translation product assembly and FLOSS score calculation. Afterwards three extra tables will be added:
        
    A table FLOSS_ref_fractions with reference fractions amongst each possible fragment length
    
            RPF |   fraction
            ----  ----------
            26  |   0.0248963500775222
            

    A table FLOSS_cutoff containing the FLOSS score cutoff value for each amount of reads
    
            nreads  |   score
            -------- ------------
            1       |   0.63122634070155
            
            
    A table TIS_(analysisID)_transcripts_FLOSS containg the amount of Ribo reads, the FLOSS score and the classification for each TIS ID. The TIS ID is a combination of the transcript ID and the supposed translation initiation site (multiple initiation sites possible for one transcript). The table ID is a chronological key specific for the table. A classification being "Good" means that the reading frame is probably coding. A classification being "Extreme" means that the reading frame is probably non-coding. The amount of reads can also be 0 ('no reads') or out of the range of the cutoff ('out of cutoff range'). Then, the calssification will be impossible.

            tableID  |   TIS_ID            |   nreads  |     FLOSS              |  classification
            --------- ---------------------- ----------  ------------------------ -----------------
            1        |   426972_121413981  |   978     |    0.0959806966454703  |     Good
            2        |   423089_34772322   |   10468   |    0.0368119161096991  |     Good



13 Executing the Proteoformer pipeline: 
---------------------------------------
	The Proteoformer scripts can be executed individually from step 1 through step 7 as described above or through a 
	wrapper bash script run_proteoformer.sh. 
	
	The input arguments to the bash script are described in the individual perl scripts.













			