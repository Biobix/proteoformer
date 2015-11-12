README Quality (SandraS)

# General prerequisites
----------------------

- Perl v.5.10
- DBI package
- Parallel ForkManager package
- SQLite
- R

- Ribo-seq reads parsed to SQLite results database (eg. HCT_HS_STAR_N_results.db)
- ENSEMBL databases (SQLite format)
	-> tables: seq_region, gene, transcript, exon_transcript, exon, translation, coord_system_id
	* Mouse: ensembl_mus_musculus_core_66_37 (ENS_mmu_66.db)
	* Human: ensembl_homo_sapiens_core_66_37 (ENS_hsa_66.db)

*******************
EXPRESSION ANALYSIS
*******************
#Constructing transcript translation SQLITE tables
-------------------------------------------------------------
>> ribo_translation.pl
	> Input variables:
		my $work_dir = $ARGV[0]; #Working directory (eg. "/data/RIBO_runs/RIBO_BioBix_Pipeline/")
		my $species = $ARGV[1]; #Species name (English, for now only "mouse" or "human")
		my $version = $ARGV[2]; #Ensembl Assembly (human: "66" or "70", mouse: "66", "70", or "72")
		my $ribo_run = $ARGV[3]; #Name of sequencing run, MUST BE THE SAME AS first part of the name of SQLITE results db (eg."HCT_HS_STAR_N")
		my $lane = $ARGV[4]; #Lane you want to analyze (eg. "Lane2" or "lane2"): depends on how it is called in the table names of SQLITE results db
		my $mapper = $ARGV[5]; #Name of the used mapper (which is also defined in the name of the SQLITE db) (eg. "STAR")
		my $cores = $ARGV[6]; #Number of cores you want to use
		
	> Output:
 		One SQLITE tables in results db:
 		- '$ribo_run'_'$lane'_transcript_translation
 		
 	> Usage:
 	$ perl ribo_translation.pl 'workdir' 'species_name(english)' 'ENSEMBL_version' 'ribo_run' 'lane' 'mapper' 'cores'
		

*******************
QUALITY ASSESSMENT
*******************
	
#1. FASTQC
-----------
>> Usage:
Command line (Aramis): $ fastqc --outdir=/some/other/dir/ lane1.fastq

#2. Gene distribution and dynamic range of ribo-seq reads
----------------------------------------------------------

>> 2.1 gene_distribution.pl
 	> Input variables: 
 		my $work_dir = $ARGV[0]; #Working directory (eg. "/data/RIBO_runs/RIBO_BioBix_Pipeline/")
 		my $species = $ARGV[1]; #Species name (English, for now only "mouse" or "human")
		my $version = $ARGV[2]; #Ensembl Assembly (human: "66" or "70", mouse: "66", "70", or "72")
		my $ribo_run = $ARGV[3]; #Name of sequencing run, MUST BE THE SAME AS first part of the name of SQLITE results db (eg."HCT_HS_STAR_N")
		my $lane = $ARGV[4]; #Lane you want to analyze (eg. "Lane2" or "lane2"): depends on how it is called in the table names of SQLITE results db
		my $mapper = $ARGV[5]; #Name of the used mapper (which is also defined in the name of the SQLITE db) (eg. "STAR")
		my $outdir = $ARGV[6]; #Output directory (eg. /data/sandras/RIBO-seq/quality/hiseq)
		
 	> Output:
 		one text file with the gene_ids, the count of the ribo-seq reads on that gene
 		- '$ribo_run'_'$lane'_gene_distribution.txt
 		
 	> Usage:
 	$ perl gene_distribution.pl 'workdir' 'species_name(english)' 'ENSEMBL_version' 'ribo_run' 'lane' 'mapper' 'outdir'
 	
>> 2.2 quality_plots.R
	> Input variables:
		workdir <- as.character(args[1]) #Directory input/output files (= same as outputfile of 2.1)
		ribo_run <- as.character(args[2]) #Name of sequencing run, MUST BE THE SAME AS first part of the name of SQLITE results db (eg."HCT_HS_STAR_N")
		lane <- as.character(args[3]) #Lane you want to analyze (eg. "Lane2" or "lane2"): depends on how it is called in the table names of SQLITE results db
		
	> Output:
		three figures (.pdf), namely the plots of the gene abundance, the cumulative gene abundance and the density
		- '$ribo_run'_'$lane'_gene_abundance.pdf
		- '$ribo_run'_'$lane'_gene_cumulative_abundance.pdf
		- '$ribo_run'_'$lane'_gene_density.pdf
		
	> Usage:
	$ Rscript quality_plots.R 'workdir' 'ribo_run' 'lane'

#3. Functional annotation of ribo-seq reads
--------------------------------------------
	
>> 3.1 metagenic_classification.pl
	> Input variables:
		my $work_dir = $ARGV[0]; #Working directory (eg. "/data/RIBO_runs/RIBO_BioBix_Pipeline")
		my $species = $ARGV[1]; #Species name (English, for now only "mouse" or "human")
		my $version = $ARGV[2]; #Ensembl Assembly (human: "66" or "70", mouse: "66", "70", or "72")
		my $ribo_run = $ARGV[3]; #Name of sequencing run, MUST BE THE SAME AS first part of the name of SQLITE results db (eg."HCT_HS_STAR_N")
		my $lane = $ARGV[4]; #Lane you want to analyze (eg. "Lane2" or "lane2"): depends on how it is called in the table names of SQLITE results db
		my $mapper = $ARGV[5]; #Name of the used mapper (which is also defined in the name of the SQLITE db) (eg. "STAR")
		my $cores = $ARGV[6]; #Number of cores you want to use
		my $outdir = $ARGV[7]; #Output directory (eg. /data/sandras/RIBO-seq/quality/hiseq)
		
 	> Output:
 		two text files in the map 'quality/ribo_table' for the corresponding ribo-table with the reads that map on coding transcripts 
 		and non-coding transcripts together with their functional annotation
 		- '$ribo_run'_'$lane'_annotation_coding.txt
 		- '$ribo_run'_'$lane'_annotation_noncoding.txt
 		
 	> Usage: 
 	$ perl metagenic_classification.pl 'workdir' 'species_name(english)' 'ENSEMBL_version' 'ribo_run' 'lane' 'mapper' 'cores' 'outdir'
 	
 >> 3.2 metagenic_piecharts.R
 	> Input variables:
		workdir <- as.character(args[1]) #Directory input/output files (= same as outputfile of 3.1)
		ribo_run <- as.character(args[2]) #Name of sequencing run, MUST BE THE SAME AS first part of the name of SQLITE results db (eg."HCT_HS_STAR_N")
		lane <- as.character(args[3]) #Lane you want to analyze (eg. "Lane2" or "lane2"): depends on how it is called in the table names of SQLITE results db

	> Output: 
		two figures (.pdf), namely the functional region annotation pie charts
		- '$ribo_run'_'$lane'_annotation_coding.pdf
		- '$ribo_run'_'$lane'_annotation_noncoding.pdf
		
 	> Usage:
 	 $ Rscript metagenic_piecharts.R 'workdir' 'ribo_run' 'lane'
 	 

 	

