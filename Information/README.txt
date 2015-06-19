Provides galaxy tools for the Ribosome profiling pipeline -  http://www.biobix.be/TODO
The ribosome profiling pipeline cannot be installed automatically from a toolshed

##############################################################
## MANUAL INSTALLATION FOR THE RIBOSOME PROFILING PIPELINE: ##
##############################################################

1  Galaxy specific configuration:
---------------------------------
To manually install the ribosome profiling pipeline a set of configuration files, wrapper scripts, datatypes and other configuration and 
location files need to be added to specific folders of your galaxy system. These are listed below and can be found within this directory.

	- Add configuration files (.xml), perl code (.pl), shell code (.sh) and location files (.loc) from tools/ribo_prof_pipeline to your galaxy installation
		
		Overview:
		---------
			- tools/ribo_prof_pipeline/1_mapping.pl
			- tools/ribo_prof_pipeline/1_mapping.xml
			- tools/ribo_prof_pipeline/2_assembly.pl
			- tools/ribo_prof_pipeline/2_assembly.xml
			- tools/ribo_prof_pipeline/blast_db.loc[*]
			- tools/ribo_prof_pipeline/ENS_db.loc[*]
			- tools/ribo_prof_pipeline/filterSAMfile.pl
			- tools/ribo_prof_pipeline/gene_distribution.pl
			- tools/ribo_prof_pipeline/gene_distribution.xml
			- tools/ribo_prof_pipeline/generate_translation_db.pl
			- tools/ribo_prof_pipeline/generate_translation_db.xml
			- tools/ribo_prof_pipeline/igenomes.loc[*]
			- tools/ribo_prof_pipeline/metagenic_classification.pl
			- tools/ribo_prof_pipeline/metagenic_classification.xml
			- tools/ribo_prof_pipeline/metagenenic_piecharts.R
			- tools/ribo_prof_pipeline/quality.R
			- tools/ribo_prof_pipeline/ribo_translation.pl
			- tools/ribo_prof_pipeline/ribo_translation.xml
			- tools/ribo_prof_pipeline/snp_calling
			- tools/ribo_prof_pipeline/snp_calling.xml
			- tools/ribo_prof_pipeline/snpdb.loc[*]
			- tools/ribo_prof_pipeline/snpIndexBuilder.pl
			- tools/ribo_prof_pipeline/splitVCFaltRecords.pl
			- tools/ribo_prof_pipeline/tis_calling.xml
			- tools/ribo_prof_pipeline/TIS_overview.pl
			- tools/ribo_prof_pipeline/tis_overview.xml
			- tools/ribo_prof_pipeline/Tiscalling_categorised.pl
			
		[*] See point 4: "Tool data location files" for more information on how to use the location (.loc) files
		
	- Add the datatype definition file from lib/galaxy/datatypes/ribo_prof.py to your galaxy installation
	
	- Add the following import line to:  lib/galaxy/datatypes/registry.py
			import riboprof # added for Ribosome profiling pipeline
			
	- Add datatypes between the <registration>   </registration> tags in:  datatypes_conf.xml
		 	
		 <!-- Start Ribo_prof_pipeline Datatypes -->
    			<datatype extension="tis" type="galaxy.datatypes.riboprof:TisOverview"/>
   			<datatype extension="sqlitedb" type="galaxy.datatypes.riboprof:SqliteDb" />
    		<!-- End Ribo_prof_pipeline Datatypes -->
    		
    - Copy the ribo_prof_tool_conf.xml file to your galaxy installation and add it to the universe_wsgi.ini file:
    		
    		# Tool config files, defines what tools are available in Galaxy.
		# Tools can be locally developed or installed from Galaxy tool sheds.
		tool_config_file = ribo_prof_tool_conf.xml,tool_conf.xml
		
	  or add the <section name="Ribosome Profiling pipeline" id="ribo_prof_pipe">...</section> from the ribo_prof_tool_conf to your tool_conf.xml
	  
	- Reorganize the integrated_tool_panel.xml and add the following lines where you want them to appear in the galaxy tool column.
		
    		<section id="ribo_prof_pipeline" name="Ribosome Profiling pipeline" version="">
        		<tool id="mapping" />
        		<tool id="transcript_calling" />
        		<tool id="tis_calling" />
        		<tool id="tis_overview" />
        		<tool id="snp_calling1" />
        		<tool id="2_assembly" />
        		<tool id="Ensembl2Swissprot Mapping" />
        		<tool id="metagenic_classification" />
        		<tool id="gene_distribution" />
    		</section>

2 Installing prerequisites: Perl packages and tool binaries
-----------------------------------------------------------

The ribosome profiling pipeline is primarily written in Perl code, using different Perl specific packages. 
It is also dependant on a set of tool binaries which should all be installed on your galaxy system before the pipeline can execute all of its commands.

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
	------
		- STAR (https://code.google.com/p/rna-star/)
		- TOPHAT2 (http://tophat.cbcb.umd.edu/)
		- BLASTP[**] ( ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) or USEARCH[***] (http://www.drive5.com/usearch/download.html)
		- R (http://www.r-project.org/)
		- SAMTOOLS 1.19 or above (http://sourceforge.net/projects/samtools/files/samtools/)
		- GATK (http://www.broadinstitute.org/gatk/download)
		- PICARD (http://sourceforge.net/projects/picard/files/picard-tools/)
		- SQLITE3 (http://www.sqlite.org/download.html)
		- FASTX toolkit (http://hannonlab.cshl.edu/fastx_toolkit/)
		
		[**]
		[***] The usearchX.Y.Z executatble should be renamed to "usearch"
		
The tool binary paths should be included in the $PATH variable. 
The path to the picard tool JAR files should be added to the $CLASSPATH variable. 
This can be done globally by altering/adding these variables to the /etc/profile file in most linux distributions or by altering the $PATH for the galaxy user.  
		
3 Installing prerequisites: Data dependencies
---------------------------------------------

	Igenomes:
	---------
	Install the igenomes in a specific igenomes_root folder and make sure it is accessible to the galaxy user.
	(http://support.illumina.com/sequencing/sequencing_software/igenome.ilmn)
	
		- gtf file: 					IGENOMES_ROOT/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf
		- reference whole genome sequence: 		IGENOMES_ROOT/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/
		- reference chromosome sequences: 		IGENOMES_ROOT/Mus_musculus/Ensembl/GRCm38/Sequence/Chromosomes/
		- PHIX-control sequences: 			IGENOMES_ROOT/Mus_musculus/Ensembl/GRCm38/Sequence/AbundantSequences/phix.fa
		- Chr size file: 				IGENOMES_ROOT/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/ChromInfo.txt
		- TopHat2 (Bowtie2) and STAR indexes: 		IGENOMES_ROOT/Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index
								IGENOMES_ROOT/Mus_musculus/Ensembl/GRCm38/Sequence/STARIndex
		- rRNA seq file	(custom made)			IGENOMES_ROOT/Mus_musculus/Ensembl/GRCm38/Sequence/AbundantSequences/rRNA_species.fa[****]
		
		[****] rRNA seq fasta files have to be compiled, you can download pre-formatted rRNA sequence fasta files from the website (www.biobix.be/TODO)
												
	Ensembl SQLite databases:
	-------------------------
	SQLite Ensembl DB with tables gene, coord_system, exon, exon_transcript, transcript, translation and seq_region. 
 		
 		- You can create these custom sqlite3 databases by using the mysqldump2sqlite3 perl script (mysql2sqlite3/mysqldump2sqlite3.pl).
 		  See script for a detailed how-to. Make sure your sqlite3 database follows the correct naming convention.
 		  		
 		  		Database naming convention:		ENS_[species short name]_[Ensembl annotation version].db
 		  		Species short name: 			Human:		hsa	
 		  							Mouse:		mmu
 		  							Fruitfly:	dme
 		  							Arabidopsis:	ath
 		  		
 		  		Examples: 				ENS_hsa_70.db
 		  							ENS_mmu_72.db
 		 							ENS_dme_74.db
 		  							ENS_ath_16.db
 		  					
 		- You can also download pre-formatted sqlite databases from the website (www.biobix.be/TODO)
														
	Blast search databases:
	-------------------------
	The blast formated databases can be generated with the makeblastdb or the usearch -makeudb_ublast commands
	
		- blastp
			makeblastdb -in <protein sequence file in fasta> -out <database name> -dbtype prot
			
		- ublast
			usearch -makeudb_ublast <protein sequence file in fasta> -output <database name>
			
		Update the path to the blast databases in the ublast_db.loc and blastp_db.loc files.
		
4 Tool data location files:
---------------------------

Reorganize the tool_data_table_conf.xml and add the following lines in between the <tables>  </tables>. 

	<!-- Location of Ribo-seq tool files  -->
    <table name="ribo_prof_ENS" comment_char="#">
        <columns>value, dbkey, name, path</columns>
        <file path="tools/ribo_prof_pipeline/ENS_db.loc" />
    </table>
    <table name="ribo_prof_igenomes" comment_char="#">
        <columns>value, dbkey, name, path</columns>
        <file path="tools/ribo_prof_pipeline/igenomes.loc" />
    </table>
    <table name="ribo_prof_ublast_db" comment_char="#">
        <columns>value, dbkey, name, path</columns>
        <file path="tools/ribo_prof_pipeline/ublast_db.loc" />
    </table>    
	<table name="ribo_prof_blastp_db" comment_char="#">
        <columns>value, dbkey, name, path</columns>
        <file path="tools/ribo_prof_pipeline/blastp_db.loc" />
    </table>
    <table name="ribo_prof_snpdb" comment_char="#">
    	<columns>value, dbkey, name, path</columns>
    	<file path="tools/ribo_prof_pipeline/snpdb.loc" />
    </table>

Add the paths to the Igenome root folder, blast DBs, SNP DBs and Ensembl SQLite databases you installed under "3 Installing prerequisites: Data dependencies" to the specific .loc files. 
	
	- Igenomes: 			See tools/ribo_prof_pipeline/igenomes.loc for how to configure
	- Ensembl sqlite dbs:		See tools/ribo_prof_pipeline/ENS_db.loc for how to configure
	- BlastP dbs:			See tools/ribo_prof_pipeline/blastp_db.loc for how to configurep
	- UBlast dbs:			See tools/ribo_prof_pipeline/ublast_db.loc for how to configure
	- SNP  dbs:			See tools/ribo_prof_pipeline/snpdb.loc for how to configure	 