Provides galaxy tools for the PROTEOFORMER pipeline -  http://www.biobix.be/PROTEOFORMER
The PROTEOFORMER pipeline cannot be installed automatically from a toolshed

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


########################################################
## MANUAL INSTALLATION FOR THE PROTEOFORMER PIPELINE: ##
########################################################

1  Galaxy specific configuration:
---------------------------------
To manually install the PROTEOFORMER pipeline a set of configuration files, wrapper scripts, datatypes and other configuration and 
location files need to be added to specific folders of your galaxy system. These are listed below and can be found within this directory.

	- Add configuration files (.xml), perl code (.pl), shell code (.sh) and location files (.loc) from tools/proteoformer to your galaxy installation
		
		Overview:
		---------
			- tools/proteoformer/1_mapping.pl
			- tools/proteoformer/1_mapping.xml
			- tools/proteoformer/2_assembly.pl
			- tools/proteoformer/2_assembly.xml
			- tools/proteoformer/blastp_db.loc[*]
			- tools/proteoformer/ENS_db.loc[*]
			- tools/proteoformer/filterSAMfile.pl
			- tools/proteoformer/gene_distribution.pl
			- tools/proteoformer/gene_distribution.xml
			- tools/proteoformer/generate_translation_db.pl
			- tools/proteoformer/generate_translation_db.xml
			- tools/proteoformer/igenomes.loc[*]
			- tools/proteoformer/metagenic_classification.pl
			- tools/proteoformer/metagenic_classification.xml
			- tools/proteoformer/metagenenic_piecharts.R
			- tools/proteoformer/quality.R
			- tools/proteoformer/ribo_translation.pl
			- tools/proteoformer/ribo_translation.xml
			- tools/proteoformer/snp_calling
			- tools/proteoformer/snp_calling.xml
			- tools/proteoformer/snpdb.loc[*]
			- tools/proteoformer/snpIndexBuilder.pl
			- tools/proteoformer/splitVCFaltRecords.pl
			- tools/proteoformer/tis_calling.xml
			- tools/proteoformer/TIS_overview.pl
			- tools/proteoformer/tis_overview.xml
			- tools/proteoformer/Tiscalling_categorised.pl
			- tools/proteoformer/ublast_db.loc[*]
			
		[*] See point 4: "Tool data location files" for more information on how to use the location (.loc) files
		
	- Add the datatype definition file from lib/galaxy/datatypes/proteoformer.py to your galaxy installation
	
	- Add the following import line to:  lib/galaxy/datatypes/registry.py
			import proteoformer # added for PROTEOFORMER pipeline
			
	- Add datatypes between the <registration>   </registration> tags in:  datatypes_conf.xml
		 	
		 <!-- Start PROTEOFORMER Datatypes -->
    			<datatype extension="tis" type="galaxy.datatypes.proteoformer:TisOverview"/>
   			<datatype extension="sqlitedb" type="galaxy.datatypes.proteoformer:SqliteDb" />
    		<!-- End PROTEOFORMER Datatypes -->
    		
    - Copy the proteoformer_tool_conf.xml file to your galaxy installation and add it to the universe_wsgi.ini file:
    		
    		# Tool config files, defines what tools are available in Galaxy.
		# Tools can be locally developed or installed from Galaxy tool sheds.
		tool_config_file = proteoformer_tool_conf.xml,tool_conf.xml
		
	  or add the <section name="PROTEOFORMER pipeline" id=“proteof”ormer>…</section> from the proteoformer_tool_conf to your tool_conf.xml
	  
	- Reorganize the integrated_tool_panel.xml and add the following lines where you want them to appear in the galaxy tool column.
		
    		<section id=“proteoformer” name="PROTEOFORMER pipeline" version="">
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

The PROTEOFORMER pipeline is primarily written in Perl code, using different Perl specific packages. 
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
		- STAR, v2.4.2a or higher (https://code.google.com/p/rna-star/)
		- TOPHAT2, v2.0.13 or higher (http://tophat.cbcb.umd.edu/)
		- BLASTP ( ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) or USEARCH[**] (http://www.drive5.com/usearch/download.html)
		- R (http://www.r-project.org/)
		- SAMTOOLS, v. 1.19 or higher (http://sourceforge.net/projects/samtools/files/samtools/)
		- GATK (http://www.broadinstitute.org/gatk/download)
		- PICARD (http://sourceforge.net/projects/picard/files/picard-tools/)
		- SQLITE3 (http://www.sqlite.org/download.html)
		- FASTX toolkit, v.0.0.13 or higher (http://hannonlab.cshl.edu/fastx_toolkit/)
		
		[**] The usearchX.Y.Z executatble should be renamed to "usearch"
		
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
								IGENOMES_ROOT/Mus_musculus/Ensembl/GRCm38/Sequence/STARIndex (this folder and the STAR indexes are created 
								automatically when ther first STAR job is launched)
                - rRNA seq file	(custom made)                   IGENOMES_ROOT/Mus_musculus/Ensembl/GRCm38/Sequence/AbundantSequences/species_rRNA.fa[***]
		- rRNA seq file	(custom made)			IGENOMES_ROOT/Mus_musculus/Ensembl/GRCm38/Sequence/AbundantSequences/species_tRNA.fa[***]
		- rRNA seq file	(custom made)			IGENOMES_ROOT/Mus_musculus/Ensembl/GRCm38/Sequence/AbundantSequences/species_sn-o-RNA.fa[***]
		
        [***] use 3-letter abbreviation for species: e.g. dme for Drosophila melanogaster, hsa for Homo sapiens, mmu for Mus musculus
              rRNA, tRNA, sn-o-RNA sequence fasta files have to be compiled. BioMart Ensembl allows you to download specific types of RNA sequences.

											
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
 		  		
 		  		Examples: 				ENS_hsa_75.db
 		  							ENS_mmu_75.db
 		 							ENS_dme_75.db
 		  							ENS_ath_16.db
 		  					
 		- You can also download pre-formatted sqlite databases from the website (www.biobix.be/PROTEOFORMER)
														
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

    <!-- Location of PROTEOFORMER tool files  -->
    <table name="PROTEOFORMER_ENS" comment_char="#">
        <columns>value, dbkey, name, path</columns>
        <file path="tools/proteoformer/ENS_db.loc" />
    </table>
    <table name="PROTEOFORMER_igenomes" comment_char="#">
        <columns>value, dbkey, name, path</columns>
        <file path="tools/proteoformer/igenomes.loc" />
    </table>
    <table name="PROTEOFORMER_BlastDB" comment_char="#">
        <columns>value, dbkey, name, path</columns>
        <file path="tools/proteoformer/blast_db.loc" />
    </table>
    <table name="PROTEOFORMER_ublast_db" comment_char="#">
        <columns>value, dbkey, name, path</columns>
        <file path="tools/proteoformer/ublast_db.loc" />
    </table>
    <table name="PROTEOFORMER_blastp_db" comment_char="#">
        <columns>value, dbkey, name, path</columns>
        <file path="tools/proteoformer/blastp_db.loc" />
    </table>
    <table name="PROTEOFORMER_snpdb" comment_char="#">
    	<columns>value, dbkey, name, path</columns>
    	<file path="tools/proteoformer/snpdb.loc" />
    </table>

Add the paths to the Igenome root folder, blast DBs, SNP DBs and Ensembl SQLite databases you installed under "3 Installing prerequisites: Data dependencies" to the specific .loc files. 
	
	- Igenomes: 			See tools/proteoformer/igenomes.loc for how to configure
	- Ensembl sqlite dbs:		See tools/proteoformer/ENS_db.loc for how to configure
	- BlastP dbs:			See tools/proteoformer/blastp_db.loc for how to configurep
	- UBlast dbs:			See tools/proteoformer/ublast_db.loc for how to configure
	- SNP  dbs:			See tools/proteoformer/snpdb.loc for how to configure	 