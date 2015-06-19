README sORF assembly intergenic (JeroenC)



TODOTODTODTODTOTDODODTODOTODTODTODTODTODTODTODTODOTDOTDOTO




1 General prerequisites
-----------------------
- Perl + v.5.10
- SQLite
- DBI package
- Parallel ForkManager package
- Getopt::Long package
- Storable 'dclone' (standard with perl 5.16)

2 Input parameters & How To
---------------------------

Prerequisites:

 - Igenome of organism (http://tophat.cbcb.umd.edu/igenomes.html) 

 - The Igenome folder structure as detailed in the Mapping.pl file.

 - SQLite Ensembl DB with tables gene coord_system exon exon_transcript transcript translation seq_region. See mysqldump2sqlite3.pl for detailed how to.
  
The input variables of the mapping script are:

./TIScalling_v3.pl --name mESC_GA_STAR_N --CHX lane1 --LTM lane2 --species mouse --dir /data/RIBO_runs/RIBO_Ingolia_GerbenM/ --min_count 10 --local_max 1 --R 0.05 --cores 22 


GetOptions(

"name=s"=>\$run_name, 				# Name of the run 						mandatory argument

"CHX=s"=>\$CHX, 				# Name of CHX reads lane 					mandatory argument

"LTM=s"=>\$LTM, 				# Name of LTM reads lane 					mandatory argument

"species=s"=>\$species, 			# Species name, eg mouse/human/(fruitfly)			mandatory argument

"dir=s"=>\$work_dir,            		# Path to the working directory                                 mandatory argument

"min_count=s"=>\$min_count, 			# Value for min count parameter 				optional argument (default = 10 => Lee et al, 2012)     [*]

"local_max=s"=>\$local_max,			# Value for local max parameter 				optional argument (default = 1 => Lee et al, 2012)      [**]

"R=s"=>\$R, 					# Value for R							optional argument (default = 0.05 => Lee et al, 2012)   [***]

"cores=s"=>\$cores, 				# Number of cores to be used					optional argument  # Check if possible -> (default = number of chromosomes)

"tmp:s" =>\$tmpfolder         			# Folder where temporary files are stored,                      optional  argument (default = $TMP env setting)

);

[*]     The min_count is defined as the minimal number of reads needed on a single position to define that position as a peak. 
[**]    The local_max (in codons) is defined as the region (position -/+ local_max) around a specific position 
        in which the position has to have the highest number of reads to be able to qualify as a peak position. 
[***] 	R = Rltm-Rchx in which Rk = (Xk/Nk) * 10 (k = LTM, CHX), Xk is the number of reads on that transcript in data k. (Lee et al. 2012)

3 When is a ribosome profile a TIS
-----------------------

- In a first step, ribosome profiles are combined into peaks. If the ribosome profile is detected on the first position of a ATG or near-cognate start codon it is a peak. If the +/- 1 position is an ATG near-cognate, the ribosome profile is also a peak and the position changes. As such, different profiles can be combined into 1 peak. If this is the case, the ribosome profile hits are added up. These so called "peak_shifts" are stored in the TIS SQLite DB(see output).   

- In a second step for (each combined) peak the three parameters (min_count, local_max and R) are calculated, the combined peaks are only retained as true TISses if they meet the three parameters. 


4 Output 
---------

SQLite DB (run_name_results.db)

    1) Overview table (run_name_TIS_overview)
    For each analysis a new row is added to the overview table. 
        Data for some of the parameters and arguments are stored with an attached ID value.
            - ID
            - CHX lane
            - LTM lane
            - min_count
            - local_max
            - R
    2) Analysis specific table (run_name_TIS_ID)
    All TISses identified during the analysis of an experiment (with specific arguments and values for parameters) are stored in a table. 
        Data available in the output for each TIS in the case of Ensembl transcripts
            - Transcript ID
            - Transcript Stable ID
            - Transcript biotype
            - Chromosome
            - Strand
            - TIS position
            - Distance between TIS position and transcript start site
            - Distance between TIS position and annotated TIS
            - Annotation of TIS (5'UTR, aTIS, CDS, 3' UTR / NA in the case of ncRNA transcript)
            - Start codon for TIS
            - Peak shift (Information on ribosome profiles used for the construction of the combined peak)
            - Number of ribosome profile hits on TIS
            - Value for Rltm-Rchx 
    
    
    
    
    
    
    
