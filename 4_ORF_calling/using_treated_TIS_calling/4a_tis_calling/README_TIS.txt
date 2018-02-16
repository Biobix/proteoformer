README TIScalling 

1 General prerequisites
-----------------------
- Perl + v.5.10
- SQLite
- DBI package
- Parallel ForkManager package
- Getopt::Long package
- Storable 'dclone' (standard with perl 5.16)
- Cwd

2 Input parameters & How To
---------------------------

Prerequisites:

 -	Igenome of organism (species and Ensembl annotation version specific) (http://tophat.cbcb.umd.edu/igenomes.html) 

 - 	The Igenome folder structure as detailed in the 1_Mapping.pl file.

 - 	SQLite Ensembl DB with tables gene coord_system exon exon_transcript transcript translation seq_region. 
 	See mysqldump2sqlite3.pl for detailed how to. Or download pre-formatted version from (www.biobix.be/ENS)
  
The input variables of the mapping script are:

./TIScalling_categorised.pl --sqlite SQLite/results_STAR12.db --cores 22 --local_max 1 --R_aTIS 0.01 --min_count_aTIS 5 --R_5 0.05 --min_count_5 10 --R_CDS 0.15 --min_count_CDS 15 --R_3 0.05 --min_count_3 10 --R_no_trans 0.05 --min_count_no_trans 10 --igenomes_root /data/igenomes/Mus_musculus/Ensembl/GRCm38/ --ens_db getcwd/SQLite/ENS_mmu_72.db --dir getwd

GetOptions(
"dir=s"=>\$work_dir,                            # Path to the working directory                                                                 optional argument
"cores=i"=>\$cores,                             # Number of cores to be used                                                                    mandatory argument  # -> (default = number of chromosomes)
"tmp:s" =>\$tmpfolder,                          # Folder where temporary files are stored,                                                      optional  argument (default = $TMP env setting)
"igenomes_root:s" =>\$IGENOMES_ROOT,            # Igenomes root folder                                                                          mandatory argument
"sqlite_db:s" =>\$sqlite_db,                    # SQLite results db with mapping and tr_translation tables                                      mandatory argument
"ens_db:s" =>\$ens_db,                          # SQLite Ensembl DB                                                                             mandatory argument
"local_max:i" =>\$local_max,                    # The range wherein the localmax is (in basepairs, e.g. 1 means +- 1 basepair)                  optional argument (default 1) [*]
"min_count_aTIS:i" =>\$min_count_aTIS,          # The minimum count of riboseq profiles mapping to the aTIS site                                optional argument (default 5) [**]
"R_aTIS:f" => \$R_aTIS,                         # The Rltm - Rchx value calculated based on both CHX and LTM data for a aTIS                    optional argument (default .01) [***]
"min_count_5:i" =>\$min_count_5,                # The minimum count of riboseq profiles mapping to the 5'UTR site                               optional argument (default 10) [**]
"R_5:f" => \$R_5,                               # The Rltm - Rchx value calculated based on both CHX and LTM data for a 5'UTR TIS               optional argument (default .05) [***]
"min_count_CDS:i" =>\$min_count_CDS,            # The minimum count of riboseq profiles mapping to the CDS site                                 optional argument (default 15) [**]
"R_CDS:f" => \$R_CDS,                           # The Rltm - Rchx value calculated based on both CHX and LTM data for a CDS TIS                 optional argument (default .15) [***]
"min_count_3:i" =>\$min_count_3,                # The minimum count of riboseq profiles mapping to the 3'UTR site                               optional argument (default 10) [**]
"R_3:f" => \$R_3,                               # The Rltm - Rchx value calculated based on both CHX and LTM data for a 3'UTR TIS               optional argument (default .05) [***]
"min_count_no_trans:i" =>\$min_count_no_trans,  # The minimum count of riboseq profiles mapping to the no translation site                      optional argument (default 10) [**]
"R_no_trans:f" => \$R_no_trans,                 # The Rltm - Rchx value calculated based on both CHX and LTM data for a no translation TIS      optional argument (default .05) [***]
"out_sqlite:s" =>\$out_sqlite                   # Galaxy specific history file location                                                         Galaxy specific
);

[*]     The min_count is defined as the minimal number of reads needed on a single position to define that position as a peak. 
[**]    The local_max (in basepairs) is defined as the region (position -/+ local_max) around a specific position 
        in which the position has to have the highest number of reads to be able to qualify as a peak position. 
[***] 	R = Rltm-Rchx in which Rk = (Xk/Nk) * 10 (k = LTM, CHX), Xk is the number of reads on that transcript in data k. (Lee et al. 2012)

3 When is a ribosome profile a TIS
-----------------------

- 	In a first step, ribosome profiles are combined into peaks. If the ribosome profile is detected on the first position of a ATG or near-cognate start codon it is a peak. 
	If the +/- 1 position is an ATG or near-cognate, the ribosome profile is also a peak, but the position is changed. As such, different profiles can be combined into 1 peak. 
	If this is the case, the ribosome profile hits are added up. These so called "peak_shifts" are stored in the TIS SQLite DB (see output).   

- 	In a second step for each (combined) peak the three parameters (min_count, local_max and R) are calculated, the combined peaks are only retained as true TISses if they meet the settings of the three already discussed parameters. 

-	The min_count and R parameters can be set specifically for each annotation clas; 5'UTR, annotated TIS (aTIS), coding sequence (CDS), 3'UTR, no translation(no_trans).

4 Output 
---------

SQLite DB (results.db)

    1) Overview table (TIS_overview)
    	For each analysis a new row is added to the overview table. 
        Data for the TIScalling parameters are stored with an attached ID value.
            
            ID  |   local_max   |   min_count_aTIS  |   R_aTis  |   min_count_5UTR  |   R_5UTR  |   min_count_CDS   |   R_CDS   |   min_count_3UTR  |   R_3UTR  |   min_count_no_trans  |   R_no_trans
			---- --------------- ------------------- ----------- ------------------- ----------- ------------------- ----------- ------------------- ----------- ----------------------- -------------
			1   |   1           |   5               |   0.01    |   10              |   0.05    |   15              |   0.15    |   10              |   0.05    |   10                  |   0.05
			2   |   1           |   5               |   0.01    |   10              |   0.05    |   15              |   0.15    |   10              |   0.05    |   10                  |   0.05
			3   |   1           |   5               |   0.01    |   10              |   0.05    |   15              |   0.15    |   10              |   0.05    |   10                  |   0.05
			...
            
            
    2) Analysis specific table (TIS_ID)
    	All TISses identified during the analysis of an experiment (with specific arguments and values for parameters) are stored in a table. 
        Data available in the output for each TIS in the case of Ensembl transcripts
			
			transcript_id   |   stable_id               |   biotype                 |   chr |   strand  |   start       |   dist_to_transcript_start    |   dist_to_aTIS    |   annotation      |   start_codon     |   peak_shift  |   count   |   Rltm_min_Rchx
			---------------- --------------------------- --------------------------- ------- ----------- --------------- ------------------------------- ------------------- ------------------- ------------------- --------------- ----------- --------------------
			356625          |   ENSMUST00000161973      |   processed_pseudogene    |   1   |   1       |   143683266   |   37                          |   NA              |   no_translation  |   AGG             |   +1 0 -1     |   151.0   |   0.822453056367707
			356625          |   ENSMUST00000161973      |   processed_pseudogene    |   1   |   1       |   143683274   |   45                          |   NA              |   no_translation  |   ATT             |   +1 0 -1     |   57.0    |   0.183378085739926
			356625          |   ENSMUST00000161973      |   processed_pseudogene    |   1   |   1       |   143683840   |   611                         |   NA              |   no_translation  |   GTG             |   0           |   23.0    |   0.116375404130359
			...
    
    
    
    
    
    
    
