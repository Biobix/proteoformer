#!/usr/bin/env python
import traceback
import getopt
import sys
import os
import sqlite3 as sqlite

__author__ = 'Steven Verbruggen'

'''

Use the RiboZINB module (author: Elvis Ndah) for determining actively translated transcript isoforms in PROTEOFORMER.


ARGUMENTS
    -w | --work_dir                     The working directory (default: CWD)
    -x | --tmpfolder                    The temporary folder (default: work_dir/tmp)
    -p | --result_db                    The results database (mandatory)
    -m | --mincount                     The minimum reads count for a transcript to be called (default: 5)
    -n | --no_of_samples                The number of iterations when generating a negative set (default: 30)
    -f | --fdr                          The false discovery rate (default: 0.05)
    -s | --default_score                Use the default score threshold (d) or estimate threshold by permutation test (p)
                                            (default: d)
    -v | --cutoff                       The default score threshold (default: 0.1)
    -a | --alpha                        Proportion of noise to generate negative set (default: 1)



EXAMPLE

python RiboZINB.py -p SQLite/results.db

'''


def main():

    # Catch command line with getopt
    try:
        myopts, args = getopt.getopt(sys.argv[1:], "w:x:p:m:n:f:g:s:v:a:",["work_dir=","tmpfolder=","result_db=","mincount=",\
            "no_of_samples=", "fdr=", "fdr_type=","default_score=","cutoff=", "alpha="])
    except getopt.GetoptError as err:
        print err
        sys.exit()

    # Catch arguments
    # o == option
    # a == argument passed to the o
    for o, a in myopts:
        if o in ('-w', '--work_dir'):
            work_dir = a
        if o in ('-p', '--result_db'):
            result_db = a
        if o in ('-x', '--tmpfolder'):
            tmpfolder = a
        if o in ('-m', '--mincount'):
            mincount = int(a)
        if o in ('-n', '--no_of_samples'):
            no_of_samples = int(a)
        if o in ('-f', '--fdr'):
            fdr = float(a)
        if o in ('-s', '--default_score'):
            default_score = a
        if o in ('-v', '--cutoff'):
            cutoff = float(a)
        if o in ('-a', '--alpha'):
            alpha = float(a)


    #Parse arguments
    try:
        work_dir
    except:
        work_dir = os.getcwd()
    os.chdir(work_dir)
    try:
        tmpfolder
    except:
        tmpfolder = work_dir+"/tmp"
    try:
        result_db
    except:
        print "Do not forget the result DB parameter!"
        sys.exit()
    try:
        mincount
    except:
        mincount = 5
    try:
        no_of_samples
    except:
        no_of_samples = 30
    try:
        fdr
    except:
        fdr = 0.05
    try:
        default_score
    except:
        default_score = 'd'
    if default_score!='d' and default_score!='p':
        print "ERROR: default score argument should be 'd' or 'p'!"
        sys.exit()
    default_score_print=''
    default_score_command=''
    if default_score=='d':
        default_score_print='default_score'
        default_score_command='Y'
    else:
        default_score_print='permutation_test'
        default_score_command='N'
    try:
        cutoff
    except:
        cutoff=0.1
    try:
        alpha
    except:
        alpha = int(1)

    # Print input arguments
    print "PROTEOFORMER: RiboZINB"
    print "----------------------"
    print
    print "Input arguments:"
    print
    print "Working directory:                                           "+work_dir
    print "Temporary files folder:                                      "+tmpfolder
    print "Results DB:                                                  "+result_db
    print "Minimum reads count:                                         "+str(mincount)
    print "Number of iterations when generating negative set:           "+str(no_of_samples)
    print "False discovery rate:                                        "+str(fdr)
    print "Default score or permutation test:                           "+default_score_print
    print "Default score threshold:                                     "+str(cutoff)
    print "Noise proportion factor alpha:                               "+str(alpha)
    print
    sys.stdout.flush()

    #Create tmp folder
    if not os.path.isdir(tmpfolder):
        os.system("mkdir " + tmpfolder)
    ribozinb_tmp = tmpfolder+"/RiboZINB"
    if not os.path.isdir(ribozinb_tmp):
        os.system("mkdir "+ribozinb_tmp)

    # Get the arguments out of resultDB
    ensDB, igenomes_root, species, ens_version, cores, exp_name = get_arguments(result_db)

    print "Input arguments from "+result_db+":"
    print
    print "Ensembl DB:                                                  "+ensDB
    print "Ensembl version:                                             "+str(ens_version)
    print "igenomes root folder:                                        "+igenomes_root
    print "Species:                                                     "+species
    print "Number of cores:                                             "+str(cores)
    print "Experiment name:                                             "+exp_name
    print
    sys.stdout.flush()

    # Conversion of species terminology
    speciesLatin = "Mus_musculus" if species == "mouse" else \
        "Homo_sapiens" if species == "human" else \
        "Arabidopsis_thaliana" if species == "arabidopsis" else \
        "Drosophila_melanogaster" if species == "fruitfly" else ""
    speciesShort = "mmu" if species == "mouse" else \
        "hsa" if species == "human" else \
        "ath" if species == "arabidopsis" else \
        "dme" if species == "fruitfly" else ""
    # Assembly
    assembly = "GRCm38" if species == "mouse" and ens_version >= 70 else \
        "NCBIM37" if species == "mouse" and ens_version < 70 else \
        "GRCh38" if species == "human" and ens_version > 75 else \
        "GRCh37" if species == "human" and ens_version <= 75 else \
        "TAIR10" if species == "arabidopsis" else \
        "BDGP5" if species == "fruitfly" else ""

    #Find the correct Ensembl coord system ID
    coord_system_id = find_coord_system_id(ensDB)

    # Get the counts table in csv format
    print "\nGet the counts table in CSV format\n"
    sys.stdout.flush()
    counts_csv = counts_table_export(result_db, "count_fastq1", "CHX", ribozinb_tmp)

    #GTF file
    gtf_file = igenomes_root+"/"+speciesLatin+"/Ensembl/"+assembly+"/Annotation/Genes/genes_"+str(ens_version)+".gtf"

    #Total mappable reads
    mapped_total = get_mapped_total(result_db, "fastq1")

    #Download all RiboZINB scripts from GitHub
    print
    sys.stdout.flush()
    main_script_folder, R_scripts_folder = set_scripts(work_dir)

    #Execute RiboZINB
    command = "perl "+main_script_folder+"/RiboZINB.pl -p "+counts_csv+" -g "+gtf_file+" -e "+exp_name+" -w "\
        +ribozinb_tmp+" -m "+str(mincount)+" -r "+str(mapped_total)+" -d N -t "+str(cores)+" -n "\
        +str(no_of_samples)+" -s "+R_scripts_folder+" -f "+str(fdr)+" -v "+str(cutoff)+" -dt "\
        +default_score_command+" -a "+str(alpha)
    print "\nExecuting RiboZINB"
    print "- - - - - - - - - - - - - "
    print "Command:"
    print "\t"+command
    sys.stdout.flush()
    os.system(command)
    print"- - - - - - - - - - - - - - "
    print
    sys.stdout.flush()
    os.chdir(work_dir)

    #Move S curve to the working directory of proteoformer
    os.system("mv "+ribozinb_tmp+"/"+exp_name+"_scurve.pdf "+work_dir)

    #Parse necessary table into Result DB table and search extra info from ensemblDB
    print "Parse results into results DB\n"
    sys.stdout.flush()
    parse_results(result_db, ensDB, ribozinb_tmp, exp_name, coord_system_id)

    print "Remove RiboZINB scripts\n"
    sys.stdout.flush()
    os.system("rm -rf "+main_script_folder)

    #print "Remove RiboZINB tmp folder\n"
    os.system("rm -rf "+ribozinb_tmp)

    print "--- DONE ---\n\n"
    sys.stdout.flush()

##########
#  SUBS  #
##########

### Parse results into a transcripts table
def parse_results(result_db, ensDB, ribozinb_tmp, exp_name, coord_system_id):

    #Init
    expressed_isoforms_file = ribozinb_tmp+"/"+exp_name+"_Ribo_expressed_isoforms.txt"
    transcript_i=1

    #Remove header
    without_header = ribozinb_tmp+"/expressed_isoforms_without_header.txt"
    os.system("sed '1d' "+expressed_isoforms_file+" > "+without_header)

    #Establish Ensembl DB connection
    try:
        con_ens = sqlite.connect(ensDB)
    except:
        print "ERROR: Could not connect to "+ensDB
        sys.exit()
    cur_ens = con_ens.cursor()

    #Define csv file
    csv_path = ribozinb_tmp + "/features_tr_" + exp_name + ".csv"
    with open(csv_path, 'w') as FW:
        #Read in results table
        with open(expressed_isoforms_file, 'r') as FR:
            next(FR)
            for line in FR:
                line = line.rstrip() #Remove trailing whitespace

                #Convert the line to a features list needed for the transcript table in the results DB
                features, above_cutoff = convert_to_results_db_list(line, cur_ens, transcript_i, coord_system_id)

                #Convert the features to a line for the csv file
                line_csv = convert_features_to_csv_line(features)

                #Write to csv if the isoform passed the score cutoff
                if above_cutoff == 'Y':
                    FW.write(line_csv)

                #Loop number aug
                transcript_i += 1

    #Clean up
    os.system("rm -rf "+without_header)

    #Dump csv file in SQLite
    dump_csv(csv_path, result_db)

    return

### Dump csv file into SQLite
def dump_csv(csv_path, result_db):

    #Init
    table_name = "tr_translation_ribozinb"

    #Open results DB connection
    try:
        con = sqlite.connect(result_db)
    except:
        print "ERROR: could not connect to "+result_db+"\n"
        sys.exit()
    cur = con.cursor()

    #Drop existing table
    query = "DROP TABLE IF EXISTS "+table_name
    cur.execute(query)

    #Create new table
    query = "CREATE TABLE IF NOT EXISTS "+table_name+" ("\
                "transcript_id VARCHAR(100) NOT NULL,"\
                "stable_id VARCHAR(100) NOT NULL,"\
                "chr char(50) NOT NULL default '',"\
                "seq_region_id VARCHAR(10) NOT NULL,"\
                "seq_region_strand VARCHAR(2) NOT NULL,"\
                "seq_region_start FLOAT NOT NULL,"\
                "seq_region_end FLOAT NOT NULL,"\
                "read_counts FLOAT NOT NULL,"\
                "normalized_counts FLOAT NOT NULL,"\
                "biotype VARCHAR(100) NOT NULL,"\
                "exon_coverage VARCHAR(5) NOT NULL,"\
                "canonical VARCHAR(5) NOT NULL,"\
                "ccds VARCHAR(20) NOT NULL,"\
                "gene_stable_id VARCHAR(100) NOT NULL,"\
                "PRIMARY KEY (stable_id, gene_stable_id)"\
                ");"
    cur.execute(query)

    #Create indices
    #index1 is on exon_coverage but not needed for ribozinb
    index2 = table_name+"_fastq1_seq_region_id"
    query = "CREATE INDEX IF NOT EXISTS "+index2+" ON "+table_name+" (seq_region_id);"
    cur.execute(query)

    #Dump csv
    os.system("sqlite3 -separator , "+result_db+" \".import "+csv_path+" "+table_name+" \"")
    #Remove tmp csv file
    os.system("rm -rf "+csv_path)

    #Add RiboZINB argument to arguments table
    query = "DELETE FROM arguments WHERE value='tr_calling';"
    cur.execute(query)
    query = "INSERT INTO arguments (variable, value) VALUES (\"tr_calling\", \"ribozinb\");"
    cur.execute(query)

    return

### Convert features to line for csv file
def convert_features_to_csv_line(features):

    #Init
    line=""

    for feature in features:
        line = line+str(feature)+","
    line = line.rstrip(",")
    line += "\n"

    return line

### Convert the features in the expressed isoforms file to the features needed for the transcript table
def convert_to_results_db_list(line, cur_ens, transcript_id, coord_system_id):

    #Init desired features list for transcripts table
    features_tr=[]

    #Split the line into the different features (list)
    features_expressed_iso = line.split("\t")
    above_cutoff = features_expressed_iso[24]

    #transcript id
    features_tr.append(search_ensembl_transcript_id(features_expressed_iso[2]), cur_ens)
    #transcript stable id
    features_tr.append(features_expressed_iso[2])
    #Chr and seq_region_id
    chr, seq_region_id = search_chr_of_transcript(features_tr[1], cur_ens, coord_system_id)
    features_tr.append(chr)
    features_tr.append(seq_region_id)
    #seq_region_strand
    if features_expressed_iso[3]=="+":
        seq_region_strand="1"
    elif features_expressed_iso[3]=="-":
        seq_region_strand="-1"
    else:
        print "Parsing ERROR for "+features_tr[1]+". Strand is not + or -!"
    features_tr.append(seq_region_strand)
    #seq_region_start and seq_region_end
    seq_region_start, seq_region_end = search_start_and_end(features_tr[1], cur_ens)
    features_tr.append(seq_region_start)
    features_tr.append(seq_region_end)
    #Read counts
    features_tr.append(int(features_expressed_iso[13]))
    #Normalized counts
    features_tr.append(float(features_tr[7]/float(features_expressed_iso[8])))
    #Biotype
    features_tr.append(features_expressed_iso[6])
    #Exon coverage
    features_tr.append("NA")
    #Canonical
    gene_stable_id = features_expressed_iso[0]
    features_tr.append(check_if_canonical(features_tr[1], gene_stable_id, cur_ens))
    #CCDS
    features_tr.append(check_CCDS(features_expressed_iso[4], features_expressed_iso[5]))
    # Gene stable ID
    features_tr.append(gene_stable_id)

    return features_tr, above_cutoff

### Search Ensembl transcript ID based on transcript stable ID
def search_ensembl_transcript_id(stable_id, cur_ens):
    
    #Init
    transcript_id=0;
    
    #Search in Ensembl db
    query = "SELECT transcript_id FROM transcript WHERE stable_id='"+stable_id+"';"
    if cur_ens.execute(query):
        transcript_id = int(cur_ens.fetchone()[0])
    else:
        print "ERROR: could not find the Ensembl transcript ID based on the stable ID for transcript "+stable_id
        print "Query: "+query
    
    return transcript_id

### Check if the transcript is the CCDS and parse CCDS number
def check_CCDS(boolean_CCDS, CCDS_number):

    #Init
    CCDS_tr = ""

    if boolean_CCDS=="Y":
        CCDS_tr = CCDS_number
    elif boolean_CCDS=="N":
        CCDS_tr="No"

    return CCDS_tr

### Check if the transcript is the canonical transcript of the gene
def check_if_canonical(tr_stable_id, gene_stable_id, cur_ens):

    #Init
    canonical=""
    canonical_transcript_id=""

    #Get the canonical transcript ID of the gene
    query = "SELECT tr.stable_id FROM transcript AS tr JOIN gene AS g ON tr.transcript_id=g.canonical_transcript_id" \
            " WHERE g.stable_id=\""+gene_stable_id+"\";"
    if cur_ens.execute(query):
        canonical_transcript_id = cur_ens.fetchone()[0]
    else:
        print "ERROR: Could not find the canonical transcript ID of gene "+gene_stable_id
        print "QUERY: "+query

    #Compare the given transcript ID with the canonical transcript ID of the gene
    if tr_stable_id==canonical_transcript_id:
        canonical = "Yes"
    else:
        canonical = "No"

    return canonical

### Search seq_region_start and seq_region_end of the transcript
def search_start_and_end(tr_stable_id, cur_ens):

    #Init
    start=0
    end=0

    query = "SELECT seq_region_start, seq_region_end FROM transcript WHERE stable_id=\""+tr_stable_id+"\";"
    if cur_ens.execute(query):
        output = cur_ens.fetchone()
        start = int(output[0])
        end = int(output[1])
    else:
        print "ERROR: Could not find start and end of transcript "+tr_stable_id

    return start, end

### Search chr of transcript ID
def search_chr_of_transcript(tr_stable_id, cur_ens, coord_system_id):

    #Init
    seq_region_id=0
    chrm=''

    #Find seq region id of the transcript
    query = "SELECT seq_region_id FROM transcript WHERE stable_id=\""+tr_stable_id+"\";"
    if cur_ens.execute(query):
        seq_region_id = cur_ens.fetchone()[0]
    else:
        print "Could not find the seq region ID of transcript "+tr_stable_id
        print "QUERY: "+query

    #Convert seq region id to chromosome name
    query = "SELECT name FROM seq_region WHERE coord_system_id="+str(coord_system_id)+" AND seq_region_id="+str(seq_region_id)+";"
    if cur_ens.execute(query):
        chrm = cur_ens.fetchone()[0]
    else:
        print "Could not find the chromosome name of seq region ID "+seq_region_id
        print "QUERY: "+query

    return chrm, seq_region_id

### Download all RiboZINB scripts from GitHub
def set_scripts(wd):

    #Init
    main_scripts_folder = wd+"/RiboZINB_scripts"
    R_scripts_folder = wd+"/RiboZINB_scripts/Rscipts"

   #Remove previous versions
    os.system("rm -rf "+main_scripts_folder)

    #Download from GitHub and unpack
    os.system("wget -q --no-check-certificate https://github.com/Biobix/RiboZINB/archive/master.zip")
    os.system("unzip -q master.zip")
    os.system("mkdir "+main_scripts_folder)
    os.system("mv RiboZINB-master/* "+main_scripts_folder)
    os.system("rm -rf master.zip")
    os.system("rm -rf RiboZINB-master")

    #Check if all scripts are present
    if (os.path.isfile(main_scripts_folder+"/RiboZINB.pl") and os.path.isfile(R_scripts_folder+"/FDR.R") and \
        os.path.isfile(R_scripts_folder+"/RiboZINB.R") and os.path.isfile(R_scripts_folder+"/merge.R") and \
        os.path.isfile(R_scripts_folder+"/s_curve.R")):
        print "RiboZINB scripts downloaded\n"
    else:
        print "Could not found all necessary scripts\n"
        sys.exit()

    return main_scripts_folder, R_scripts_folder

### Get total mapped reads
def get_mapped_total(db, treat):
    try:
        con = sqlite.connect(db)
    except:
        print "Could not connect to "+db
        sys.exit()

    #Init
    mapped_total = 0

    with con:
        cur = con.cursor()
        if cur.execute("SELECT mapped_T FROM statistics WHERE sample LIKE '%"+treat+"%' AND type='genomic';"):
            mapped_total = int(cur.fetchone()[0])
        else:
            print "Could not fetch the total mapped reads in "+db

    return mapped_total


### Export the counts (CHX) table to CSV format
def counts_table_export(result_db, table_name, treat, prepath):
    output_name = prepath+"/counts_"+treat+".csv"
    command = "sqlite3 -csv "+result_db+" \"SELECT * FROM "+table_name+";\" > "+output_name
    os.system(command)
    return output_name

def find_coord_system_id(ensDB):

    #Connect to Ensembl DB
    try:
        con = sqlite.connect(ensDB)
    except:
        print "ERROR: Could not connect to "+ensDB
        sys.exit()

    cur = con.cursor()
    query = "SELECT coord_system_id FROM coord_system WHERE name=\"chromosome\" AND rank=1;"
    cur.execute(query)
    coord_system_id = cur.fetchone()[0]

    return coord_system_id

### Get arguments out of resultsDB
def get_arguments(db):
    try:
        con = sqlite.connect(db)
    except:
        print "Could not connect to "+db
        sys.exit()

    #Init
    ens_db=''
    igenomes_root=''
    species=''
    ens_v=''
    cores=''
    run_name=''

    with con:
        cur = con.cursor()

        if cur.execute("SELECT value FROM arguments WHERE variable='ens_db';"):
            ens_db = str(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the Ensembl database argument in "+db
            sys.exit()
        if cur.execute("SELECT value FROM arguments WHERE variable='igenomes_root';"):
            igenomes_root = str(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the igenomes root argument in "+db
            sys.exit()
        if cur.execute("SELECT value FROM arguments WHERE variable='species';"):
            species = str(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the species argument in " + db
            sys.exit()
        if cur.execute("SELECT value FROM arguments WHERE variable='ensembl_version';"):
            ens_v = int(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the ensembl version argument in " + db
            sys.exit()
        if cur.execute("SELECT value FROM arguments WHERE variable='nr_of_cores';"):
            cores = int(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the cores argument in " + db
            sys.exit()
        if cur.execute("SELECT value FROM arguments WHERE variable='run_name';"):
            run_name = str(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the run name argument in " + db
            sys.exit()

    return ens_db, igenomes_root, species, ens_v, cores, run_name

#########Set MAIN####################
if __name__ == "__main__":
    try:
        main()
    except Exception, e:
        traceback.print_exc()
#####################################