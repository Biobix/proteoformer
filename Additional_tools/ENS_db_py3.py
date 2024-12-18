__author__ = 'vladie'

##############################################
## Ensembl sqlite3 database creation script ##
##############################################

###########
## USAGE ##
###########

'''
ARGUMENTS:

-v | --version                              Ensembl annotation version to download (supported versions: from 74)
-s | --species                              Specify the desired species for which gene annotation files should be downloaded
                                            currently supported species:

                                            human                       |   homo_sapiens
                                            mouse                       |   mus_musculus
                                            horse						|   equus_caballus
                                            chinese_hamster             |   Cricetulus griseus
                                            arctic_squirrel             |   urocitellus_parryii
                                            fruitfly                    |   drosophila_melanogaster
                                            saccharomyces_cerevisiae    |   yeast
                                            caenorhabditis_elegans      |   c.elegans
                                            SL1344                      |   Salmonella enterica subsp. enterica serovar Typhimurium str. SL1344
                                            MYC_ABS_ATCC_19977          |   Mycobacterium abscessus atcc 19977
                                            CNECNA3                     |   Cryptococcus_neoformans_var_grubii_h99_gca_000149245
                                            earthmoss                   |   physcomitrium_patens

-h | --help                                 Print this useful help message


The Ensembl database will be installed in the current directory ( directory where the script is run from )

EXAMPLE:

python ENS_db.py -v 98 -s human

DEPENDENCIES:

This program depends upon wget and gzip, which are normally pre-installed on any unix system

'''
import os
import shutil
import sqlite3
import csv
import getopt
import sys
import stat
import traceback

def main():

    try:
        myopts, args = getopt.getopt(sys.argv[1:],"s:v:h",["version=","species=", "help"])
    except getopt.GetoptError as err:
        print(err)
        sys.exit()

    ###############################
    # o == option
    # a == argument passed to the o
    ###############################

    #
    # Catch arguments
    #

    for o, a in myopts:
        if o in ('-h', '--help'):
            print_help()
            sys.exit()
        if o in ('-s','--species'):
            species=a
        if o in ('-v','--version'):
            ens_v=a

    #
    # Check for correct argument, output argument and parse
    #

    if(species == ''):
        print("Error: do not forget to pass the species argument !")
        sys.exit()
    elif(ens_v == ''):
        print("Error: do not forget to pass the ensembl version argument")
        sys.exit()
    print(("Ensembl version used : " + ens_v))
    print(("Selected species     : " + species))

    directory = os.getcwd()
    os.chdir(directory)


    os.system('mkdir tmp')

    if os.path.isdir("tmp/ENS"):
        shutil.rmtree('tmp/ENS')

    os.system('mkdir tmp/ENS')

    #
    # Dictionary for 3letter abbreviation for species
    #

    speciesdict = {'zebrafish': 'dre','danio_rerio': 'dre','human': 'hsa', 'mouse': 'mmu', 'horse': 'eca', 'chinese_hamster': 'cgr', 'fruitfly': 'dme',
                   'yeast': 'sce', 'c.elegans': 'cel',
                   'homo_sapiens': 'hsa', 'mus_musculus': 'mmu', 'equus_caballus': 'eca', 'arctic_squirrel': 'upa', 'drosophila_melanogaster': 'dme',
                   'saccharomyces_cerevisiae': 'sce', 'caenorhabditis_elegans': 'cel','rat': 'rnv',
                   'rattus_norvegicus': 'rnv',
                   'MYC_ABS_ATCC_19977': 'MYC_ABS_ATCC_19977', 'mycobacterium_abscessus_atcc_19977': 'MYC_ABS_ATCC_19977',
                   'SL1344': 'SL1344', 'CNECNA3' : 'CNECNA3', 'Cryptococcus_neoformans_var_grubii_h99_gca_000149245' : 'CNECNA3', 'earthmoss': 'ppa', 'physcomitrium_patens': 'ppa'};


    #
    # Check whether database already exists
    #

    if os.path.isfile(directory + '/ENS_' + speciesdict[species] + '_' + ens_v + '.db'):
        print(("There already exists a file called " + 'ENS_' + speciesdict[species] + '_' + ens_v + '.db in: ' + directory))
        print("Please delete the existing file if you want to create the database anew.")
        sys.exit()

    # create temporary folder when necessary

    # check whether ENS folder exist within tmp folder
    # if true, delete folder and make anew
    # else create ENS folder

    #
    # Function to download necessary files from ftp directory
    #

    def download(ftp_link,core):
        print("Downloading Ensembl database files")
        try:
            os.system('wget -q ' + ftp_link + core +' -P tmp/ENS')
            os.system('wget -q ' + ftp_link + '/coord_system.txt.gz -P tmp/ENS')
            os.system('wget -q ' + ftp_link + '/exon.txt.gz -P tmp/ENS')
            os.system('wget -q ' + ftp_link + '/exon_transcript.txt.gz -P tmp/ENS')
            os.system('wget -q ' + ftp_link + '/gene.txt.gz -P tmp/ENS')
            os.system('wget -q ' + ftp_link + '/seq_region.txt.gz -P tmp/ENS')
            os.system('wget -q ' + ftp_link + '/transcript.txt.gz -P tmp/ENS')
            os.system('wget -q ' + ftp_link + '/translation.txt.gz -P tmp/ENS')
            os.system('wget -q ' + ftp_link + '/xref.txt.gz -P tmp/ENS')
            os.system('wget -q ' + ftp_link + '/object_xref.txt.gz -P tmp/ENS')

        except OSError as e:
            if e.errno == os.errno.ENOENT:
                print(" wget was not found on your system, please install wget in order to run this program")
            else:
                print("something went wrong trying to run wget")
                raise
    #
    # Getting right ftp directory based on arguments
    #

    canEns_v=str(int(ens_v)+53) #Plant and Bacteria Ensembl releases are 53 less than the other species.
    if (species=='human' or species=='homo_sapiens'):
        if(int(ens_v) >= 76 and int(ens_v) <= 113):
             core='/homo_sapiens_core_' + ens_v + '_38.sql.gz'
             download('ftp://ftp.ensembl.org/pub/release-'+ ens_v +'/mysql/homo_sapiens_core_' + ens_v + '_38/',core)
        elif(int(ens_v) >= 74 and int(ens_v) <= 75):
             core='/homo_sapiens_core_' + ens_v + '_37.sql.gz'
             download('ftp://ftp.ensembl.org/pub/release-'+ ens_v +'/mysql/homo_sapiens_core_' + ens_v + '_37/',core)
        else:
            print(("ERROR: unsupported ensembl version: " + ens_v))
            print("supported ensembl versions: human from 74 till 113")
            sys.exit()
    elif (species=='rat' or species=='rattus_norvegicus'):
        if(int(ens_v) >= 74 and int(ens_v) <= 104):
            core='/rattus_norvegicus_core_'+str(ens_v)+'_6.sql.gz'
            download('ftp://ftp.ensembl.org/pub/release-'+ ens_v + '/mysql/rattus_norvegicus_core_' + ens_v +'_6/',core)
        else:
            print(("ERROR: unsupported ensembl version: " + ens_v))
            print("supported ensembl versions: from 75 till 104")
            sys.exit()
    elif (species=='mouse' or species=='mus_musculus'):
        if(int(ens_v) >= 74 and int(ens_v) < 103):
            core='/mus_musculus_core_' + ens_v + '_38.sql.gz'
            download('ftp://ftp.ensembl.org/pub/release-'+ ens_v + '/mysql/mus_musculus_core_' + ens_v +'_38/',core)
        elif(int(ens_v) >= 103 and int(ens_v) <= 113):
            core='/mus_musculus_core_' + ens_v + '_39.sql.gz'
            download('ftp://ftp.ensembl.org/pub/release-'+ ens_v + '/mysql/mus_musculus_core_' + ens_v +'_39/',core)
        else:
            print(("ERROR: unsupported ensembl version: " + ens_v))
            print("supported ensembl versions: from 75 till 113")
            sys.exit()
    elif (species=='horse' or species=='equus_caballus'):
        if(int(ens_v) >= 74 and int(ens_v) <= 100):
            core='/equus_caballus_core_' + ens_v + '_3.sql.gz'
            download('ftp://ftp.ensembl.org/pub/release-'+ ens_v + '/mysql/equus_caballus_core_' + ens_v +'_3/',core)
        else:
            print(("ERROR: unsupported ensembl version: " + ens_v))
            print("supported ensembl versions: from 75 till 100")
            sys.exit()
    elif (species=='chinese_hamster'):
        if(int(ens_v) >=96 and int(ens_v)<=109):
            core='/cricetulus_griseus_picr_core_'+ens_v+'_1.sql.gz'
            download('ftp://ftp.ensembl.org/pub/release-'+ens_v+'/mysql/cricetulus_griseus_picr_core_'+ens_v+'_1/', core)
        else:
            print("ERROR: unsupported ensembl version: "+ens_v)
            print("Supported ensembl versions: from 96 to 109")
            sys.exit()
    elif (species=='arctic_squirrel'):
        if(int(ens_v) >=96 and int(ens_v)<=103):
            core='/urocitellus_parryii_core_'+ens_v+'_1.sql.gz'
            download('ftp://ftp.ensembl.org/pub/release-'+ens_v+'/mysql/urocitellus_parryii_core_'+ens_v+'_1/', core)
        else:
            print(("ERROR: unsupported ensembl version: " + ens_v))
            print("supported ensembl versions: from 96 till 103")
            sys.exit()
    elif (species=='fruitfly' or species=='drosophila_melanogaster'):
        if(int(ens_v) >= 110):
            core='/drosophila_melanogaster_core_' + ens_v + '_10.sql.gz'
            download('ftp://ftp.ensembl.org/pub/release-' + ens_v +'/mysql/drosophila_melanogaster_core_' + ens_v +'_10/',core)
        elif(int(ens_v) >= 109 and int(ens_v) <= 103):
            core='/drosophila_melanogaster_core_' + ens_v + '_9.sql.gz'
            download('ftp://ftp.ensembl.org/pub/release-' + ens_v +'/mysql/drosophila_melanogaster_core_' + ens_v +'_9/',core)
        elif(int(ens_v) >= 102 and int(ens_v) <= 99):
            core='/drosophila_melanogaster_core_' + ens_v + '_8.sql.gz'
            download('ftp://ftp.ensembl.org/pub/release-' + ens_v +'/mysql/drosophila_melanogaster_core_' + ens_v +'_8/',core)
        elif(int(ens_v) >= 98 and int(ens_v) <= 96):
            core='/drosophila_melanogaster_core_' + ens_v + '_7.sql.gz'
            download('ftp://ftp.ensembl.org/pub/release-' + ens_v +'/mysql/drosophila_melanogaster_core_' + ens_v +'_7/',core)
        elif(int(ens_v) >= 74 and int(ens_v) <= 95):
            core='/drosophila_melanogaster_core_' + ens_v + '_6.sql.gz'
            download('ftp://ftp.ensembl.org/pub/release-' + ens_v +'/mysql/drosophila_melanogaster_core_' + ens_v +'_6/',core)
        else:
            print(("ERROR: unsupported ensembl version: " + ens_v))
            print("supported ensembl versions: from 74 till 88 ")
            sys.exit()
    elif (species == 'saccharomyces_cerevisiae' or species == 'yeast'):
        if(int(ens_v) >= 74 and int(ens_v) <= 113):
            core='/saccharomyces_cerevisiae_core_' + ens_v + '_4.sql.gz'
            download('ftp://ftp.ensembl.org/pub/release-' + ens_v +'/mysql/saccharomyces_cerevisiae_core_' + ens_v +'_4/',core)
        else:
            print(("ERROR: unsupported ensembl version: " + ens_v))
            print("supported ensembl versions: from 74 till 113")
            sys.exit()
    elif (species=='caenorhabditis_elegans' or species =="c.elegans"):
        if (int(ens_v) >= 91 and int(ens_v) <= 92):
            core = '/caenorhabditis_elegans_core_' + ens_v + '_260.sql.gz'
            download('ftp://ftp.ensembl.org/pub/release-' + ens_v + '/mysql/caenorhabditis_elegans_core_' + ens_v + '_260/',
                     core)
        elif (int(ens_v) >= 74 and int(ens_v) <= 90):
            core='/caenorhabditis_elegans_core_' + ens_v + '_250.sql.gz'
            download('ftp://ftp.ensembl.org/pub/release-' + ens_v +'/mysql/caenorhabditis_elegans_core_' + ens_v +'_250/',core)
        else:
            print(("ERROR: unsupported ensembl version: " + ens_v))
            print("supported ensembl versions: from 74 till 90")
            sys.exit()
    elif (species=='danio_rerio' or species =="zebrafish"):
        if(int(ens_v) >= 74 and int(ens_v) <= 90):
            core='/danio_rerio_core_' + ens_v + '_10.sql.gz'
            download('ftp://ftp.ensembl.org/pub/release-' + ens_v +'/mysql/danio_rerio_core_' + ens_v +'_10/',core)
        else:
            print(("ERROR: unsupported ensembl version: " + ens_v))
            print("supported ensembl versions: from 74 till 88")
            sys.exit()
    elif (species=='SL1344' or species=='salmonella_eentrica_serovar_typhimurium'):
        core = 'bacteria_23_collection_core_' + ens_v +'_' + canEns_v + '_1.sql.gz'
        download('ftp://ftp.ensemblgenomes.org/pub/release-' + ens_v +'/bacteria//mysql/bacteria_23_collection_core_' + ens_v +'_' + canEns_v + '_1/', core)
    elif (species=='MYC_ABS_ATCC_19977' or species=='mycobacterium_abscessus_atcc_19977'):
        core = 'bacteria_16_collection_core_' + ens_v +'_' + canEns_v + '_1.sql.gz'
        download('ftp://ftp.ensemblgenomes.org/pub/release-' + ens_v +'/bacteria//mysql/bacteria_16_collection_core_' + ens_v +'_' + canEns_v + '_1/', core)
    elif (species=='CNECNA3' or species=='Cryptococcus_neoformans_var_grubii_h99_gca_000149245'):
        core = 'fungi_basidiomycota1_collection_core_' + ens_v +'_' + canEns_v + '_1.sql.gz'
        download('ftp://ftp.ensemblgenomes.org/pub/release-' + ens_v +'/fungi/mysql/fungi_basidiomycota1_collection_core_' + ens_v +'_' + canEns_v + '_1/', core)
    elif (species=='physcomitrium_patens' or species=='earthmoss'):
        if (int(ens_v) <= 51):
            core = 'physcomitrium_patens_core_'+ens_v+'_'+canEns_v+'_2.sql.gz'
            download('ftp://ftp.ensemblgenomes.org/pub/release-' + ens_v +'/plants/mysql/physcomitrium_patens_core_' + ens_v +'_' + canEns_v + '_2/', core)

    else:
        print(("Error: unsupported species: " +species ))
        print("Supported species: human, fruitfy, mouse, horse, arctic_squirrel, saccharomyces_cerevisiae, caenorhabditis_elegans, SL1344, MYC_ABS_ATCC_19977, CNECNA3, Cryptococcus_neoformans_var_grubii_h99_gca_000149245")
        sys.exit()

    #
    # Change directory to ENS folder in tmp
    # unzip downloaded files
    #

    os.chdir(directory + '/tmp/ENS')
    try:
        os.system('gzip -d *.gz')
    except OSError as e:
        if e.errno == os.errno.ENOENT:
             print(" wget was not found on your system, please install gzip in order to run this program")
        else:
            print("something went wrong trying to run gzip")
            raise

    #From ensembl 109 onwards: MySQL Row formatting option that does not fit SQLite -> delete these options from dump file
    if(int(ens_v)>=109):
        os.system("sed -i 's/ ROW_FORMAT=DYNAMIC//g' *_"+ens_v+"_*.sql")
        os.system("sed -i 's/ ROW_FORMAT=FIXED//g' *_"+ens_v+"_*.sql")

    #Edit collating sequence
    os.system("sed -i 's/COLLATE latin1_bin/COLLATE BINARY/g' *_"+ens_v+"_*.sql")

    #load db schema
    file = open("mysql2sqlite.sh", "w")
    file.write("""
    #!/usr/bin/env bash
    cat $1 |
    grep -v 'LOCK' |
    grep -v ' KEY ' |
    grep -v ' UNIQUE KEY ' |
    grep -v ' PRIMARY KEY ' |
    perl -pe 's/ ENGINE[ ]*=[ ]*[A-Za-z_][A-Za-z_0-9]*(.*DEFAULT)?/ /gi' |
    perl -pe 's/ CHARSET[ ]*=[ ]*[A-Za-z_][A-Za-z_0-9]*/ /gi' |
    perl -pe 's/ CHECKSUM=[A-Za-z_0-9]*/ /gi' |
    perl -pe 's/ MAX_ROWS=[A-Za-z_0-9]*/ /gi' |
    perl -pe 's/ AVG_ROW_LENGTH=[A-Za-z_0-9]*/ /gi' |
    perl -pe 's/ [ ]*AUTO_INCREMENT=[0-9]* / /gi' |
    perl -pe 's/ unsigned / /g' |
    perl -pe 's/ set[(][^)]*[)] / varchar(255) /gi' |
    perl -pe 's/ auto_increment/ primary key autoincrement/gi' |
    perl -pe 's/ smallint[(][0-9]*[)] / integer /gi' |
    perl -pe 's/ tinyint[(][0-9]*[)] / integer /gi' |
    perl -pe 's/ int[(][0-9]*[)] / integer /gi' |
    perl -pe 's/ character set [^ ]* / /gi' |
    perl -pe 's/ enum[(][^)]*[)] / varchar(255) /gi' |
    perl -pe 's/ on update [^,]*//gi' |
    perl -e 'local $/;$_=<>;s/,\\n\)/\\n\)/gs;print "begin;\\n";print;print "commit;\\n"' |
    perl -pe '
    if (/^(INSERT.+?)\(/) {
       $a=$1;
       s/\\\\'\\''/'\\'\\''/g;
       s/\\\\n/\\n/g;
       s/\),\(/\);\\n$a\(/g;
    }
    '
        """)
    file.close()

    #set permission for executeable file
    st = os.stat('mysql2sqlite.sh')
    os.chmod('mysql2sqlite.sh', st.st_mode | stat.S_IEXEC)

    try:
        os.system('./mysql2sqlite.sh *.sql | sqlite3 '+directory + '/ENS_' + speciesdict[species] + '_' + ens_v + '.db')
    except OSError as e:
        print("ERROR: could not import database architecture")


    #
    # create SQLite3 DB and estabilish connection
    #

    (' Start importing and parsing data \n this can take a while')
    print('Creating sqlite3 database')
    conn = sqlite3.connect(directory + '/ENS_' + speciesdict[species] + '_' + ens_v + '.db')
    cur = conn.cursor()

    with open(directory + '/tmp/ENS/coord_system.txt','rt') as inputfile:
        reader= csv.reader(inputfile,delimiter="\t")
        first_row=next(reader)
        columns=len(first_row)
        columns_string="(?"
        for j in range (0,columns-1):
            columns_string=columns_string+",?"
        columns_string=columns_string+')'
        inputfile.seek(0)
        for row in reader:
            cur.execute("INSERT INTO coord_system VALUES"+columns_string,row)
    print("Finished importing coord_system table")

    with open(directory + '/tmp/ENS/exon.txt','rt') as inputfile:
        reader= csv.reader(inputfile,delimiter="\t")
        first_row=next(reader)
        columns=len(first_row)
        columns_string="(?"
        for j in range (0,columns-1):
            columns_string=columns_string+",?"
        columns_string=columns_string+')'
        inputfile.seek(0)
        for row in reader:
            cur.execute("INSERT INTO exon VALUES"+columns_string,row)
    print("Finished importing exon table")

    with open(directory + '/tmp/ENS/exon_transcript.txt','rt') as inputfile:
        reader= csv.reader(inputfile,delimiter="\t")
        first_row=next(reader)
        columns=len(first_row)
        columns_string="(?"
        for j in range (0,columns-1):
            columns_string=columns_string+",?"
        columns_string=columns_string+')'
        inputfile.seek(0)
        for row in reader:
            cur.execute("INSERT INTO exon_transcript VALUES"+columns_string,row)
    print("Finished importing exon_transcript table")

    with open(directory + '/tmp/ENS/gene.txt','rt') as inputfile:
        reader= csv.reader(inputfile,delimiter="\t")
        first_row=next(reader)
        columns=len(first_row)
        columns_string="(?"
        for j in range (0,columns-1):
            columns_string=columns_string+",?"
        columns_string=columns_string+')'
        inputfile.seek(0)
        for row in reader:
            cur.execute("INSERT INTO gene VALUES"+columns_string,row)
    print("Finished importing gene table")

    with open(directory + '/tmp/ENS/seq_region.txt','rt') as inputfile:
        reader= csv.reader(inputfile,delimiter="\t")
        first_row=next(reader)
        columns=len(first_row)
        columns_string="(?"
        for j in range (0,columns-1):
            columns_string=columns_string+",?"
        columns_string=columns_string+')'
        inputfile.seek(0)
        for row in reader:
            cur.execute("INSERT INTO seq_region VALUES"+columns_string,row)
    print("Finished importing seq_region table")

    with open(directory + '/tmp/ENS/transcript.txt','rt') as inputfile:
        reader= csv.reader(inputfile,delimiter="\t")
        first_row=next(reader)
        columns=len(first_row)
        columns_string="(?"
        for j in range (0,columns-1):
            columns_string=columns_string+",?"
        columns_string=columns_string+')'
        inputfile.seek(0)
        for row in reader:
            cur.execute("INSERT INTO transcript VALUES"+columns_string,row)
    print("Finished importing transcript table")

    with open(directory + '/tmp/ENS/translation.txt','rt') as inputfile:
        reader= csv.reader(inputfile,delimiter="\t")
        first_row=next(reader)
        columns=len(first_row)
        columns_string="(?"
        for j in range (0,columns-1):
            columns_string=columns_string+",?"
        columns_string=columns_string+')'
        inputfile.seek(0)
        for row in reader:
            cur.execute("INSERT INTO translation VALUES"+columns_string,row)
    print("Finished importing translation table")


    conn.text_factory = str

    with open(directory + '/tmp/ENS/object_xref.txt','rt') as inputfile:
        reader= csv.reader(inputfile,delimiter="\t")
        first_row=next(reader)
        columns=len(first_row)
        columns_string="(?"
        for j in range (0,columns-1):
            columns_string=columns_string+",?"
        columns_string=columns_string+')'
        inputfile.seek(0)
        for row in reader:
            cur.execute("INSERT INTO object_xref VALUES"+columns_string,row)
    print("Finished importing object_xref table")

    with open(directory + '/tmp/ENS/xref.txt','rt') as inputfile:
        reader= csv.reader(inputfile,delimiter="\t")
        first_row=next(reader)
        columns=len(first_row)
        columns_string="(?"
        for j in range (0,columns-1):
            columns_string=columns_string+",?"
        columns_string=columns_string+')'
        inputfile.seek(0)
        for row in reader:
            if columns==8:
                if len(row)==1:
                    break;
                elif len(row)==6:
                    row.append(" ")
                    row.append(" ")
                    cur.execute("INSERT INTO xref VALUES"+columns_string,row)
                else:
                    cur.execute("INSERT INTO xref VALUES"+columns_string,row)
            else:
                cur.execute("INSERT INTO xref VALUES"+columns_string,row)
    print("Finished importing xref table")

    conn.commit()

    print('Ensembl gene annotation database creation successful')
    shutil.rmtree(directory + '/tmp/ENS')

    return

def print_help():

    help="""
ARGUMENTS:

-v | --version                              Ensembl annotation version to download (supported versions: from 74)
-s | --species                              Specify the desired species for which gene annotation files should be downloaded
                                            currently supported species:

                                            human                       |   homo_sapiens
                                            mouse                       |   mus_musculus
                                            rat                         |   rattus_norvegicus
                                            horse                       |   equus_caballus
                                            chinese_hamster             |   Cricetulus griseus
                                            arctic_squirrel             |   urocitellus_parryii
                                            fruitfly                    |   drosophila_melanogaster
                                            saccharomyces_cerevisiae    |   yeast
                                            caenorhabditis_elegans      |   c.elegans
                                            SL1344                      |   Salmonella enterica subsp. enterica serovar Typhimurium str. SL1344
                                            MYC_ABS_ATCC_19977          |   Mycobacterium abscessus atcc 19977
                                            CNECNA3                     |   Cryptococcus_neoformans_var_grubii_h99_gca_000149245
                                            earthmoss                   |   physcomitrium_patens

-h | --help                                 Print this useful help message


The Ensembl database will be installed in the current directory ( directory where the script is run from )

EXAMPLE:

python ENS_db.py -v 98 -s human

DEPENDENCIES:

This program depends upon wget and gzip, which are normally pre-installed on any unix system

    """

    print(help)
    print()

    return



###### DIRECT TO MAIN ############
if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        traceback.print_exc()
##################################
