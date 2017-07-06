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
                                            fruitfly                    |   drosophila_melanogaster
                                            saccharomyces_cerevisiae    |   yeast
                                            caenorhabditis_elegans      |   c.elegans


The Ensembl database will be installed in the current directory ( directory where the script is run from )

EXAMPLE:

python ENS_db.py -v 78 -s human

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

try:
    myopts, args = getopt.getopt(sys.argv[1:],"s:v:",["version=","species="])
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
print("Ensembl version used : " + ens_v)
print("Selected species     : " + species)

directory = os.getcwd()
os.chdir(directory)

if not os.path.exists("tmp"):
    os.system('mkdir tmp')

if os.path.isdir("tmp/ENS"):
    shutil.rmtree('tmp/ENS')

os.system('mkdir tmp/ENS')

#
# Dictionary for 3letter abbreviation for species
#

speciesdict = {'human': 'hsa', 'mouse': 'mmu', 'fruitfly': 'dme', 'yeast': 'sce', 'c.elegans': 'cel',
               'homo_sapiens': 'hsa', 'mus_musculus': 'mmu', 'drosophila_melanogaster': 'dme',
               'saccharomyces_cerevisiae': 'sce', 'caenorhabditis_elegans': 'cel'};


#
# Check whether database already exists
#

if os.path.isfile(directory + '/ENS_' + species + '_v_' + ens_v + '.db'):
    print("There already exists a file called " + 'ENS_' + speciesdict[species] + '_' + ens_v + '.db in: ' + directory)
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

if (species=='human' or species=='homo_sapiens'):
    if(int(ens_v) >= 76 and int(ens_v) <= 88):
         core='/homo_sapiens_core_' + ens_v + '_38.sql.gz'
         download('ftp://ftp.ensembl.org/pub/release-'+ ens_v +'/mysql/homo_sapiens_core_' + ens_v + '_38/',core)
    elif(int(ens_v) >= 74 and int(ens_v) <= 75):
         core='/homo_sapiens_core_' + ens_v + '_37.sql.gz'
         download('ftp://ftp.ensembl.org/pub/release-'+ ens_v +'/mysql/homo_sapiens_core_' + ens_v + '_37/',core)
    else:
        print("ERROR: unsupported ensembl version: " + ens_v)
        print("supported ensembl versions: from 74 till 88")
        sys.exit()
elif (species=='mouse' or species=='mus_musculus'):
    if(int(ens_v) >= 74 and int(ens_v) <= 88):
        core='/mus_musculus_core_' + ens_v + '_38.sql.gz'
        download('ftp://ftp.ensembl.org/pub/release-'+ ens_v + '/mysql/mus_musculus_core_' + ens_v +'_38/',core)
    else:
        print("ERROR: unsupported ensembl version: " + ens_v)
        print("supported ensembl versions: from 75 till 88")
        sys.exit()
elif (species=='fruitfly' or species=='drosophila_melanogaster'):
    if(int(ens_v) >= 74 and int(ens_v) <= 88):
        core='/drosophila_melanogaster_core_' + ens_v + '_546.sql.gz'
        download('ftp://ftp.ensembl.org/pub/release-' + ens_v +'/mysql/drosophila_melanogaster_core_' + ens_v +'_546/',core)
    else:
        print("ERROR: unsupported ensembl version: " + ens_v)
        print("supported ensembl versions: from 74 till 88 ")
        sys.exit()
elif (species == 'saccharomyces_cerevisiae' or species == 'yeast'):
    if(int(ens_v) >= 74 and int(ens_v) <= 88):
        core='/saccharomyces_cerevisiae_core_' + ens_v + '_4.sql.gz'
        download('ftp://ftp.ensembl.org/pub/release-' + ens_v +'/mysql/saccharomyces_cerevisiae_core_' + ens_v +'_4/',core)
    else:
        print("ERROR: unsupported ensembl version: " + ens_v)
        print("supported ensembl versions: from 74 till 88")
        sys.exit()
elif (species=='caenorhabditis_elegans' or species =="c.elegans"):
    if(int(ens_v) >= 74 and int(ens_v) <= 88):
        core='/caenorhabditis_elegans_core' + ens_v + '_245.sql.gz'
        download('ftp://ftp.ensembl.org/pub/release-' + ens_v +'/mysql/caenorhabditis_elegans_core_' + ens_v +'_245/',core)
    else:
        print("ERROR: unsupported ensembl version: " + ens_v)
        print("supported ensembl versions: from 74 till 88")
        sys.exit()
elif (species=='danio_rerio' or species =="zebrafish"):
    if(int(ens_v) >= 74 and int(ens_v) <= 88):
        core='/danio_rerio_core_' + ens_v + '_9.sql.gz'
        download('ftp://ftp.ensembl.org/pub/release-' + ens_v +'/mysql/danio_rerio_core_' + ens_v +'_9/',core)
    else:
        print("ERROR: unsupported ensembl version: " + ens_v)
        print("supported ensembl versions: from 74 till 88")
        sys.exit()
else:
    print("Error: unsupported species: " +species )
    print("Supported species: human, fruitfy, mouse, saccharomyces_cerevisiae, caenorhabditis_elegans")
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

#load db schema
file = open("mysql2sqlite.sh", "w+r")
file.write("""
#!/usr/bin/env bash
cat $1 |
grep -v 'LOCK' |
grep -v ' KEY ' |
grep -v ' UNIQUE KEY ' |
grep -v ' PRIMARY KEY ' |
perl -pe 's/ ENGINE[ ]*=[ ]*[A-Za-z_][A-Za-z_0-9]*(.*DEFAULT)?/ /gi' |
perl -pe 's/ CHARSET[ ]*=[ ]*[A-Za-z_][A-Za-z_0-9]*/ /gi' |
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

with open(directory + '/tmp/ENS/coord_system.txt','rb') as inputfile:
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

with open(directory + '/tmp/ENS/exon.txt','rb') as inputfile:
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

with open(directory + '/tmp/ENS/exon_transcript.txt','rb') as inputfile:
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

with open(directory + '/tmp/ENS/gene.txt','rb') as inputfile:
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

with open(directory + '/tmp/ENS/seq_region.txt','rb') as inputfile:
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

with open(directory + '/tmp/ENS/transcript.txt','rb') as inputfile:
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

with open(directory + '/tmp/ENS/translation.txt','rb') as inputfile:
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

with open(directory + '/tmp/ENS/object_xref.txt','rb') as inputfile:
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

with open(directory + '/tmp/ENS/xref.txt','rb') as inputfile:
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


