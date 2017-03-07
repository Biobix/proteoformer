__author__ = 'gurb'
#Based on previous work of 'Steven Verbruggen, cfr. ORFbasedCounts.py'

'''

ARGUMENTS:

    -w | --work_dir						 The working directory (default = current directory)
    -c | --sqliteC						 The sqlite database holding the count data. (mandatory)
    -e | --sqliteE						 The sqlite database holding the Ensembl annotation (mandatory)
    -f | --feature				         The annotation feature at what level to summarize. Can be set to 'exon','transcript','gene', 'promoter' (mandatory)
    -t | --transcripts                   The type of transcripts to process. Can be set to 'all','canonical' (default = 'all')
    -d | --data                          The type of data used to generate the counts. Can be set to 'ribo', 'rna', 'rrbs' (mandatory)
    -m | --mapping                       The type of read mapping. Can be set to 'unique','multi' (default = 'unique' for ribo or rna, not specified for rrbs)
    -r | --rrbs                          The table holding the rrbs count data (mandatory if feature='promoter')
	-o | --five_prime_offset			 The five_prime_offset for feature boundaries, used for gene features (default = 0)

EXAMPLE:

python summarize_feature.py --sqliteC SQLite/results.db --sqliteE SQLite/ENS_hs_82.db --feature transcript --transcripts all

DEPENDENCIES:
sqlite3


'''



import sqlite3 as sqlite
import csv
import os
import getopt
from multiprocessing import Pool
from collections import defaultdict
import traceback
import sys
import re
import timeit
import math

def main():

    print "		        \t------------------------------------------"
    print "				| summarize counts for annotated regions |"
    print "				------------------------------------------"
    print "\n"
    print "\t\tThis program fetches counts data for annotated regions\n"
    print "\t\tThe count data can be based on RNA-seq, RIBO-seq, RRBS, ...\n"
    print "\t\tThe annotation  - is fetched from an Ensembl SQLite DB\n"
    print "\t\t                - can be gene, transcript, exon, ...\n"
    ### Gets and Sets

    # #Timer
    startTime = timeit.default_timer()

    #Init
    sqliteC=''
    sqliteE=''
    feature=''

    #Catch command line with getopt
    try:
        myopts, args = getopt.getopt(sys.argv[1:],"w:c:e:f:t:m:r:d:o:",["work_dir=","sqliteC=","sqliteE=","feature=","transcripts=","mapping=","rrbs=","data=","five_prime_offset="])
    except getopt.GetoptError as err:
        print err
        sys.exit()

    # Catch arguments
    # o == option
    #  a == argument passed to the o
    for o, a in myopts:
        if o in ('-w','--work_dir'):
            workdir=a
        elif o in ('-c','--sqliteC'):
            sqliteC=a
        elif o in ('-e','--sqliteE'):
            sqliteE=a
        elif o in ('-f','--feature'):
            feature=a
        elif o in ('-t','--transcripts'):
            transcripts=a
        elif o in ('-m','--mapping'):
            mapping=a
        elif o in ('-r','--rrbs'):
            rrbs=a
        elif o in ('-d','--data'):
            data=a
        elif o in ('-o','--five_prime_offset'):
            five_prime_offset=a
    try:
        workdir
    except:
        workdir=''
    try:
        sqliteC
    except:
        sqliteC=''
    try:
        sqliteE
    except:
        sqliteE=''
    try:
        feature
    except:
        feature=''
    try:
        transcripts
    except:
        transcripts=''
    try:
        mapping
    except:
        mapping=''
    try:
        rrbs
    except:
        rrbs=''
    try:
        data
    except:
        data=''
    try:
        five_prime_offset
    except:
        five_prime_offset=0


    # Check for correct arguments and parse
    if workdir == '':
        workdir = os.getcwd()
    if workdir !='':
        os.chdir(workdir)
    if sqliteC=='':
        print "Error: do not forget to pass the SQLite Counts argument!"
        sys.exit()
    if sqliteE=='':
        print "Error: do not forget to pass the SQLite Ensembl argument!"
        sys.exit()
    if data=='':
        print "Error: do not forget to pass the sequencing data type: 'rna', 'ribo', 'rrbs'!"
        sys.exit()
    if feature=='':
        print "Error: do not forget to pass the feature argument: 'exon', 'transcript' or 'gene'!"
        sys.exit()
    if feature != 'gene':
        if transcripts == '':
            transcripts = 'all'
        if transcripts != '' and transcripts != 'canonical' and transcripts != 'all':
            print "Error: the following values are allowed for --transcripts : 'all', 'canonical'!"
            sys.exit()
    if mapping == '' and data != 'rrbs': #Set default value
        mapping = 'unique'
    if mapping != '' and mapping != 'unique' and mapping != 'multi' and data != 'rrbs': #Check validity of mapping value for non-rrbs data
        print "Error: the following values are allowed for --mapping : 'unique', 'multi' (for RNA- or RIBO-seq data)!"
    if mapping != '' and data == 'rrbs':
        print "Error: the  --mapping should not be specified for RRBS data!"
        sys.exit()
    if data == 'rrbs' and rrbs == '':
        print "Error: the table holding the RRBS data should be given --rrbs : 'table_name'!"
        sys.exit()

    #Get arguments from arguments table
    igenomes_root, species, ens_version, cores, readtype = get_arguments(sqliteC)



    # Print used arguments and parameters
    print "Arguments:"
    print "  Working directory:            "+workdir
    print "  SQLite Ensembl:               "+sqliteE
    print "  SQLite Count data:            "+sqliteC
    print "  Sequencing data:              "+data
    print "  RRBS data table:              "+rrbs
    print "  Transcripts:                  "+transcripts
    print "  Mapping:                      "+mapping
    print "  Feature:                      "+feature
    print "  5-prime offset:			   "+str(five_prime_offset)


    directory = os.getcwd()
    print directory
    os.chdir(directory)

    #Make tmp folder if it not exists
    if(not os.path.isdir(directory+"/tmp")):
        os.mkdir(directory+"/tmp",)

    #Conversion of species terminology
    speciesLatin="Mus_musculus" if species=="mouse" else \
        "Homo_sapiens" if species=="human" else \
        "Arabidopsis_thaliana" if species=="arabidopsis" else \
        "Drosophila_melanogaster" if species=="fruitfly" else ""
    speciesShort="mmu" if species=="mouse" else \
        "hsa" if species=="human" else \
        "ath" if species=="arabidopsis" else \
        "dme" if species=="fruitfly" else ""
    #Assembly
    assembly="GRCm38" if species=="mouse" and ens_version>=70 else \
        "NCBIM37" if species=="mouse" and ens_version<70 else \
        "GRCh38" if species=="human" and ens_version>75 else \
        "GRCh37" if species=="human" and ens_version<=75 else \
        "TAIR10" if species=="arabidopsis" else \
        "BDGP5" if species=="fruitfly" else ""

    #Get chromosomes
    chrs={}
    chromosomeSizesFile = igenomes_root+"/"+speciesLatin+"/Ensembl/"+assembly+"/Annotation/Genes/ChromInfo.txt"
    if os.path.isfile(chromosomeSizesFile):
        chrs = get_chrs(chromosomeSizesFile, species)
    else:
        chromosomeSizesFile = igenomes_root+"/"+speciesLatin+"/Ensembl/"+assembly+"/Sequence/WholeGenomeFasta/GenomeSize.xml"
        if os.path.isfile(chromosomeSizesFile):
            chrs = get_chrsXml(chromosomeSizesFile, species)
        else:
            print "ERROR: chromosome sizes file could not be found."
            sys.exit()
    #Get coordinate system id for chromosomes
    coordSystemId = get_coord_system_id_chr(sqliteE, assembly)

    #Get sequencing type
    #The readtype can be ribo, PE_polyA, SE_polyA, PE_total or SE_total
    pattern = re.compile('SE|PE')
    m = pattern.search(readtype)
    if m:
        seqtype = 'rna'
    else:
        seqtype = 'ribo'
    if data == 'rrbs':
        seqtype = 'rrbs'

    if feature == 'gene':
        table_name = feature+"_"+seqtype+"_"+mapping+"_count"
    elif feature in ('transcript','exon') and data in ('rna','ribo'):
        table_name = feature+"_"+seqtype+"_"+mapping+"_"+transcripts+"_count"
    elif feature in ('transcript','exon') and data in ('rrbs'):
        table_name = feature+"_"+seqtype+"_"+transcripts+"_count"
    elif feature == 'promoter':
        table_name = feature+"_"+seqtype+"_"+transcripts+"_count"


    ##############
    ###  MAIN  ###
    ##############
    #chrs = ['Y']
    pool = Pool(processes=cores)
    [pool.apply_async(feature_count_chr, args=(chrom,sqliteE,sqliteC,directory,feature,transcripts,table_name,mapping,rrbs,data,five_prime_offset)) for chrom in chrs.keys()]
    #[pool.apply_async(feature_count_chr, args=(chrom,sqliteE,sqliteC,directory,feature,transcripts,table_name,mapping,rrbs,five_prime_offset)) for chrom in chrs]
    pool.close()
    pool.join()
    
    #sys.exit()
    #Dump in SQLite
    print
    print "Dump data into SQLite database"
    print
    
	
    dumpSQLite(chrs.keys(), directory, sqliteC,table_name,data)
    #dumpSQLite(['Y'], directory, sqliteC,table_name,data)

    #Timer
    stopTime = timeit.default_timer()
    runtime = math.ceil(stopTime - startTime)
    m,s = divmod(runtime, 60)
    h,m = divmod(m,60)
    runtime = "%d:%02d:%02d" % (h,m,s)
    print "\n   ------ Program COMPLETED ------"
    print "   The program run time: "+runtime+"\n\n"

    return

##############
###  SUBS  ###
##############

## Get counts for features for chr ##
def feature_count_chr(chrom,sqliteE,sqliteC,directory,feature,transcripts,table_name,mapping,rrbs,data,five_prime_offset):
    try:

        #print(chrom)
        #Get features per chr
        features = get_feature_chr(chrom, sqliteE, feature, transcripts)
        #print_dict(features)

        #Get all counts per chr (both forward and reverse)
        countsforw, countsrev = counts_chr(chrom, sqliteC, mapping, rrbs, data)
        #print_dict(countsforw)


        #Match counts to feature region
        features = calculate_count_feature(features, countsforw, countsrev, feature,five_prime_offset)
        #print_dict(features)

        #Write to csv file for each chromosome
        write_to_csv(chrom, features, directory, feature, table_name)
        print "*) Process for chromosome "+str(chrom)+" completed"
    except:
            traceback.print_exc()
    return

## Get exons per chromosome ##
def get_feature_chr(chrom, sqliteE,feature,transcripts):

    try:
        #print(chrom)
        if feature == 'exon':
            query = "select ex.stable_id feature_id, ex.seq_region_start feature_start,ex.seq_region_end feature_end, sr.name chrom, ex.seq_region_strand strand " \
                    ", ge.stable_id gene_id, tr.biotype, tr.stable_id transcript_id, ext.rank exon_rank "\
                    "from transcript tr inner join gene ge on ge.gene_id = tr.gene_id "

            if transcripts == 'canonical':
                    query = query + "and ge.canonical_transcript_id = tr.transcript_id "

            query = query + "inner join seq_region sr on sr.seq_region_id = tr.seq_region_id " \
                            "inner join exon_transcript ext on ext.transcript_id = tr.transcript_id " \
                            "inner join exon ex on ex.exon_id = ext.exon_id " \
                            "where sr.name = '"+str(chrom)+"'"

        #The transcript scripts also fetches a list of start/stops of the exons included in the transcript
        #The group_concat sorts in descending order (for SQLite, that is...)
        elif feature == 'transcript':
            query = "select tr.stable_id feature_id, tr.seq_region_start feature_start,tr.seq_region_end feature_end, sr.name chrom, tr.seq_region_strand strand " \
                    ", ge.stable_id gene_id, tr.biotype " \
                    ", group_concat(ex.seq_region_start) starts, group_concat(ex.seq_region_end) ends, count(ext.rank) exon_count " \
                    "from transcript tr inner join gene ge on ge.gene_id = tr.gene_id "

            if transcripts == 'canonical':
                    query = query + "and ge.canonical_transcript_id = tr.transcript_id "

            query = query + "inner join seq_region sr on sr.seq_region_id = tr.seq_region_id " \
                            "inner join exon_transcript ext on ext.transcript_id = tr.transcript_id " \
                            "inner join exon ex on ex.exon_id = ext.exon_id " \
                            "where sr.name = '"+str(chrom)+"'" \
                            " group by tr.stable_id"


        elif feature == 'gene':
            query = "select ge.stable_id feature_id, ge.seq_region_start feature_start, ge.seq_region_end feature_end, sr.name chrom, ge.seq_region_strand strand " \
                    "from gene ge " \
                    "inner join seq_region sr on sr.seq_region_id = ge.seq_region_id " \
                    "where sr.name = '"+str(chrom)+"'"

        elif feature == 'promoter':
            query = "select tr.stable_id feature_id, tr.seq_region_start feature_start,tr.seq_region_end feature_end, sr.name chrom, tr.seq_region_strand strand " \
                    ", ge.stable_id gene_id, tr.biotype " \
                    "from transcript tr inner join gene ge on ge.gene_id = tr.gene_id "

            if transcripts == 'canonical':
                    query = query + "and ge.canonical_transcript_id = tr.transcript_id "

            query = query + " inner join seq_region sr on sr.seq_region_id = tr.seq_region_id " \
                            "where sr.name = '"+str(chrom)+"'"

        #print(query)
        features = fetchDict(sqliteE, query)


        #Initialize the read count already
        for featId in features:
            features[featId]['counts_sum']=0 #Gives the total count
            features[featId]['counts_n']=0   #Gives the number of positions covered

    except:
        traceback.print_exc()

    return features

## Get counts from count database ##
def counts_chr(chrom, sqliteC,mapping,rrbs,data):

    #Initialize
    countsfor = {}
    countsrev = {}

    if data in ('rna','ribo'):
        try:

            #Count table based on unique/multi mapper reads
            if mapping == 'unique':
                count_table = 'count_fastq1_unique'
            elif mapping == 'multi':
                count_table = 'count_fastq1'

            #Fetch data from count tables
            query = "SELECT start, count FROM "+count_table+" WHERE chr='"+str(chrom)+"' AND strand = '1';"
            countsfor = fetchDict(sqliteC, query)
            query = "SELECT start, count FROM "+count_table+" WHERE chr='"+str(chrom)+"' AND strand = '-1';"
            countsrev = fetchDict(sqliteC, query)

        except:
            traceback.print_exc()

    if data == 'rrbs':
        try:

            #CpG_forw is aliased to 'start'
            #meth_perc is aliased to 'count'
            #This to make the generic (i.e. similar to other counts data (ribo/rna))
            query = "select CpG_forw start, meth_perc count " \
            "from "+rrbs+" where RRBS_coverage >= 10 and chr = '"+str(chrom)+"'"
            countsfor = fetchDict(sqliteC, query)

        except:
            traceback.print_exc()

    return countsfor, countsrev


#Get the counts for the features ##
def calculate_count_feature(features, countsforw,countsrev,feature,five_prime_offset):

    try:


        if feature == 'promoter':

            features_for=[]
            features_rev=[]
            for id in features:
                if(features[id]['strand']==1):
                    features_for.append(id)
                elif(features[id]['strand']==-1):
                    features_rev.append(id)
            features_for = sorted(features_for, key=lambda x: features[x]['feature_start'])
            features_rev = sorted(features_rev, key=lambda x: features[x]['feature_end'])

            #This makes 'sense' (forward count positions)
            #Init
            window = []

            #Go through all sorted read positions (only forward counts since RRBS DNA-methylation data is sense-independent)
            for readPosition in sorted(countsforw):
                #print(readPosition)
                for featId in features_for[:]: #Go through sorted features list. [:] means working on a copy of the features_for list, because if you would not, features will move forward in the list if you remove the first feature and in the next iteration, you examine the third feature (which moved to the second place in the list) and you will have skipped the second feature (no being on the first place). Therefore, work on a copy of the original list but remove in the original list. For each new read position however, the copy will be replaced as you start a new iteration over features_for
                    #print(featId)
                    if features[featId]['feature_start']-1000<=readPosition:
                        #print " append "+ featId+" "+ str(features[featId]['feature_start']-1000)
                        window.append(featId) #Bring the feature ID over to the window if the read >= feature_start
                        features_for.remove(featId) #Remove from sorted features so it will not be scanned next time
                    else:
                        break #Because features are sorted, all the other features will be more downstream
                for featId in window[:]:#For all features in the window
                    if features[featId]['feature_start']+200<readPosition:
                        #print" remove "+ featId+" "+ str(features[featId]['feature_start']+200)
                        window.remove(featId) #Remove features from window which are completely upstream of the read
                    else: #Only features with the read in between start and stop remain in the window
                        features[featId]['counts_sum']  = features[featId]['counts_sum'] + countsforw[readPosition]['count']
                        features[featId]['counts_n'] += 1
                        #print" sum "+ featId+" "+ str(countsforw[readPosition]['count'])+ " added to sum="+str(features[featId]['counts_sum'])+ " count="+str(features[featId]['counts_n'])

            #Empty window
            window=[]

            #This makes 'antisense'
            #Go through all sorted read positions (only forward counts since RRBS DNA-methylation data is sense-independent)
            for readPosition in sorted(countsforw):
                #print(readPosition)
                for featId in features_rev[:]: #Go through sorted features list. [:] means working on a copy of the features_rev list, because if you would not, features will move forward in the list if you remove the first feature and in the next iteration, you examine the third feature (which moved to the second place in the list) and you will have skipped the second feature (no being on the first place). Therefore, work on a copy of the original list but remove in the original list. For each new read position however, the copy will be replaced as you start a new iteration over features_rev
                    #print(featId)
                    if features[featId]['feature_end']-200<=readPosition:
                        #print " append "+ featId+" "+ str(features[featId]['feature_end'])
                        window.append(featId) #Bring the feature ID over to the window if the read >= feature_end
                        features_rev.remove(featId) #Remove from sorted features so it will not be scanned next time
                    else:
                        break #Because features are sorted, all the other features will be more downstream
                for featId in window[:]:#For all features in the window
                    if features[featId]['feature_end']+1000<readPosition:
                        #print" remove "+ featId+" "+ str(features[featId]['feature_end'])
                        window.remove(featId) #Remove features from window which are completely upstream of the read
                    else: #Only features with the read in between start and stop remain in the window
                        features[featId]['counts_sum']  = features[featId]['counts_sum'] + countsforw[readPosition]['count']
                        features[featId]['counts_n'] += 1
                        #print" sum "+ featId+" "+ str(countsforw[readPosition]['count'])+ " added to sum="+str(features[featId]['counts_sum'])+ " count="+str(features[featId]['counts_n'])

            #Empty window
            window=[]

        elif feature in ('gene','transcript','exon'):

            features_for=[]
            features_rev=[]
            for id in features:
                if(features[id]['strand']==1):
                    features_for.append(id)
                elif(features[id]['strand']==-1):
                    features_rev.append(id)
            features_for = sorted(features_for, key=lambda x: features[x]['feature_start'])
            features_rev = sorted(features_rev, key=lambda x: features[x]['feature_start'])

            #This makes 'sense'
            #Init
            window = []
            #print 'sense'

            #Go through all sorted read positions
            for readPosition in sorted(countsforw):
                #print(readPosition)
                for featId in features_for[:]: #Go through sorted features list. [:] means working on a copy of the features_sorted list, because if you would not, features will move forward in the list if you remove the first feature and in the next iteration, you examine the third feature (which moved to the second place in the list) and you will have skipped the second feature (no being on the first place). Therefore, work on a copy of the original list but remove in the original list. For each new read position however, the copy will be replaced as you start a new iteration over features_sorted
                    #print(featId)
                    if features[featId]['feature_start']-int(five_prime_offset)<=readPosition:
                        #print " append "+ featId+" "+ str(features[featId]['feature_start'])
                        window.append(featId) #Bring the feature ID over to the window if the read >= feature_start
                        features_for.remove(featId) #Remove from sorted features so it will not be scanned next time
                    else:
                        break #Because features are sorted, all the other features will be more downstream
                for featId in window[:]:#For all features in the window
                    if features[featId]['feature_end']<readPosition:
                        #print" remove "+ featId+" "+ str(features[featId]['feature_end'])
                        window.remove(featId) #Remove features from window which are completely upstream of the read
                    else: #Only features with the read in between start and stop remain in the window

                        # Only features within the exons are counted
                        if feature == 'transcript':
                            starts = map(int,re.split(",",features[featId]['starts']))
                            ends = map(int,re.split(",",features[featId]['ends']))
                            exon_count = features[featId]['exon_count']
                            #print(starts)
                            #print(ends)
                            #print(exon_count)
                            #print featId
                            for exonRank in range(exon_count):
                                if(readPosition>=starts[exonRank] and readPosition<=ends[exonRank]):#Check if read is in one of the exons
                                    #print "exon_"+str(exonRank)+" "
                                    features[featId]['counts_sum']  = features[featId]['counts_sum'] + countsforw[readPosition]['count']
                                    features[featId]['counts_n'] += 1
                                    #print readPosition
                                    #print" sum "+ featId+" "+ str(countsforw[readPosition]['count'])+ " added to sum="+str(features[featId]['counts_sum'])+ " count="+str(features[featId]['counts_n'])

                        # Use start/end of feature (i.e. gene or exon) to include count data in summarization
                        else:
                            features[featId]['counts_sum']  = features[featId]['counts_sum'] + countsforw[readPosition]['count']
                            features[featId]['counts_n'] += 1
                            #print readPosition
                            #print" sum "+ featId+" "+ str(countsforw[readPosition]['count'])+ " added to sum="+str(features[featId]['counts_sum'])+ " count="+str(features[featId]['counts_n'])

            #Empty window
            window=[]
            #print 'antisense'
            #This makes 'antisense'
            #Go through all sorted read positions
            for readPosition in sorted(countsrev):
                #print(readPosition)
                for featId in features_rev[:]: #Go through sorted features list. [:] means working on a copy of the features_sorted list, because if you would not, features will move forward in the list if you remove the first feature and in the next iteration, you examine the third feature (which moved to the second place in the list) and you will have skipped the second feature (no being on the first place). Therefore, work on a copy of the original list but remove in the original list. For each new read position however, the copy will be replaced as you start a new iteration over features_sorted
                    #print(featId)
                    if features[featId]['feature_start']<=readPosition:
                        #print " append "+ featId+" "+ str(features[featId]['feature_start'])
                        window.append(featId) #Bring the feature ID over to the window if the read >= feature_start
                        features_rev.remove(featId) #Remove from sorted features so it will not be scanned next time
                    else:
                        break #Because features are sorted, all the other features will be more downstream
                for featId in window[:]:#For all features in the window
                    if features[featId]['feature_end']+int(five_prime_offset)<readPosition:
                        #print" remove "+ featId+" "+ str(features[featId]['feature_end'])
                        window.remove(featId) #Remove features from window which are completely upstream of the read
                    else: #Only features with the read in between start and stop remain in the window

                        # Only features within the exons are counted
                        if feature == 'transcript':
                            starts = map(int,re.split(",",features[featId]['starts']))
                            ends = map(int,re.split(",",features[featId]['ends']))
                            exon_count = features[featId]['exon_count']
                            for exonRank in range(exon_count):
                                if(readPosition>=starts[exonRank] and readPosition<=ends[exonRank]):#Check if read is in one of the exons
                                    #print "exon_"+str(exonRank)+" "
                                    features[featId]['counts_sum']  = features[featId]['counts_sum'] + countsrev[readPosition]['count']
                                    features[featId]['counts_n'] += 1
                                    #print readPosition
                                    #print" sum "+ featId+" "+ str(countsrev[readPosition]['count'])+ " added to sum="+str(features[featId]['counts_sum'])+ " count="+str(features[featId]['counts_n'])

                        # Use start/end of feature (i.e. gene or exon) to include count data in summarization
                        else:
                            features[featId]['counts_sum']  = features[featId]['counts_sum'] + countsrev[readPosition]['count']
                            features[featId]['counts_n'] += 1
                            #print readPosition
                            #print" sum "+ featId+" "+ str(countsrev[readPosition]['count'])+ " added to sum="+str(features[featId]['counts_sum'])+ " count="+str(features[featId]['counts_n'])

            # #Empty window
            window=[]

    except:
        traceback.print_exc()

    return features

##Write exon information and count data to csv file for that chromosome
def write_to_csv(chrom, features, directory, feature, table_name):
    try:
        #Define a csv writer object
        csvWriter = csv.writer(open(directory+"/tmp/"+table_name+"chr"+str(chrom)+"_tmp.csv","wb"))

        #Write all exons and their info to the csv file
        for featId in features:

            if features[featId]['counts_n'] <> 0:
                counts_sum  = features[featId]['counts_sum']
                counts_n  = features[featId]['counts_n']

                csvWriter.writerow([featId, counts_sum, counts_n])

    except:
        traceback.print_exc()

    return

## Dump counts of all chromosomes in SQLite DB ##
def dumpSQLite(chrs, directory, sqliteC,table_name,data):
    try:
        #Make DB connection
        try:
            con = sqlite.connect(sqliteC)
        except:
            print "Could not connect to "+sqliteC
            sys.exit()

        with con:
            cur = con.cursor()
            #Remove possible existing table
            cur.execute("DROP TABLE if exists "+table_name)

            if data in ('rna','ribo'):
                cur.execute("CREATE TABLE "+table_name+ \
                            "(feature_stable_id varchar(128), counts_sum float, counts_n float)")
            if data == 'rrbs':
                cur.execute("CREATE TABLE "+table_name+ \
                            "(feature_stable_id varchar(128), meth_perc_sum float, meth_perc_n float)")


            for chrom in chrs:
                #Convert integer to string (chromosome ids)
                if not isinstance(chrom, basestring):
                    chrom = str(chrom)


                #Store info from chromosomal csv files in SQLite DB
                try:
                    os.system("sqlite3 -separator , "+sqliteC+" \".import "+directory+"/tmp/"+table_name+"chr"+str(chrom)+"_tmp.csv "+table_name+"\"")
                    print "CSV to SQLite dump success for chromosome "+chrom
                except:
                    print "CSV to SQLite dump failed for chromosome "+ chrom

        #Remove tmp files
        for chrom in chrs:
            os.system("rm -rf "+directory+"/tmp/"+table_name+"chr"+str(chrom)+"_tmp.csv")

    except:
        traceback.print_exc()

    return


## Fetch hash from sqlite query ##
# The first column should contain the ID
def fetchDict(dbpathE, query):
    #Init
    output = defaultdict(lambda: defaultdict())

    #Make connection to DB
    try:
        con = sqlite.connect(dbpathE)
        con.text_factory = str
    except:
        print "Could not connect to "+dbpathE
        sys.exit()

    with con:
        cur = con.cursor()
        if cur.execute(query):
            #Fetch all data from query
            query_output = cur.fetchall()
            colnames = list(map(lambda x: x[0], cur.description))

        #Restructure data into a nested dictionary
        for i in range(0, len(query_output)):
            id = query_output[i][0]
            for j in range(1,len(colnames)):
                col = colnames[j]
                value = query_output[i][j]
                output[id][col]=value
    return output


## Data Dumper for dictionaries and defaultDicts ##
def print_dict(dictionary, ident = '', braces=0):
    """ Recursively prints nested dictionaries."""

    for key, value in dictionary.iteritems():
        if isinstance(value, dict):
            print '%s%s%s%s' %(ident,braces*'[',key,braces*']')
            print_dict(value, ident+'  ', braces)
        else:
            print ident+'%s = %s' %(key, value)


## Get coordinate system ID ##
def get_coord_system_id_chr(eDB, assem):
    try:
        con=sqlite.connect(eDB)
    except:
        print "Could not connect to Ensembl DB"
        sys.exit()

    #Init
    csid = 0

    with con:
        cur = con.cursor()

        if cur.execute("SELECT coord_system_id FROM coord_system WHERE name = 'chromosome' AND version = '"+assem+"';"):
            csid = cur.fetchone()[0]
        else:
            print "Could not find coordinate system id"
            sys.exit()

    return csid

## Get arguments ##
def get_arguments(dbpath):
    try:
        con=sqlite.connect(dbpath)
    except:
        print "Could not connect to results DB"
        sys.exit()

    #Init
    igro=''
    sp=''
    vs=0
    cs=0
    rt=''

    with con:
        cur = con.cursor()


        if cur.execute("SELECT value FROM arguments WHERE variable='igenomes_root';"):
            igro = str(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the igenomes root from the arguments table in the SQLite DB"
            sys.exit()

        if cur.execute("SELECT value FROM arguments WHERE variable='species';"):
            sp = str(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the species from the arguments table in the SQLite DB"
            sys.exit()

        if cur.execute("SELECT value FROM arguments WHERE variable='ensembl_version';"):
            vs = int(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the ensembl version from the arguments table in the SQLite DB"
            sys.exit()

        if cur.execute("SELECT value FROM arguments WHERE variable='readtype';"):
            rt = str(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the readtype from the arguments table in the SQLite DB"
            sys.exit()

        if cur.execute("SELECT value FROM arguments WHERE variable='nr_of_cores';"):
            cs = int(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the number of cores from the arguments table in the SQLite DB"
            sys.exit()

    return igro, sp, vs, cs, rt

## Get chromosomes ##
def get_chrs(chrSizeFile, species):
    #Init
    chrs = {}

    #open file
    try:
        FR = open(chrSizeFile, 'r')
    except:
        print "ERROR with opening chromosome sizes file"
    #parse
    for line in FR:
        parts = re.split('\W+',line)
        if species=="fruitfly" and parts[0]=="M":
            chrs["dmel_mitochondrion_genome"] = parts[1]
        else:
            chrs[parts[0]] = parts[1]

    return chrs

## Get chromosomes in XML ##
def get_chrsXml(chrSizeFile, species):
    #Init
    chrs = {}

    #open file
    try:
        FR = open(chrSizeFile, 'r')
    except:
        print "ERROR with opening chromosome sizes xml file"
    #parse
    for line in FR:
        pattern = re.compile('contigName=\"(\w{1,2})\" totalBases\W\"(\d+)\"')
        m = pattern.search(line)
        if m:
            if species=="fruitfly" and m.group(1)=="M":
                chrs["dmel_mitochondrion_genome"] = m.group(2)
            else:
                chrs[m.group(1)] = m.group(2)


    return chrs

#######Set Main##################
if __name__ == "__main__":		#
    try:						#
        main()					#
    except Exception, e:		#
        traceback.print_exc()	#
#################################