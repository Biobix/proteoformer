__author__ = 'Steven Verbruggen'

#############################################################
# Fetches read counts for each ORF in PROTEOFORMER pipeline #
#############################################################

'''

ARGUMENTS:

	-w | --work_dir						 The working directory (default = current directory)
	-s | --sqlitedb						 The sqlite database with data from former steps. (mandatory)
	-t | --tis_ids						 The TIS ID or a comma seperated list of IDs for which the ORF based counts need to be determined. If you set 'all' as input, the ORF based counts of all the IDs will be determined (mandatory)
	-c | --trim							 The first and the last part of the ORF can have seperate counts because of RIBOseq artefacts. These extreme parts can have 'absolute' lengths or 'relative' in function of the ORF length. 'none' will set all pre and post counts to zero (default = relative)
											Attention should be taken for 'absolute' as for too small orfs the will be taken relative because of extreme parts becoming bigger than the actual ORF itself!!
	-n | --nt_trim						 If you chose 'absolute' trim, you have to insert here the length in nt of the extreme parts (<299nt, but suggestions are in the order of magnitude of 15nt). If you chose 'relative', you have to insert here a decimal number (< 0.5). (default=15 for absolute, default=0.025 for relative)
	-l | --ltm							 Whether the program should also make an ORF based counts table of the 'treated' (e.g. LTM) counts: 'y' or 'n' (default 'n')

EXAMPLE:

python ORFbasedCounts.py --sqlitedb SQLite/results.db --tis_ids 1 --trim absolute --nt_trim 15 --ltm n --rna /data2/steven/eIF1/NTmRNA/SQLite/results.db

DEPENDENCIES:
sqlite3


'''

import traceback
import sys
import os
import getopt
import sqlite3 as sqlite
from multiprocessing import Pool
import timeit
import re
from collections import defaultdict
import math
import csv

def main():
	print
	print "				------------------------"
	print "				| ORF BASED RPF COUNTS |"
	print "				------------------------"
	print "\t\tThis program fetches the read counts for each found ORF\n"


	#######################
	###  GETS and SETs  ###
	#######################

	#Timer
	startTime = timeit.default_timer()

	#Init
	workdir=''
	sqlitedb=''

	#Catch command line with getopt
	try:
		myopts, args = getopt.getopt(sys.argv[1:],"w:s:t:c:n:l:r:",["work_dir=","sqlitedb=","tis_ids=","trim=","nt_trim=","ltm=","rna="])
	except getopt.GetoptError as err:
		print err
		sys.exit()

	# Catch arguments
	# o == option
	# a == argument passed to the o
	for o, a in myopts:
		if o in ('-w','--work_dir'):
			workdir=a
		elif o in ('-s','--sqlitedb'):
			sqlitedb=a
		elif o in ('-t','--tis_ids'):
			tis_id=a
		elif o in ('-c','--trim'):
			trim=a
		elif o in ('-n','--nt_trim'):
			nt_trim=a
		elif o in ('-l', '--ltm'):
			ltm=a
		elif o in ('r', '--rna'):
			rna='y'
			rnadb=a

	try:
		workdir
	except:
		workdir=''
	try:
		sqlitedb
	except:
		sqlitedb=''
	try:
		tis_ids
	except:
		tis_ids=''
	try:
		trim
	except:
		trim=''
	try:
		nt_trim
	except:
		nt_trim=''
	try:
		ltm
	except:
		ltm=''
	try:
		rnadb
	except:
		rna='n'
		rnadb=''

	# Check for correct arguments and parse
	if(workdir == ''):
		workdir = os.getcwd()
	if(workdir !=''):
		os.chdir(workdir)
	if(sqlitedb==''):
		print "Error: do not forget to pass the SQLite DB argument!"
		sys.exit()
	if(tis_id==''):
		print "Error: do not forget to pass the TIS IDs argument!"
		sys.exit()
	tis_id_list = []
	if(tis_id=='all'):
		tis_id_list = get_all_tis_ids(sqlitedb)
	elif(tis_id.isdigit()):
		tis_id_list.extend(tis_id)
	elif(bool(re.search(',',tis_id))):
		tis_id_list = re.split(',', tis_id)
	else:
		print "Error: TIS IDs argument should be a number, a comma seperated list of numbers or 'all'!"
		sys.exit()
	if(trim==''):
		trim = 'relative'
	elif(trim!='absolute' and trim!='relative' and trim!='none'):
		print "Error: the trim argument should be 'absolute', 'relative' or 'none'"
		sys.exit()
	if(nt_trim==''):
		if(trim=='absolute'):
			nt_trim=15
		elif(trim=='relative'):
			nt_trim=0.025
		elif(trim=='none'):
			nt_trim=0
	else:
		if(trim=='absolute'):
			try:
				nt_trim=int(nt_trim)
			except:
				print "Error: could not parse the nt_trim argument into an integer!"
				sys.exit()
			if(nt_trim<1):
				print "Error: nt_trim should be an integer starting from 1!"
				sys.exit()
			if(not isinstance(nt_trim, int)):
				print "Error: nt_trim should be an integer starting from 1 for absolute extreme parts!"
				sys.exit()
		elif(trim=='relative'):
			try:
				nt_trim=float(nt_trim)
			except:
				print "Error: could not parse the nt_trim argument into a float!"
				sys.exit()
			if(nt_trim<=0 or nt_trim>=0.5):
				print "Error: nt_trim should be a decimal number inbetween 0 and 0.5 for relative extreme parts!"
				sys.exit()
			if(not isinstance(nt_trim, float)):
				print "Error: nt_trim should be a decimal number inbetween 0 and 0.5 for relative extreme parts!"
				sys.exit()
		elif(trim=='none'):
			if(nt_trim!=0):
				print "nt_trim was set to 0 instead of "+str(nt_trim)+" because 'none' was selected in the trim argument"
				nt_trim=0
		else:
			print "Error: the trim argument should be 'absolute', 'relative' or 'none'"
			sys.exit()
	if (ltm==''):
		ltm='n'
	elif(ltm!='y' and ltm!='n'):
		print "Error: the ltm argument should be 'y' or 'n'"
		sys.exit()

	#Get arguments from arguments table
	ensDB, igenomes_root, species, ens_version, cores = get_arguments(sqlitedb)

	# Print used arguments and parameters
	print "Arguments:"
	print "  Working directory:                         "+workdir
	print "  SQLite DB:                                 "+sqlitedb
	print "  TIS IDs:                                   "+tis_id
	print "  Extreme parts length:                      "+trim
	print "      Threshold:                             "+str(nt_trim)
	print "  LTM:                                       "+ltm
	print "  Ensembl DB:                                "+ensDB
	print "  Ensembl version:                           "+str(ens_version)
	print "  Igenomes root folder:                      "+igenomes_root
	print "  Species:                                   "+species
	print "  Cores:                                     "+str(cores)
	print "  Include RNA:                               "+rna
	if(rna=='y'):
		print "      RNA DB:                                "+rnadb
	print

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
		"BDGP5" if species=="fruitfly" and ens_version<79 else \
		"BDGP6" if species=="fruitfly" and ens_version>=79 else ""
	
	#Make tmp folder if it not exists
	if(not os.path.isdir(workdir+"/tmp")):
		os.system("mkdir "+workdir+"/tmp")
	if(not os.path.isdir(workdir+"/tmp/ORFbasedCounts")):
		os.system("mkdir "+workdir+"/tmp/ORFbasedCounts");

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
	coordSystemId = get_coord_system_id_chr(ensDB, assembly)


	######## FOR TESTING #########
	# chrsY={}
# 	chrsY['Y']=chrs['Y']
# 	chrs={}
# 	chrs=chrsY
	##############################


	##############
	###  MAIN  ###
	##############

	for id in tis_id_list:
		#Get ORF based counts per chr
		print "Determine ORF based counts for each chromosome seperately by starting up multiple processes"
		print
		pool = Pool(processes=cores)
		[pool.apply_async(ORF_based_counts_per_chr, args=(chr,sqlitedb, id, trim, nt_trim, ltm, rna, rnadb, workdir)) for chr in chrs.keys()]
		pool.close()
		pool.join()

		#Dump in SQLite
		print
		print "Dump data into SQLite database"
		print
		dumpSQLite(id, chrs, ltm, rna, workdir, sqlitedb)

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

## Get the ORF based counts per chromosome ##
def ORF_based_counts_per_chr(chr, sqlitedb, tis_id, trim, nt_trim, ltm, rna, rnadb, workdir):
	try:
		#Get the exon structure for each ORF
		orfs, exonBoundaries = get_exon_structure(chr, sqlitedb, tis_id, trim, nt_trim, ltm, rna)

		#Get reads
		CHXfor, CHXrev = getUntreatedReads(chr, sqlitedb)
		if (ltm=='y'):
			LTMfor, LTMrev = getTreatedReads(chr, sqlitedb)
		else:
			LTMfor=''
			LTMrev=''
		if (rna=='y'):
			RNAfor, RNArev = getRNAReads(chr, rnadb)
		else:
			RNAfor=''
			RNArev=''

		#Match the counts to the orfs
		orfs = match_counts_to_orfs(orfs, exonBoundaries, CHXfor, CHXrev, ltm, LTMfor, LTMrev, rna, RNAfor, RNArev)

		#Write to csv file for each chromosome
		write_to_csv(chr, orfs, exonBoundaries, ltm, rna, workdir)

		print "*) Process for chromosome "+chr+" completed"
	except:
		traceback.print_exc()

	return

## Dump counts of all chromosomes in SQLite DB ##
def dumpSQLite(tisid, chrs, ltm, rna, workdir, sqlitedb):
	try:
		#Make DB connection
		try:
			con = sqlite.connect(sqlitedb)
		except:
			print "Could not connect to "+sqlitedb
			sys.exit()

		with con:
			cur = con.cursor()
			tableName = "TIS_"+tisid+"_ORFbasedCounts"
			#Remove possible existing table
			dropQuery = "DROP TABLE IF EXISTS "+tableName
			cur.execute(dropQuery)

			#Create new table for ORF based counts
			if(rna=='n'):
				createQuery = "CREATE TABLE IF NOT EXISTS '"+tableName+"' ("\
					"'ORF_ID' varchar(128) NOT NULL default '',"\
					"'tr_stable_id' varchar(128) NOT NULL default '',"\
					"'start' int(10) NOT NULL default '',"\
					"'chr' char(50) NOT NULL default '',"\
					"'strand' int(2) NOT NULL default '',"\
					"'trim' int(4) NOT NULL default '',"\
					"'start_main_orf' int(10) NOT NULL default '',"\
					"'end_main_orf' int(10) NOT NULL default '',"\
					"'pre_count' float default NULL,"\
					"'count' float default NULL,"\
					"'post_count' float default NULL)"
				cur.execute(createQuery)
			elif(rna=='y'):
				createQuery = "CREATE TABLE IF NOT EXISTS '"+tableName+"' ("\
					"'ORF_ID' varchar(128) NOT NULL default '',"\
					"'tr_stable_id' varchar(128) NOT NULL default '',"\
					"'start' int(10) NOT NULL default '',"\
					"'chr' char(50) NOT NULL default '',"\
					"'strand' int(2) NOT NULL default '',"\
					"'trim' int(4) NOT NULL default '',"\
					"'start_main_orf' int(10) NOT NULL default '',"\
					"'end_main_orf' int(10) NOT NULL default '',"\
					"'pre_count' float default NULL,"\
					"'count' float default NULL,"\
					"'post_count' float default NULL,"\
					"'RNApre_count' float default NULL,"\
					"'RNA_count' float default NULL,"\
					"'RNApost_count' float default NULL)"
				cur.execute(createQuery)


			#Store info from chromosomal csv files in SQLite DB
			for chr in chrs:
				try:
					os.system("sqlite3 -separator , "+sqlitedb+" \".import "+workdir+"/tmp/ORFbasedCounts/ORFbasedCounts_chr"+chr+"_tmp.csv "+tableName+"\"")
				except:
					print "CSV to SQLite dump failed for chromosome "+chr

			#Remove tmp files
			for chr in chrs:
				os.system("rm -rf "+workdir+"/tmp/ORFbasedCounts/ORFbasedCounts_chr"+chr+"_tmp.csv")

			#Same for treated counts
			if(ltm=='y'):
				tableName = "TIS_"+tisid+"_ORFbasedCounts_LTM"
				#Remove possible existing table
				dropQuery = "DROP TABLE IF EXISTS "+tableName
				cur.execute(dropQuery)

				#Create new table for ORF based counts
				if(rna=='n'):
					createQuery = "CREATE TABLE IF NOT EXISTS '"+tableName+"' ("\
						"'ORF_ID' varchar(128) NOT NULL default '',"\
						"'tr_stable_id' varchar(128) NOT NULL default '',"\
						"'start' int(10) NOT NULL default '',"\
						"'chr' char(50) NOT NULL default '',"\
						"'strand' int(2) NOT NULL default '',"\
						"'trim' int(4) NOT NULL default '',"\
						"'start_main_orf' int(10) NOT NULL default '',"\
						"'end_main_orf' int(10) NOT NULL default '',"\
						"'pre_count' float default NULL,"\
						"'count' float default NULL,"\
						"'post_count' float default NULL)"
					cur.execute(createQuery)
				elif(rna=='y'):
					createQuery = "CREATE TABLE IF NOT EXISTS '"+tableName+"' ("\
						"'ORF_ID' varchar(128) NOT NULL default '',"\
						"'tr_stable_id' varchar(128) NOT NULL default '',"\
						"'start' int(10) NOT NULL default '',"\
						"'chr' char(50) NOT NULL default '',"\
						"'strand' int(2) NOT NULL default '',"\
						"'trim' int(4) NOT NULL default '',"\
						"'start_main_orf' int(10) NOT NULL default '',"\
						"'end_main_orf' int(10) NOT NULL default '',"\
						"'pre_count' float default NULL,"\
						"'count' float default NULL,"\
						"'post_count' float default NULL,"\
						"'RNApre_count' float default NULL,"\
						"'RNA_count' float default NULL,"\
						"'RNApost_count' float default NULL)"
					cur.execute(createQuery)

				#Store info from chromosomal csv files in SQLite DB
				for chr in chrs:
					try:
						os.system("sqlite3 -separator , "+sqlitedb+" \".import "+workdir+"/tmp/ORFbasedCounts/ORFbasedCounts_chr"+chr+"_LTM_tmp.csv "+tableName+"\"")
					except:
						print "CSV to SQLite dump failed for chromosome "+chr

				#Remove tmp files
				for chr in chrs:
					os.system("rm -rf "+workdir+"/tmp/ORFbasedCounts/ORFbasedCounts_chr"+chr+"_LTM_tmp.csv")

			os.system("rm -rf "+workdir+"/tmp/ORFbasedCounts")

	except:
		traceback.print_exc()
	return


##Write ORF information and count data to csv file for that chromosome
def write_to_csv(chr, orfs, exonBoundaries, ltm, rna, workdir):
	try:
		#Define a csv writer object
		csvWriter = csv.writer(open(workdir+"/tmp/ORFbasedCounts/ORFbasedCounts_chr"+chr+"_tmp.csv","wb"))

		#Write all ORFs and their info to the csv file
		if(rna=='n'):
			for id in orfs:
				csvWriter.writerow([id, orfs[id]['tr_stable_id'], orfs[id]['start'], chr, orfs[id]['strand'], orfs[id]['trim'], exonBoundaries[id]['start_main_orf'], exonBoundaries[id]['end_main_orf'], orfs[id]['pre_count'], orfs[id]['count'], orfs[id]['post_count']])
		elif(rna=='y'):
			for id in orfs:
				csvWriter.writerow([id, orfs[id]['tr_stable_id'], orfs[id]['start'], chr, orfs[id]['strand'], orfs[id]['trim'], exonBoundaries[id]['start_main_orf'], exonBoundaries[id]['end_main_orf'], orfs[id]['pre_count'], orfs[id]['count'], orfs[id]['post_count'], orfs[id]['RNApre_count'], orfs[id]['RNAcount'], orfs[id]['RNApost_count']])

		if(ltm=='y'):
			#The same for treated counts
			csvWriterLTM = csv.writer(open(workdir+"/tmp/ORFbasedCounts/ORFbasedCounts_chr"+chr+"_LTM_tmp.csv","wb"))

			if(rna=='n'):
				for id in orfs:
					csvWriterLTM.writerow([id, orfs[id]['tr_stable_id'], orfs[id]['start'], chr, orfs[id]['strand'], orfs[id]['trim'], exonBoundaries[id]['start_main_orf'], exonBoundaries[id]['end_main_orf'], orfs[id]['LTMpre_count'], orfs[id]['LTMcount'], orfs[id]['LTMpost_count']])
			elif(rna=='y'):
				for id in orfs:
					csvWriterLTM.writerow([id, orfs[id]['tr_stable_id'], orfs[id]['start'], chr, orfs[id]['strand'], orfs[id]['trim'], exonBoundaries[id]['start_main_orf'], exonBoundaries[id]['end_main_orf'], orfs[id]['LTMpre_count'], orfs[id]['LTMcount'], orfs[id]['LTMpost_count'], orfs[id]['RNApre_count'], orfs[id]['RNAcount'], orfs[id]['RNApost_count']])

	except:
		traceback.print_exc()
	return

## Match the counts to the orfs
def match_counts_to_orfs(orfs, exonBoundaries, CHXfor, CHXrev, ltm, LTMfor, LTMrev, rna, RNAfor, RNArev):
	try:
		#Init
		window = []

		#Make 2 lists of the ORF IDs, sorted by their start position (absolute, so stop for antisense). One for sense ORFs and one for antisense ORFs
		orfs_for=[]
		orfs_rev=[]
		for id in orfs:
			if(orfs[id]['strand']==1):
				orfs_for.append(id)
			elif(orfs[id]['strand']==-1):
				orfs_rev.append(id)
		orfs_for = sorted(orfs_for, key=lambda x: orfs[x]['start'])
		orfs_rev = sorted(orfs_rev, key=lambda x: orfs[x]['stop'])

		#For sense
		#Go through all sorted read positions
		for readPosition in sorted(CHXfor):
			for OrfId in orfs_for[:]: #Go through sorted sense ORFs list. [:] means working on a copy of the orfs_for list, because if you would not, orfs will move forward in the list if you remove the first orf and in the next iteration, you examine the third orf (which moved to the second place in the list) and you will have skipped the second orf (no being on the first place). Therefore, work on a copy of the original list but remove in the original list. For each new read position however, the copy will be replaced as you start a new iteration over orfs_for
				if(orfs[OrfId]['start']<=readPosition):
					window.append(OrfId) #Bring the ORF ID over to the window if the read >= start position ORF
					orfs_for.remove(OrfId) #Remove from sorted sense ORFs so it will not be scanned next time
				else:
					break #Because ORFs are sorted, all the other ORFs will be more downstream
			for OrfId in window:#For all ORFs in the window
				if(orfs[OrfId]['stop']<readPosition):
					window.remove(OrfId) #Remove ORFs from window which are completely upstream of the read
				else: #Only ORFs with the read in between start and stop remain in the window
					for exonRank in range(1,exonBoundaries[OrfId]['number_of_exons']+1):
						if(readPosition>=exonBoundaries[OrfId][exonRank]['start'] and readPosition<=exonBoundaries[OrfId][exonRank]['end']):#Check if read is in one of the exons
							if(readPosition<exonBoundaries[OrfId]['start_main_orf']):#Check if count is in upstream trimmed part of the ORF
								orfs[OrfId]['pre_count'] = orfs[OrfId]['pre_count'] + CHXfor[readPosition]['count']
							elif(readPosition>exonBoundaries[OrfId]['end_main_orf']):#Check if count is in downstream trimmed part of the ORF
								orfs[OrfId]['post_count'] = orfs[OrfId]['post_count'] + CHXfor[readPosition]['count']
							else: #Else, the read is in the main part of the ORF
								orfs[OrfId]['count'] = orfs[OrfId]['count'] + CHXfor[readPosition]['count']

		#Empty window
		window=[]

		#Do similar for antisense
		#For antisense ORFS: start (TIS position) and stop are relative. The starts and stops lists however are absolute and similar to the ones for sense ORFs.
		for readPosition in sorted(CHXrev):
			for OrfId in orfs_rev[:]:
				if(orfs[OrfId]['stop']<=readPosition):
					window.append(OrfId)
					orfs_rev.remove(OrfId)
				else:
					break
			for OrfId in window:
				if(orfs[OrfId]['start']<readPosition):
					window.remove(OrfId)
				else:
					for exonRank in range(1,exonBoundaries[OrfId]['number_of_exons']+1):#Also the rank number are here strand independent (see get_exon_structure)
						if(readPosition>=exonBoundaries[OrfId][exonRank]['start'] and readPosition<=exonBoundaries[OrfId][exonRank]['end']):
							if(readPosition<exonBoundaries[OrfId]['start_main_orf']):
								orfs[OrfId]['post_count'] = orfs[OrfId]['post_count'] + CHXrev[readPosition]['count'] #Start_main_orf and end_main_orf are absolute, so change post and pre count here for antisense
							elif(readPosition>exonBoundaries[OrfId]['end_main_orf']):
								orfs[OrfId]['pre_count'] = orfs[OrfId]['pre_count'] + CHXrev[readPosition]['count']
							else:
								orfs[OrfId]['count'] = orfs[OrfId]['count'] + CHXrev[readPosition]['count']

		#Empty window
		window=[]


		#Then, do the same for the treated counts
		if(ltm=='y'):
			#Make 2 lists of the ORF IDs, sorted by their start position (absolute, so stop for antisense). One for sense ORFs and one for antisense ORFs
			orfs_for=[]
			orfs_rev=[]
			for id in orfs:
				if(orfs[id]['strand']==1):
					orfs_for.append(id)
				elif(orfs[id]['strand']==-1):
					orfs_rev.append(id)
			orfs_for = sorted(orfs_for, key=lambda x: orfs[x]['start'])
			orfs_rev = sorted(orfs_rev, key=lambda x: orfs[x]['stop'])

			#For sense
			#Go through all sorted read positions
			for readPosition in sorted(LTMfor):
				for OrfId in orfs_for[:]: #Go through sorted sense ORFs list. [:] means working on a copy of the orfs_for list, because if you would not, orfs will move forward in the list if you remove the first orf and in the next iteration, you examine the third orf (which moved to the second place in the list) and you will have skipped the second orf (no being on the first place). Therefore, work on a copy of the original list but remove in the original list. For each new read position however, the copy will be replaced as you start a new iteration over orfs_for
					if(orfs[OrfId]['start']<=readPosition):
						window.append(OrfId) #Bring the ORF ID over to the window if the read >= start position ORF
						orfs_for.remove(OrfId) #Remove from sorted sense ORFs so it will not be scanned next time
					else:
						break #Because ORFs are sorted, all the other ORFs will be more downstream
				for OrfId in window:#For all ORFs in the window
					if(orfs[OrfId]['stop']<readPosition):
						window.remove(OrfId) #Remove ORFs from window which are completely upstream of the read
					else: #Only ORFs with the read in between start and stop remain in the window
						for exonRank in range(1,exonBoundaries[OrfId]['number_of_exons']+1):
							if(readPosition>=exonBoundaries[OrfId][exonRank]['start'] and readPosition<=exonBoundaries[OrfId][exonRank]['end']):#Check if read is in one of the exons
								if(readPosition<exonBoundaries[OrfId]['start_main_orf']):#Check if count is in upstream trimmed part of the ORF
									orfs[OrfId]['LTMpre_count'] = orfs[OrfId]['LTMpre_count'] + LTMfor[readPosition]['count']
								elif(readPosition>exonBoundaries[OrfId]['end_main_orf']):#Check if count is in downstream trimmed part of the ORF
									orfs[OrfId]['LTMpost_count'] = orfs[OrfId]['LTMpost_count'] + LTMfor[readPosition]['count']
								else: #Else, the read is in the main part of the ORF
									orfs[OrfId]['LTMcount'] = orfs[OrfId]['LTMcount'] + LTMfor[readPosition]['count']

			#Empty window
			window=[]

			#Do similar for antisense
			#For antisense ORFS: start (TIS position) and stop are relative. The starts and stops lists however are absolute and similar to the ones for sense ORFs.
			for readPosition in sorted(LTMrev):
				for OrfId in orfs_rev[:]:
					if(orfs[OrfId]['stop']<=readPosition):
						window.append(OrfId)
						orfs_rev.remove(OrfId)
					else:
						break
				for OrfId in window:
					if(orfs[OrfId]['start']<readPosition):
						window.remove(OrfId)
					else:
						for exonRank in range(1,exonBoundaries[OrfId]['number_of_exons']+1):#Also the rank number are here strand independent (see get_exon_structure)
							if(readPosition>=exonBoundaries[OrfId][exonRank]['start'] and readPosition<=exonBoundaries[OrfId][exonRank]['end']):
								if(readPosition<exonBoundaries[OrfId]['start_main_orf']):
									orfs[OrfId]['LTMpost_count'] = orfs[OrfId]['LTMpost_count'] + LTMrev[readPosition]['count'] #Start_main_orf and end_main_orf are absolute, so change post and pre count here for antisense
								elif(readPosition>exonBoundaries[OrfId]['end_main_orf']):
									orfs[OrfId]['LTMpre_count'] = orfs[OrfId]['LTMpre_count'] + LTMrev[readPosition]['count']
								else:
									orfs[OrfId]['LTMcount'] = orfs[OrfId]['LTMcount'] + LTMrev[readPosition]['count']

			#Empty window
			window=[]

		#The same for RNA reads but without trimming
		if (rna=='y'):
			#Make 2 lists of the ORF IDs, sorted by their start positon (absolute). One for sense, one for antisense.
			orfs_for=[]
			orfs_rev=[]
			for id in orfs:
				if(orfs[id]['strand']==1):
					orfs_for.append(id)
				elif(orfs[id]['strand']==-1):
					orfs_rev.append(id)
			orfs_for = sorted(orfs_for, key=lambda x: orfs[x]['start'])
			orfs_rev = sorted(orfs_rev, key=lambda x: orfs[x]['stop'])

			#For sense
			#Go through all sorted read positions
			for readPosition in sorted(RNAfor):
				for OrfId in orfs_for[:]:
					if(orfs[OrfId]['start']<=readPosition):
						window.append(OrfId)
						orfs_for.remove(OrfId)
					else:
						break
				for OrfId in window:
					if(orfs[OrfId]['stop']<readPosition):
						window.remove(OrfId)
					else:
						for exonRank in range(1, exonBoundaries[OrfId]['number_of_exons']+1):
							if(readPosition>=exonBoundaries[OrfId][exonRank]['start'] and readPosition<=exonBoundaries[OrfId][exonRank]['end']):
								if(readPosition<exonBoundaries[OrfId]['start_main_orf']):
									orfs[OrfId]['RNApre_count'] = orfs[OrfId]['RNApre_count'] + RNAfor[readPosition]['count']
								elif(readPosition>exonBoundaries[OrfId]['end_main_orf']):
									orfs[OrfId]['RNApost_count'] = orfs[OrfId]['RNApost_count'] + RNAfor[readPosition]['count']
								else:
									orfs[OrfId]['RNAcount'] = orfs[OrfId]['RNAcount'] + RNAfor[readPosition]['count']

			#Empty window
			window=[]


			#Do similar for antisense
			for readPosition in sorted(RNArev):
				for OrfId in orfs_rev[:]:
					if(orfs[OrfId]['stop']<=readPosition):
						window.append(OrfId)
						orfs_rev.remove(OrfId)
					else:
						break
				for OrfId in window:
					if(orfs[OrfId]['start']<readPosition):
						window.remove(OrfId)
					else:
						for exonRank in range(1,exonBoundaries[OrfId]['number_of_exons']+1):
							if(readPosition>=exonBoundaries[OrfId][exonRank]['start'] and readPosition<=exonBoundaries[OrfId][exonRank]['end']):
								if(readPosition<exonBoundaries[OrfId]['start_main_orf']):
									orfs[OrfId]['RNApost_count'] = orfs[OrfId]['RNApost_count'] + RNArev[readPosition]['count']
								elif(readPosition>exonBoundaries[OrfId]['end_main_orf']):
									orfs[OrfId]['RNApre_count'] = orfs[OrfId]['RNApre_count'] + RNArev[readPosition]['count']
								else:
									orfs[OrfId]['RNAcount'] = orfs[OrfId]['RNAcount'] + RNArev[readPosition]['count']

			#Empty window
			window=[]

	except:
		traceback.print_exc()

	return orfs

## Get RNA reads ##
def getRNAReads(chr, dbpath):

	#Fetch data from count tables
	query = "SELECT start, count FROM count_fastq1 WHERE chr='"+chr+"' AND strand = '1';"
	RNAfor = fetchDict(dbpath, query)
	query = "SELECT start, count FROM count_fastq1 WHERE chr='"+chr+"' AND strand = '-1';"
	RNArev = fetchDict(dbpath, query)

	return RNAfor, RNArev

## Get treated reads ##
def getTreatedReads(chr, dbpath):

	#Fetch data from count tables
	query = "SELECT start, count FROM count_fastq2 WHERE chr='"+chr+"' AND strand = '1';"
	LTMfor = fetchDict(dbpath, query)
	query = "SELECT start, count FROM count_fastq2 WHERE chr='"+chr+"' AND strand = '-1';"
	LTMrev = fetchDict(dbpath, query)

	return LTMfor, LTMrev

## Get untreated reads ##
def getUntreatedReads(chr, dbpath):

	#Fetch data from count tables
	query = "SELECT start, count FROM count_fastq1 WHERE chr='"+chr+"' AND strand = '1';"
	CHXfor = fetchDict(dbpath, query)
	query = "SELECT start, count FROM count_fastq1 WHERE chr='"+chr+"' AND strand = '-1';"
	CHXrev = fetchDict(dbpath, query)

	return CHXfor, CHXrev

## Get exon structures of all ORFs ##
def get_exon_structure(chr, dbpath, analysis_id, trim, nt_trim, ltm, rna):
	try:
		#Fetch data of all ORFs from assembly table
		query = "SELECT tr_stable_id||'_'||start as transcript_start, tr_stable_id, strand, start, stop, starts_list, ends_list, aa_seq FROM TIS_"+str(analysis_id)+"_transcripts WHERE chr='"+str(chr)+"';"
		orfs = fetchDict(dbpath, query)
		#Determine the genomic length of the orf. Initialize the read count already
		for orf in orfs:
			orfs[orf]['length']=len(orfs[orf]['aa_seq'])*3
			orfs[orf]['count']=0
			orfs[orf]['pre_count']=0
			orfs[orf]['post_count']=0
			if(rna=='y'):
				orfs[orf]['RNAcount']=0
				orfs[orf]['RNApre_count']=0
				orfs[orf]['RNApost_count']=0
			if(ltm=='y'):
				orfs[orf]['LTMcount']=0
				orfs[orf]['LTMpre_count']=0
				orfs[orf]['LTMpost_count']=0

		#Dict structure: ID -> exon rank -> start/end : value
		#					-> start_main_orf: value
		#					-> end_main_orf: value
		exonBoundaries = defaultdict(lambda: defaultdict(lambda: defaultdict()))
		for id in orfs.keys():
			starts = map(int,re.split("_",orfs[id]['starts_list']))#strand independent
			ends = map(int,re.split("_",orfs[id]['ends_list']))
			exonBoundaries[id]['number_of_exons'] = len(starts)
			for i in range(1, len(starts)+1):# +1 because of defenition of the range function in Python. Rank numbers are defined strand independent here.
				exonBoundaries[id][i]['start']=starts[i-1]
				exonBoundaries[id][i]['end']=ends[i-1]

			#Define start border of the main part of the ORF. Absolute, similar to starts_list and ends_list, not relative like start (TIS) and stop. (Start and stop position are already part of the main ORF)
			if(trim=='absolute'):
				if(orfs[id]['length']<=300):#For absolute, starting from 300nt (sORFs border), keep the trimming relatively constant
					nt_trim_rel = nt_trim/300.0 #Too short: switch to relative
					nt_trim_abs = int(math.floor(nt_trim_rel*orfs[id]['length'])) #Round down to prevent too long trims. Scale in function of ORF length
					orfs[id]['trim']=nt_trim_abs
					rankTrim=1 #To go over exons
					if(nt_trim_abs==0):#Due to round down, it can be that no trimming will be done
						exonBoundaries[id]['start_main_orf']=exonBoundaries[id][rankTrim]['start']
					else:
						while(nt_trim_abs>0):
							if(exonBoundaries[id][rankTrim]['end']<exonBoundaries[id][rankTrim]['start']+nt_trim_abs-1): #To keep exons in mind
								nt_trim_abs = nt_trim_abs - (exonBoundaries[id][rankTrim]['end'] - exonBoundaries[id][rankTrim]['start'] + 1)#Length of trim for the next exon
								rankTrim=rankTrim+1
							elif(exonBoundaries[id][rankTrim]['end']==exonBoundaries[id][rankTrim]['start']+nt_trim_abs-1):#If trim goes to the end of an exon
								exonBoundaries[id]['start_main_orf'] = exonBoundaries[id][rankTrim+1]['start']#Main ORF will start at the start of the next exon
								nt_trim_abs=0#To close while loop
							else:
								exonBoundaries[id]['start_main_orf']=exonBoundaries[id][rankTrim]['start']+nt_trim_abs#Trim stops somewhere in this exon
								nt_trim_abs=0
				else:#For bigger ORFs
					nt_trim_abs = nt_trim
					orfs[id]['trim']=nt_trim_abs
					rankTrim=1
					if(nt_trim_abs==0):
						exonBoundaries[id]['start_main_orf']=exonBoundaries[id][rankTrim]['start']
					else:
						while(nt_trim_abs>0):
							if(exonBoundaries[id][rankTrim]['end']<exonBoundaries[id][rankTrim]['start']+nt_trim_abs-1):
								nt_trim_abs = nt_trim_abs - (exonBoundaries[id][rankTrim]['end'] - exonBoundaries[id][rankTrim]['start'] + 1)
								rankTrim=rankTrim+1
							elif(exonBoundaries[id][rankTrim]['end']==exonBoundaries[id][rankTrim]['start']+nt_trim_abs-1):
								exonBoundaries[id]['start_main_orf'] = exonBoundaries[id][rankTrim+1]['start']
								nt_trim_abs=0
							else:
								exonBoundaries[id]['start_main_orf']=exonBoundaries[id][rankTrim]['start']+nt_trim_abs
								nt_trim_abs=0
			elif(trim=='relative'):
				nt_trim_abs = int(math.floor(nt_trim*orfs[id]['length'])) #Scale in function of ORF length
				orfs[id]['trim'] = nt_trim_abs
				#To keep exons in mind
				rankTrim=1
				if(nt_trim_abs==0):
					exonBoundaries[id]['start_main_orf']=exonBoundaries[id][rankTrim]['start']
				else:
					while(nt_trim_abs>0):
						if(exonBoundaries[id][rankTrim]['end']<exonBoundaries[id][rankTrim]['start']+nt_trim_abs-1):
							nt_trim_abs = nt_trim_abs - (exonBoundaries[id][rankTrim]['end'] - exonBoundaries[id][rankTrim]['start'] + 1)
							rankTrim=rankTrim+1
						elif(exonBoundaries[id][rankTrim]['end']==exonBoundaries[id][rankTrim]['start']+nt_trim_abs-1):
							exonBoundaries[id]['start_main_orf'] = exonBoundaries[id][rankTrim+1]['start']
							nt_trim_abs=0
						else:
							exonBoundaries[id]['start_main_orf']=exonBoundaries[id][rankTrim]['start']+nt_trim_abs
							nt_trim_abs=0
			elif(trim=='none'):
				exonBoundaries[id]['start_main_orf']=exonBoundaries[id][1]['start']#Main ORF simply starts at start of the first exon (=TIS)
				orfs[id]['trim']=0

			#Define end border of the main part, similar
			if(trim=='absolute'):
				if(orfs[id]['length']<=300):
					nt_trim_rel = nt_trim/300.0
					nt_trim_abs = int(math.floor(nt_trim_rel*orfs[id]['length']))
					orfs[id]['trim']=nt_trim_abs
					rankTrim=exonBoundaries[id]['number_of_exons']#get number of last exon
					if(nt_trim_abs==0):
						exonBoundaries[id]['end_main_orf']=exonBoundaries[id][rankTrim]['end']
					else:
						while(nt_trim_abs>0):
							if(exonBoundaries[id][rankTrim]['start']>exonBoundaries[id][rankTrim]['end']-nt_trim_abs+1):
								nt_trim_abs = nt_trim_abs - (exonBoundaries[id][rankTrim]['end'] - exonBoundaries[id][rankTrim]['start'] + 1)
								rankTrim=rankTrim-1
							elif(exonBoundaries[id][rankTrim]['start']==exonBoundaries[id][rankTrim]['end']-nt_trim_abs+1):
								exonBoundaries[id]['end_main_orf'] = exonBoundaries[id][rankTrim-1]['end']
								nt_trim_abs=0
							else:
								exonBoundaries[id]['end_main_orf']=exonBoundaries[id][rankTrim]['end'] - nt_trim_abs
								nt_trim_abs=0
				else:
					nt_trim_abs=nt_trim
					orfs[id]['trim']=nt_trim_abs
					rankTrim=len(starts)#determine last exon: see higher
					if(nt_trim_abs==0):
						exonBoundaries[id]['end_main_orf']=exonBoundaries[id][rankTrim]['end']
					else:
						while(nt_trim_abs>0):
							if(exonBoundaries[id][rankTrim]['start']>exonBoundaries[id][rankTrim]['end']-nt_trim_abs+1):
								nt_trim_abs = nt_trim_abs - (exonBoundaries[id][rankTrim]['end'] - exonBoundaries[id][rankTrim]['start'] + 1)
								rankTrim=rankTrim-1
							elif(exonBoundaries[id][rankTrim]['start']==exonBoundaries[id][rankTrim]['end']-nt_trim_abs+1):
								exonBoundaries[id]['end_main_orf'] = exonBoundaries[id][rankTrim-1]['end']
								nt_trim_abs=0
							else:
								exonBoundaries[id]['end_main_orf']=exonBoundaries[id][rankTrim]['end'] - nt_trim_abs
								nt_trim_abs=0
			elif(trim=='relative'):
				nt_trim_abs = int(math.floor(nt_trim*orfs[id]['length']))
				orfs[id]['trim']=nt_trim_abs
				rankTrim=len(starts)#determine last exon: see higher
				if(nt_trim_abs==0):
					exonBoundaries[id]['end_main_orf']=exonBoundaries[id][rankTrim]['end']
				else:
					while(nt_trim_abs>0):
						if(exonBoundaries[id][rankTrim]['start']>exonBoundaries[id][rankTrim]['end']-nt_trim_abs+1):
							nt_trim_abs = nt_trim_abs - (exonBoundaries[id][rankTrim]['end'] - exonBoundaries[id][rankTrim]['start'] + 1)
							rankTrim=rankTrim-1
						elif(exonBoundaries[id][rankTrim]['start']==exonBoundaries[id][rankTrim]['end']-nt_trim_abs+1):
							exonBoundaries[id]['end_main_orf'] = exonBoundaries[id][rankTrim-1]['end']
							nt_trim_abs=0
						else:
							exonBoundaries[id]['end_main_orf']=exonBoundaries[id][rankTrim]['end'] - nt_trim_abs
							nt_trim_abs=0
			elif(trim=='none'):
				exonBoundaries[id]['end_main_orf']=exonBoundaries[id][len(starts)]['end']
				orfs[id]['trim']=0

	except Exception, e:
		traceback.print_exc()

	return orfs, exonBoundaries

## Get arguments ##
def get_arguments(dbpath):
	try:
		con=sqlite.connect(dbpath)
	except:
		print "Could not connect to results DB"
		sys.exit()

	#Init
	ensDB=''
	igro=''
	sp=''
	vs=0
	cs=0

	with con:
		cur = con.cursor()

		if cur.execute("SELECT value FROM arguments WHERE variable='ens_db';"):
			eDB = str(cur.fetchone()[0])
		else:
			print "ERROR: could not fetch Ensembl DB from arguments table in the SQLite DB"
			sys.exit()

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

		if cur.execute("SELECT value FROM arguments WHERE variable='nr_of_cores';"):
			cs = int(cur.fetchone()[0])
		else:
			print "ERROR: could not fetch the number of cores from the arguments table in the SQLite DB"
			sys.exit()

	return eDB, igro, sp, vs, cs

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

## Get all TIS IDs from TIS overview table ##
def get_all_tis_ids(dbpath):

	try:
		con=sqlite.connect(dbpath)
	except:
		print "Error: could not connect to results DB"
		sys.exit()

	with con:
		cur=con.cursor()

		if cur.execute("SELECT ID FROM TIS_overview"):
			tis_ids_list = [ item[0] for item in cur.fetchall()]

	return tis_ids_list

## Data Dumper for dictionaries and defaultDicts ##
def print_dict(dictionary, ident = '', braces=0):
	""" Recursively prints nested dictionaries."""

	for key, value in dictionary.iteritems():
		if isinstance(value, dict):
			print '%s%s%s%s' %(ident,braces*'[',key,braces*']')
			print_dict(value, ident+'  ', braces)
		else:
			print ident+'%s = %s' %(key, value)

## Fetch hash from sqlite query ##
# The first column should contain the ID
def fetchDict(dbpath, query):
	#Init
	output = defaultdict(lambda: defaultdict())

	#Make connection to DB
	try:
		con = sqlite.connect(dbpath)
		con.text_factory = str
	except:
		print "Could not connect to "+dbpath
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


#######Set Main##################
if __name__ == "__main__":		#
	try:						#
		main()					#
	except Exception, e:		#
		traceback.print_exc()	#
#################################
