import traceback
import sys
import os
import argparse
from collections import defaultdict
import re
import sqlite3 as sqlite


def main():


    #Parse arguments from command line
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS,
                                     description="This tool is part of the PROTEOFORMER pipeline. It converts a combined"
                                                 "fasta file of different TIS ids to an SQLite DB for easier consultation.",
                                     add_help=False)

    man_args = parser.add_argument_group("Mandatory parameters")
    man_args.add_argument("--fasta", "-f", action="store", required=True, nargs="?", metavar="PATH",
                          type=str, help="The combined fasta file (mandatory)")

    opt_args = parser.add_argument_group("Optional parameters")
    opt_args.add_argument("--help", "-h", action="help", help="Show help message and exit")
    opt_args.add_argument("--workdir", "-w", action="store", required=False, nargs="?", metavar="FOLDER", default="",
                          type=str, help="Working directory (default: CWD)")
    opt_args.add_argument("--tmp_csv", "-t", action="store", required=False, nargs="?", metavar="PATH",
                          default="combfasta.csv", type=str,
                          help="Temporary csv file for storing table structure "
                               "(default: combfasta.csv)")
    opt_args.add_argument("--output", "-o", action="store", required=False, nargs="?", metavar="PATH",
                          default="combfasta.db", type=str,
                          help="Ouptut SQLite database (default: combfasta.db)")

    args = parser.parse_args()

    #default workdir is CWD, default tmp, default fdr
    if args.workdir == "":
        args.workdir = os.getcwd()

    ########
    # MAIN #
    ########


    input_data = read_proteoformer_fasta(args.fasta)
    parsed_data = parse_data(input_data)
    store_in_db(parsed_data, args.tmp_csv, args.output)

    return

##########
#  SUBS  #
##########

#Store in output
def store_in_db(parsed_data, tmp_csv, db):

    #Connect to sqlite db
    table_name = "combfasta"
    con = sqlite.connect(db)

    with con:
        cur = con.cursor()

        #Delete existing table
        drop_query = "DROP TABLE IF EXISTS '"+table_name+"';"
        cur.execute(drop_query)

        #Create new table
        create_query = "CREATE TABLE IF NOT EXISTS '"+table_name+"' (" \
                        "'tr_stable_id' varchar(128) NOT NULL default ''," \
                        "'chr' char(50) NOT NULL default ''," \
                        "'start' int(10) NOT NULL default ''," \
                        "'annotation' varchar(128) NOT NULL default ''," \
                        "'bincode' varchar(128) NOT NULL default '',"\
                        "'main_db' int(10) NOT NULL default '',"\
                        "'gene_stable_id' varchar(128) NOT NULL default ''," \
                        "'start_codon' varchar(128) NOT NULL default ''," \
                        "'aTIS_call' varchar(128) NOT NULL default ''," \
                        "'side_accessions' TEXT NOT NULL default ''," \
                        "'sequence' TEXT NOT NULL default '');"
        cur.execute(create_query)

        #Construct csv
        with open(tmp_csv, 'w') as FW:
            for id in parsed_data:
                line = parsed_data[id]['tr_stable_id']+","+ \
                       parsed_data[id]['chr'] + "," + \
                       parsed_data[id]['start'] + "," + \
                       parsed_data[id]['annotation'] + "," + \
                       parsed_data[id]['bincode'] + "," + \
                       parsed_data[id]['main_db'] + "," + \
                       parsed_data[id]['gene_stable_id'] + "," + \
                       parsed_data[id]['start_codon'] + "," + \
                       parsed_data[id]['aTIS_call'] + "," + \
                       parsed_data[id]['side_accessions'] + "," + \
                       parsed_data[id]['sequence'] +"\n"
                FW.write(line)

        # Dump
        try:
            os.system("sqlite3 -separator , " + db + " \".import " + tmp_csv + " " + table_name + "\"")
            os.system("rm -rf "+tmp_csv)
        except:
            print "ERROR: dumping of file " + tmp_csv + " into database " + db + " failed!"
            sys.exit()


    return

#Parse data
def parse_data(input_data):

    #Init
    parsed_data = defaultdict(lambda: defaultdict())

    id=1
    for accession in input_data:
        m1 = re.search('^>generic\|(ENST\d+?)\_(\S+?)\_(\d+?)\_(\S+?)\_(\d+?)db(\d+?)\|(ENSG\d+?) (\S+?) (\S+?) ', accession)
        if m1:
            parsed_data[id]['tr_stable_id'] = m1.group(1)
            parsed_data[id]['chr'] = m1.group(2)
            parsed_data[id]['start'] = m1.group(3)
            parsed_data[id]['annotation'] = m1.group(4)
            parsed_data[id]['bincode'] =  m1.group(5)
            parsed_data[id]['main_db'] = m1.group(6)
            parsed_data[id]['gene_stable_id'] = m1.group(7)
            parsed_data[id]['start_codon'] = m1.group(8)
            parsed_data[id]['aTIS_call'] = m1.group(9)
            parsed_data[id]['sequence'] = input_data[accession]
            m2 = re.search('\[(.+)\]$', accession)
            if m2:
                parsed_data[id]['side_accessions'] = m2.group(1)
            else:
                parsed_data[id]['side_accessions'] = ""
            id+=1
        else:
            print accession


    return parsed_data


#read proteoformer fasta
def read_proteoformer_fasta(fasta):

    #Init
    input_data = defaultdict(lambda: defaultdict())

    #Open fasta file
    with open(fasta, 'r') as FR:
        lines = FR.readlines()
        lines = map(lambda x: x.rstrip("\n"), lines)
        lines = map(lambda x: x.rstrip("\r"), lines)
        #Init
        accession = ""
        sequence = ""
        for i in range(0, len(lines)):
            # Check if accession -> new entry
            if ((re.search('^>', lines[i]))):
                if(i!=0):
                    # Save previous entry
                    input_data[accession] = sequence
                # Get new entry accession and empty sequence
                accession = lines[i]
                sequence = ""
            else:
                # Concat next line of sequence
                sequence += lines[i]
        # Save the last entry
        if (sequence != ""):
            input_data[accession] = sequence
            accession = ""
            sequence = ""


    return input_data



###### DIRECT TO MAIN ############
if __name__ == "__main__":
    try:
        main()
    except Exception, e:
        traceback.print_exc()
##################################