import traceback
import getopt
from collections import defaultdict
import sqlite3
import os
import sys
import pandas as pd
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import hex2color
from matplotlib import cm
import seaborn as sns
from matplotlib.ticker import ScalarFormatter
import re
import time



__author__ = 'Steven Verbruggen'

'''

Make mappingQC plots

ARGUMENTS

    -w | --work_dir                         The working directory
                                                (default = current working directory)
    -r | --result_db                        The results db with data to plot
                                                (mandatory)
    -s | --input_samfile                    The sam file used as input
                                                (default STAR/fastq(1/2)/(un)treat.sam)
    -t | --treated                          Untreated or treated sample
                                                (default: untreated) (should be untreated/treated)
    -m | --mapping_unique                   Whether mapping was executed unique or not (Y/N)
                                                (default: Y)
    -f | --firstRankMultiMap                Whether mapping was executed with first rank multimapping (Y/N)
                                                (default: Y)
    -u | --unique                           MappingQC unique (Y/N)
                                                (default: Y)
    -x | --plotrpftool                      The tool to plot the RPF phase figure (grouped2D/pyplot3D/mayavi)
                                                (default: grouped2D)
    -o | --outfolder                        The output folder for plots
                                                (default images_plots)
    -h | --outhtml                          The output html file
                                                (default mappingQC_out_(un)treated.html)
    -z | --outzip                           The output zip file name
                                                (default mappingQC_(un)treated.zip)
    -p | --plastid_option                   Origin of offsets (plastid, standard or from_file)
                                                (default: plastid)
    -i | --plastid_img                      Path to the plastid offset image
                                                (mandatory if plastid option equals 'plastid')
    -e | --ensembl_db                       Ensembl database
                                                (mandatory)

EXAMPLE

python mappingQC.py --result_db SQLite/results.db

'''

def main():

    # Catch command line with getopt
    try:
        myopts, args = getopt.getopt(sys.argv[1:], "w:r:s:t:m:f:u:x:o:h:z:p:i:e:", ["work_dir=", "result_db=", "input_samfile=", "treated=", "mapping_unique=", "firstRankMultiMap=", "unique=", "plotrpftool=", "outfile=", "outhtml=", "outzip=", "plastid_option=", "plastid_img=" ,"ensembl_db="])
    except getopt.GetoptError as err:
        print err
        sys.exit()

    # Catch arguments
    # o == option
    # a == argument passed to the o
    for o, a in myopts:
        if o in ('-w', '--work_dir'):
            workdir = a
        if o in ('-r', '--result_db'):
            result_db = a
        if o in ('-s', '--input_samfile'):
            samfile = a
        if o in ('-t', '--treated'):
            treated = a
        if o in ('-m', '--mapping_unique'):
            mapping_unique = a
        if o in ('-f', '--firstRankMultiMap'):
            firstRankMultiMap = a
        if o in ('-u', '--unique'):
            unique = a
        if o in ('-x', '--plotrpftool'):
            plotrpftool = a
        if o in ('-o', '--outfolder'):
            outfolder = a
        if o in ('-h', '--outhtml'):
            outhtml = a
        if o in ('-z', '--outzip'):
            outzip = a
        if o in ('-p', '--plastid_option'):
            plastid_option = a
        if o in ('-i', '--plastid_img'):
            plastid_img = a
        if o in ('-e', '--ensembl_db'):
            ens_db = a

    try:
        workdir
    except:
        workdir = ''
    try:
        result_db
    except:
        result_db = ''
    try:
        samfile
    except:
        samfile = ''
    try:
        treated
    except:
        treated=''
    try:
        mapping_unique
    except:
        mapping_unique = ''
    try:
        firstRankMultiMap
    except:
        firstRankMultiMap = ''
    try:
        unique
    except:
        unique = ''
    try:
        plotrpftool
    except:
        plotrpftool = ''
    try:
        outfolder
    except:
        outfolder = ''
    try:
        outhtml
    except:
        outhtml = ''
    try:
        outzip
    except:
        outzip = ''
    try:
        plastid_option
    except:
        plastid_option = ''
    try:
        plastid_img
    except:
        plastid_img = ''
    try:
        ens_db
    except:
        ens_db = ''

    # Check for correct arguments and aprse
    if workdir == '':
        workdir = os.getcwd()
    if workdir != '':
        os.chdir(workdir)
    tmpfolder = workdir+"/tmp"
    if result_db == '':
        print "ERROR: do not forget to mention the result DB!"
        sys.exit()
    if samfile == '':
        if treated=="untreated":
            samfile="STAR/fastq1/untreat.sam"
        else:
            samfile="STAR/fastq2/treat.sam"
    if treated == '':
        treated="untreat"
    elif treated != "untreated" and treated != "treated":
        print "ERROR: treated option should be 'treated' or 'untreated'!"
        sys.exit()
    if mapping_unique == '':
        mapping_unique = 'Y'
    if firstRankMultiMap == '':
        firstRankMultiMap = 'Y'
    if unique == '':
        unique = 'Y'
    if plotrpftool == '':
        plotrpftool = "grouped2D"
    if outfolder == '':
        outfolder = workdir+"/images_plots/"
    if outhtml == '':
        outhtml = outfolder+"/mappingQC.html"
        outhtml_short = "mappingQC.html"
    else:  # only take the last part of the name
        backslash_test = re.search('/', outhtml)
        if backslash_test:
            m = re.search('.*/(.+?\.(html|dat))$', outhtml)
            if m:
                outhtml_short = m.group(1)
            else:
                print "Could not extract html file name out of given path ("+outhtml+")"
                sys.exit()
        else:
            m = re.search('(.+?\.(html|dat))$', outhtml)
            if m:
                outhtml_short = m.group(1)
            else:
                print "Could not extract html file name out of given path ("+outhtml+")"
                sys.exit()
    if outzip == '':
        outzip = workdir+"/mappingQC.zip"
        outzip_short = "mappingQC.zip"
    else: # Only take the last part of the name
        backslash_test = re.search('/', outzip)
        if backslash_test:
            m = re.search('.*/(.+?\.(zip|dat))$', outzip)
            if m:
                outzip_short = m.group(1)
            else:
                print "Could not extract zip file name out of given path ("+outzip+")"
        else:
            m = re.search('(.+?\.(zip|dat))$', outzip)
            if m:
                outzip_short = m.group(1)
            else:
                print "Could not extract zip file name out of given path (" + outzip + ")"
    if plastid_option == '':
        plastid_option='plastid'
    elif plastid_option!='standard' and plastid_option!='plastid' and plastid_option!='from_file':
        print "ERROR: plastid option should be 'plastid', 'standard' or 'from_file'!"
        sys.exit()
    if plastid_option == 'plastid':
        if plastid_img == '':
            print "ERROR: do not forget to give path to plastid image if offset option equals 'plastid'!"
    if ens_db=='':
        print "ERROR: do not forget to mention the Ensembl db!"
        sys.exit()

    species = get_species(result_db)

    ########
    # MAIN #
    ########

    #Make plot  directory
    if not os.path.exists(outfolder):
        os.system("mkdir -p "+outfolder)

    #Download biobix image and mappingqc images
    os.system("wget --quiet \"http://galaxy.ugent.be/static/BIOBIX_logo.png\"")
    os.system("mv BIOBIX_logo.png "+outfolder)
    os.system("wget --no-check-certificate --quiet \"https://github.com/Biobix/mQC/raw/master/mqc_tools/logo_mqc2_whitebg.png\"")
    os.system("mv logo_mqc2_whitebg.png " + outfolder)

    #Download codon refs for human or mouse
    codon_ref_file = ""
    if species == 'human':
        os.system("wget --quiet --no-check-certificate \"https://raw.githubusercontent.com/Biobix/mQC/master/mqc_tools/codon_refs/codon_reference_human.csv\"")
        os.system("mv codon_reference_human.csv "+tmpfolder+"/mappingqc_"+treated)
        codon_ref_file = tmpfolder+"/mappingqc_"+treated+"/codon_reference_human.csv"
    if species == 'mouse':
        os.system("wget --quiet --no-check-certificate \"https://raw.githubusercontent.com/Biobix/mQC/master/mqc_tools/codon_refs/codon_reference_mouse.csv\"")
        os.system("mv codon_reference_mouse.csv " + tmpfolder + "/mappingqc_"+treated)
        codon_ref_file = tmpfolder + "/mappingqc_"+treated+"/codon_reference_mouse.csv"

    #Get plot data out of results DB
    phase_distr, total_phase_distr, triplet_distr, run_name, totmaps, prefix_gene_distr, ensembl_version, species = get_plot_data(result_db, treated)

    #Make total phase distribution plot
    outfile = outfolder+"/tot_phase.png"
    plot_total_phase(total_phase_distr, outfile)

    #Make RPF-phase distribution plot
    outfile = outfolder+"/rpf_phase.png"
    if plotrpftool == "grouped2D":
        plot_rpf_phase_grouped2D(phase_distr, outfile)
    elif plotrpftool == "pyplot3D":
        plot_rpf_phase_pyplot3D(phase_distr, outfile)
    elif plotrpftool == "mayavi":
        plot_rpf_phase_mayavi(phase_distr, outfile)

    #Make phase position distribution
    phase_position_distr(tmpfolder, outfolder, treated)

    #Make triplet identity plots
    triplet_plots(triplet_distr, outfolder)

    #Make codon usage plot
    if species=='human' or species=='mouse':
        codon_usage_plot(tmpfolder, codon_ref_file, outfolder, run_name, treated, triplet_distr)

    #Write to output html file
    offsets_file = tmpfolder+"/mappingqc_"+treated+"/mappingqc_offsets_"+treated+".csv"
    if(treated=="untreated"):
        prefix_gene_distr += "_untreated_"
    else:
        prefix_gene_distr += "_treated_"
    # Copy offsets image to output folder
    offset_img = "offsets.png"
    if plastid_option=="plastid":
        os.system("cp " + plastid_img + " " + outfolder + "/" + offset_img)
    # Copy metagenic pie charts to output folder
    tmp_rankedgenes = tmpfolder + "/mappingqc_"+treated+"/rankedgenes.png"
    tmp_cumulative = tmpfolder + "/mappingqc_"+treated+"/cumulative.png"
    tmp_density = tmpfolder + "/mappingqc_"+treated+"/density.png"
    tmp_metagenic_plot_c = tmpfolder + "/mappingqc_"+treated+"/annotation_coding.png"
    tmp_metagenic_plot_nc = tmpfolder + "/mappingqc_"+treated+"/annotation_noncoding.png"
    os.system("cp " + tmp_cumulative + " " + outfolder)
    os.system("cp " + tmp_rankedgenes + " " + outfolder)
    os.system("cp " + tmp_density + " " + outfolder)
    os.system("cp " + tmp_metagenic_plot_c + " " + outfolder)
    os.system("cp " + tmp_metagenic_plot_nc + " " + outfolder)
    write_out_html(outhtml, samfile, run_name, totmaps, plastid_option, offsets_file, offset_img, prefix_gene_distr, ensembl_version, species, ens_db, treated, mapping_unique, firstRankMultiMap, unique)

    ##Archive and collect output
    #Make output archive
    output_arch = "mappingQC_archive/"
    #Try to fetch alternative name out of zip file name
    m = re.search('(.+)\.zip$', outzip_short)
    if m:
        output_arch = m.group(1)
    os.system("mkdir " + workdir + "/" + output_arch)
    #Bring output html and output images folder to archive
    os.system("cp -r " + outhtml + " " + outfolder)
    os.system("cp -r " + outfolder + " " + output_arch)
    #zip output archive
    tmpZip = workdir+"/tmp.zip"
    os.system("zip -r -q "+tmpZip+" "+output_arch)
    os.system("rm -rf " + output_arch)
    os.system("mv "+tmpZip+" "+outzip)

    return

############
### SUBS ###
############


## Write output html file
def write_out_html(outfile, samfile, run_name, totmaps, plastid, offsets_file, offsets_img, prefix_gene_distr, ensembl_version, species, ens_db, treated, mapping_unique, firstRankMultiMap, unique):

    #Load in offsets
    offsets = pd.read_csv(offsets_file, sep=',', header=None, names=["RPF", "offset"])
    max_rpf = offsets["RPF"].max()
    min_rpf = offsets["RPF"].min()
    html_table=""
    for ofs in range(min_rpf, max_rpf+1, 1):
        html_table += """<tr>
        <td>"""+str(ofs)+"""</td>
        <td>"""+str(int(offsets.loc[offsets["RPF"] == ofs]["offset"]))+"""</td>
        </tr>
        """

    #Prepare additional pieces of HTML code for eventual plastid analysis
    plastid_nav_html=""
    plastid_html=""
    if(plastid=="plastid"):
        plastid_nav_html = "<li><a href=\"#section2\">Plastid offset analysis</a></li>"
        plastid_html = """<span class="anchor" id="section2"></span>
        <h2 id="plastid">Plastid offset analysis</h2>
        <p>
        <table id="offset_table">
            <tr>
                <th id="table_header">RPF length</th>
                <th id="table_header">Offset</th>
            </tr>
            """+html_table+"""
        </table>
        <div class="img" id="plastid_img">
            <img src=\""""+offsets_img+"""\" alt="Plastid analysis" id="plastid_plot">
        </div>
        </p>
        """
    else:
        plastid_nav_html = "<li><a href=\"#section2\">Offsets overview</a></li>"
        plastid_html = """<span class="anchor" id="section2"></span>
                <h2 id="plastid">Offsets overview</h2>
                <p>
                <table id="offset_table">
                    <tr>
                        <th id="table_header">RPF length</th>
                        <th id="table_header">Offset</th>
                    </tr>
                    """ + html_table + """
                </table>
                </p>
                """

    #Prepare additional pieces for codon usage plot
    codon_usage_nav = ""
    codon_usage_part = ""
    if species=='human' or species=='mouse':
        codon_usage_nav = """<li><a href="#section9">Codon usage plot</a></li>\n"""
        codon_usage_part = """
        <span class="anchor" id="section9"></span>
        <h2 id="codon_usage">Codon usage plot</h2>
        <p>
            <div class="img">
            <img src=\"codon_usage.png" alt="codon_usage_plot" id="codon_usage_img">
            </div>
        </p>
        """

    #Structure of html file
    html_string = """<!DOCTYPE html>
<html>
<head>
   <title>Mapping QC Report """+run_name+"""</title>
   <meta charset="utf-8"></meta>
   <meta name="description" content="Overview HTML of all mappingQC results"></meta>
   <link href="https://fonts.googleapis.com/css?family=Indie+Flower" rel="stylesheet">
   <style media="screen">
        *{
            box-sizing: border-box;
            font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
        }


        nav{
            float:left;
            padding: 15px;
            width: 17%;
            position: fixed;
            height: 100%;
            overflow: auto;
            border-right: ridge;
            border-color: lightgrey;
            margin-top: 60px;
            margin-left: -20px;
            padding-left: 20px;
            background-color: white;
            z-index:1;
        }

        nav ul {
            list-style-type: none;
            margin: 0px;
            padding: 5px;
            padding-top: 15px;

        }

        nav li{
            padding: 8px;
            margin-bottom: 8px;
            background-color: #33b5e5;
            color: #ffffff;
        }

        nav li:hover {
            background-color: #0099cc;
        }

        #content{
            position: absolute;
            margin-left:19%;
            height: 76%;
        }

        #rpf_phase{
            width:100%;
        }

        #header{
            background-color: grey;
            color: white;
            position:fixed;
            height: 2.7cm;
            width:110%;
            padding: 15px;
            padding-top: 10px;
            margin-left: -10px;
            margin-top: -30px;
            margin-right: -10px;
            margin-bottom: 10px;
            overflow: visible;
            z-index: 2;
        }

        #mappingqc{
            font-family: 'Indie Flower', cursive;
            font-size: 44px;
            padding-left: 120px;
            position: relative;
            z-index: 4;
        }
        #run_name{
            position:relative;
            z-index: 4;
            display:table;
            margin:0 auto;
            margin-top:-50px;
        }
        #mqc_logo{
            height: 60%;
            position: absolute;
            left: 20px;
            top: 30px;
        }
        #biobix_logo{
            height:60%;
            position: absolute;
            right: 200px;
            top: 30px;
        }

        a {
            color: inherit;
            text-decoration: none;
        }

        .anchor{
            display: block;
            height: 14%; /*same height as header*/
            margin-top: -10px; /*same height as header*/
            visibility: hidden;
        }

        #analysis_info_table{
            border-style: none;
            border-width: 0px;
        }

        th {
            border-style: solid;
            border-width: 0px;
            border-color: white;
            border-collapse: collapse;
            padding: 5px;
            background-color: #33b5e5;
            color: #ffffff;
        }

        td {
            border-style: solid;
            border-width: 0px;
            border-color: white;
            border-collapse: collapse;
            background-color: #f2f2f2;
            padding: 5px;
        }

        img {
          max-width: 98%;
          height: auto;
          width: auto\9; /* ie8 */
        }

        #phase_relpos_distr_img {
          width:98%
        }

        #ranked_genes, #cumulative, #genes_density, #annotation_coding, #annotation_noncoding {
            width: 20cm;
        }

        #offset_table {
            float: left;
            display: block;
            margin-right: 120px;
        }

        #plastid_img {
            float: left;
            display: block;
            max-width: 600px;
        }

        #section3 {
            clear: left;
        }


        #footer{
            background-color: grey;
            color: white;
            position: fixed;
            bottom: 0cm;
            padding-left: 30px;
            margin-left: -30px;
            height: 0.7cm;
            width: 110%;
            z-index: 2;
        }
        #footer_content{
            position: fixed;
            bottom: -0.3cm;
        }
   </style>

   <style media="print">
       *{
           box-sizing: border-box;
           font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
       }


       nav{
           visibility: hidden;
           float:left;
           padding: 15px;
           width: 17%;
           position: fixed;
           height: 100%;
           overflow: auto;
           border-right: ridge;
           border-color: lightgrey;
           margin-top: 60px;
           margin-left: -20px;
           padding-left: 20px;
           background-color: white;
           z-index:1;
       }

       nav ul {
           list-style-type: none;
           margin: 0px;
           padding: 5px;
           padding-top: 15px;
    }

       nav li{
           padding: 8px;
           margin-bottom: 8px;
           background-color: #33b5e5;
           color: #ffffff;
       }

       nav li:hover {
           background-color: #0099cc;
       }

       #content{
           margin-top: 70px;
           position: absolute;
           margin-left:0%;
           height: 76%;
       }

       #rpf_phase{
           width:100%;
       }

       #header{
           visibility: hidden;
           background-color: grey;
           color: white;
           position:fixed;
           height: 2.7cm;
           width:110%;
           padding: 15px;
           padding-top: 10px;
           margin-left: -10px;
           margin-top: -30px;
           margin-right: -10px;
           margin-bottom: 10px;
           overflow: visible;
           z-index: 2;
       }

       #mappingqc{
           font-family: 'Indie Flower', cursive;
           font-size: 44px;
           padding-left: 10px;
           position: relative;
           z-index: 4;
       }
       #run_name{
           padding-left: 43px;
           position: relative;
           z-index: 4;
       }

       #biobix_logo{
           height:60%;
           position: absolute;
           right: 200px;
           top: 30px;
       }

       a {
           color: inherit;
           text-decoration: none;
       }

       .anchor{
           display: block;
           height: 14%; /*same height as header*/
           margin-top: -10px; /*same height as header*/
           visibility: hidden;
       }

       #analysis_info_table{
           border-style: none;
           border-width: 0px;
       }

       th {
           border-style: solid;
           border-width: 0px;
           border-color: white;
           border-collapse: collapse;
           padding: 5px;
           background-color: #33b5e5;
           color: #ffffff;
       }

   td {
       border-style: solid;
       border-width: 0px;
       border-color: white;
       border-collapse: collapse;
       background-color: #f2f2f2;
       padding: 5px;
   }

   img {
       max-width: 98%;
       height: auto;
       width: auto\9; /* ie8 */
   }

       #phase_relpos_distr_img {
          width:98%
       }

       #ranked_genes, #cumulative, #genes_density, #annotation_coding, #annotation_noncoding {
           width: 20cm;
       }

       #offset_table {
           float: left;
           display: block;
           margin-right: 120px;
       }

       #plastid_img {
           float: left;
           display: block;
           max-width: 600px;
       }

       #section3 {
           clear: left;
       }


       #footer{
           visibility: hidden;
           background-color: grey;
           color: white;
           position: fixed;
           bottom: 0cm;
           padding-left: 30px;
           margin-left: -30px;
           height: 0.7cm;
           width: 110%;
           z-index: 2;
       }
       #footer_content{
           position: fixed;
           bottom: -0.3cm;
       }
   </style>

</head>

<body>
    <div id="header">
        <img src="logo_mqc2_whitebg.png" alt="mqc_logo" id="mqc_logo">
        <h1><span id="mappingqc">Mapping QC</span><span id="run_name">"""+run_name+"""</span></h1>
        <img src=\"BIOBIX_logo.png\" alt="biobix_logo" id="biobix_logo">
    </div>

    <nav id="navigator">
        <ul>
            <li><a href="#section1">Analysis information</a></li>
            """+plastid_nav_html+"""
            <li><a href="#section3">Gene distributions</a></li>
            <li><a href="#section4">Metagenic classification</a></li>
            <li><a href="#section5">Total phase distribution</a></li>
            <li><a href="#section6">RPF phase distribution</a></li>
            <li><a href="#section7">Phase - relative position distribution</a></li>
            <li><a href="#section8">Triplet identity plots</a></li>
            """+codon_usage_nav+"""
        </ul>
    </nav>

    <div id="content">
        <span class="anchor" id="section1"></span>
        <h2 id="info">Analysis information</h2>
        <p>
        <table id="analysis_info_table">
            <tr>
                <th id="table_header">Feature</th>
                <th id="table_header">Value</th>
            </tr>
            <tr>
                <td>Species</td>
                <td>"""+species+"""</td>
            </tr>
            <tr>
                <td>Input sam file</td>
                <td>"""+samfile+"""</td>
            </tr>
            <tr>
                <td>Ensembl version</td>
                <td>"""+ensembl_version+"""</td>
            </tr>
            <tr>
                <td>Ensembl database</td>
                <td>"""+ens_db+"""</td>
            </tr>
            <tr>
                <td>Sample treatment</td>
                <td>"""+treated+"""</td>
            </tr>
            <tr>
                <td>Mapping unique?</td>
                <td>"""+mapping_unique+"""</td>
            </tr>
            <tr>
                <td>Mapping first rank?</td>
                <td>"""+firstRankMultiMap+"""</td>
            </tr>
            <tr>
                <td>MappingQC unique?</td>
                <td>"""+unique+"""</td>
            </tr>
            <tr>
                <td>Selected offset source</td>
                <td>"""+plastid+"""</td>
            </tr>
            <tr>
                <td>Mapped genomic sequences</td>
                <td>"""+'{0:,}'.format(totmaps).replace(',',' ')+"""</td>
            </tr>
            <tr>
                <td>Analysis date</td>
                <td>"""+time.strftime("%A %d %b %Y")+"""</td>
            </tr>
            <tr>
                <td>Analysis time</td>
                <td>"""+time.strftime("%H:%M:%S")+"""</td>
            </tr>
        </table>
        </p>

        """+plastid_html+"""

        <span class="anchor" id="section3"></span>
        <h2 id="gene_distributions">Gene distributions</h2>
        <p>
            <div class="img">
            <img src=\"rankedgenes.png\" alt="Ranked genes" id="ranked_genes">
            </div>
        </p>
        <p>
            <div class="img">
            <img src=\"cumulative.png\" alt="Cumulative genes" id="cumulative">
            </div>
        </p>
        <p>
            <div class="img">
            <img src=\"density.png\" alt="Genes density" id="genes_density">
            </div>
        </p>

        <span class="anchor" id="section4"></span>
        <h2 id="metagenic_classification">Metagenic classification</h2>
        <p>
            <div class="img">
            <img src=\"annotation_coding.png\" alt="Metagenic classification coding" id="annotation_coding">
            </div>
        </p>
        <p>
            <div class="img">
            <img src=\"annotation_noncoding.png\" alt="Noncoding classification" id="annotation_noncoding">
            </div>
        </p>

        <span class="anchor" id="section5"></span>
        <h2 id="tot_phase">Total phase distribution</h2>
        <p>
            <div class="img">
            <img src=\"tot_phase.png" alt="total phase plot" id="tot_phase_img">
            </div>
        </p>

        <span class="anchor" id="section6"></span>
        <h2 id="phase_rpf_distr">RPF phase distribution</h2>
        <p>
            <div class="img">
            <img src=\"rpf_phase.png" alt="rpf phase plot" id="rpf_phase_img">
            </div>
        </p>

        <span class="anchor" id="section7"></span>
        <h2 id="phase_relpos_distr">Phase - relative position distribution</h2>
        <p>
            <div class="img">
            <img src=\"phase_relpos_distr.png" alt="phase relpos distr" id="phase_relpos_distr_img">
            </div>
        </p>

        <span class="anchor" id="section8"></span>
        <h2 id="triplet_identity">Triplet identity plots</h2>
        <p>
            <div class="img">
            <img src=\"triplet_id.png" alt="triplet identity plots" id="triplet_id_img">
            </div>
        </p>

        """+codon_usage_part+"""

        <br><br>
    </div>

    <div id="footer">
        <p id="footer_content">Generated with PROTEOFORMER - BioBix lab Ghent (Belgium) - Steven Verbruggen</p>
    </div>

</body>
</html>
    """

    #Generate html file
    html_file = open(outfile, 'w')
    html_file.write(html_string)
    html_file.close()

    return

## Codon usage plot
def codon_usage_plot(tmpfolder, codon_ref_file, outfolder, exp_name, treated, triplet_distr):

    #Output file
    output_file = outfolder + "/codon_usage.png"

    # Read reference
    reference = read_ref(codon_ref_file)

    # Read input codon
    name = exp_name
    codon_perc = read_codon_count(triplet_distr)

    # Sort triplets based on reference percentages
    sorted_triplets = []
    for i in sorted(reference, key=lambda k: reference[k], reverse=True):
        sorted_triplets.append(i)

    # Plot
    plot_codon_perc(output_file, sorted_triplets, reference, codon_perc, name)

    return

#Plot codon percentages
def plot_codon_perc(output_file, sorted_triplets, reference, codon_percs, name):

    #Parse data
    sorted_ref_values = []
    for i in sorted_triplets:
        sorted_ref_values.append(reference[i])
    xpos = range(1, len(sorted_ref_values)+1)

    sorted_codon_values = []
    for i in sorted_triplets:
        sorted_codon_values.append(codon_percs[i])

    #Labels
    labels=[]
    codontable = get_codontable()
    for i in sorted_triplets:
        labels.append(i+" ("+codontable[i]+")")

    #Plot
    fig = plt.figure(figsize=(36,32))
    ax = plt.axes()
    points = ax.plot(xpos, sorted_ref_values, marker='o', linewidth=15, color='#3BBE71', markersize=20, alpha=0.7, label='Reference')
    bars = ax.bar(xpos, sorted_codon_values, color='#228EDA', label=name)
    plt.xlim([0,len(sorted_ref_values)+1])
    [y1, y2] = ax.get_ylim()
    plt.ylim([0, y2])
    ax.set_ylabel('Percentage codon usage [in %]', fontsize=40)
    ax.set_xlabel('Triplet (amino acid)', fontsize=40)
    plt.yticks(fontsize=30)
    ax.set_xticks(xpos)
    ax.set_xticklabels(labels, rotation='vertical', fontsize=30)
    ax.set_title('Codon usage', fontsize=80)
    plt.legend(fontsize=40)

    plt.tight_layout()

    fig.savefig(output_file)

    return

#Read codon count
def read_codon_count(triplet_distr):

    #Init
    codon_perc = defaultdict()
    counts_per_triplet = defaultdict()
    total_sum=0

    #Total sum
    for triplet in triplet_distr.keys():
        for phase in triplet_distr[triplet].keys():
            total_sum = total_sum + int(triplet_distr[triplet][phase])

    #Parse
    for triplet in triplet_distr.keys():
        counts_per_triplet[triplet]=0
        for phase in triplet_distr[triplet].keys():
            counts_per_triplet[triplet] = counts_per_triplet[triplet] + int(triplet_distr[triplet][phase])
        codon_perc[triplet] = float(counts_per_triplet[triplet])/total_sum*100

    return codon_perc

#Read reference codon usage
def read_ref(input_ref):

    #Init
    ref = defaultdict()

    with open(input_ref, 'r') as FR:
        lines = FR.readlines()
        for line in lines:
            line.rstrip("\n")
            (triplet, perc) = re.split(',', line)
            ref[triplet] = float(perc)*100

    return ref

## Plot triplet identity data
def triplet_plots(data, outputfolder):
    outfile = outputfolder+"/triplet_id.png"

    sns.set_palette('terrain') #Palette
    fig = plt.figure(figsize=(36, 32))
    grid = GridSpec(8, 9) #Construct grid for subplots
    grid_i = -1 #Grid coordinates
    grid_j = 0
    for triplet in sorted(data.keys(), key=lambda e: (get_AA(e), e)):
        grid_i += 1
        if grid_i == 8:
            grid_i = 0
            grid_j += 1
        ax = plt.subplot(grid[grid_j,grid_i]) #Define subplot axes element in grid
        df = pd.DataFrame.from_dict(data[triplet], orient="index") #By index for right orientation
        df = df.sort_index(axis=0)
        labels_list = df[0].values.tolist()
        labels_list = map(format_thousands, labels_list)
        df.plot(kind="pie", subplots="True", autopct='%.1f%%', ax=ax, legend=None, labels=labels_list, pctdistance=0.7, fontsize=16, textprops={'fontsize':16})
        #Individual legends off, subplots="True" for evading selecting y column error, autopercentage
        ax.set(ylabel='') #Do not plot y column label
        if triplet=="ATG":
            title_color=hex2color('#00ff00')
        elif triplet=="TGA" or triplet=="TAG" or triplet=="TAA":
            title_color=hex2color("#ff0000")
        else:
            title_color="k"
        ax.set_title(triplet+": "+get_AA(triplet), {'fontsize': 38}, color=title_color) #Put triplet in title
    handles, labels = ax.get_legend_handles_labels() #Get legend information of last ax object
    leg = fig.legend(handles, ["Phase 0", "Phase 1", "Phase 2"], bbox_to_anchor=(1, 0.53), fontsize=38)#Define legend
    leg.get_frame().set_edgecolor('b')
    plt.tight_layout(rect=(0.013,0,1,1)) #Prevent overlapping elements
    fig.savefig(outfile)

    return


def phase_position_distr(tmpfolder, outfolder, treated):

    #Input data, read in to pandas data frame
    inputdata_adress = tmpfolder+"/mappingqc_"+treated+"/pos_table_all.csv"
    inputdata = pd.read_csv(inputdata_adress, sep=',', header=None, names=["phase", "rel_position"])

    #Split data based on phase
    data0 = inputdata[inputdata["phase"]==0]["rel_position"]
    data1 = inputdata[inputdata["phase"] == 1]["rel_position"]
    data2 = inputdata[inputdata["phase"] == 2]["rel_position"]

    #Define 20 bins
    freq0, bin_edges0 = np.histogram(data0, bins=20, range=(0,1))
    freq1, bin_edges1 = np.histogram(data1, bins=20, range=(0,1))
    freq2, bin_edges2 = np.histogram(data2, bins=20, range=(0,1))

    #Plot data
    fig, ax = plt.subplots(1, 1, figsize=(36,32))
    bar1 = ax.bar(bin_edges0[:-1]+0.00625, freq0, 0.0125, color='#228EDA', edgecolor='none')
    bar2 = ax.bar(bin_edges1[:-1]+0.0125+0.00625, freq1, 0.0125, color='#3BBE71', edgecolor='none')
    bar3 = ax.bar(bin_edges2[:-1]+0.025+0.00625, freq2, 0.0125, color='#B7E397', edgecolor='none')
    try:
        ax.set_facecolor("#f2f2f2")
    except:
        ax.set_axis_bgcolor("#f2f2f2")
    lgd = ax.legend((bar1[0], bar2[0], bar3[0]),('Phase 0','Phase 1', 'Phase 2'), loc='center left', bbox_to_anchor=(1, 0.5), fontsize=38)

    #Axis info
    plt.ylabel("Counts", fontsize=38)
    plt.xlabel("Relative position in sequence", fontsize=38)
    plt.xlim([0, 1])
    ax.tick_params(labelsize=34)

    #Save output
    plt.tight_layout()
    fig.savefig(outfolder+"/phase_relpos_distr.png", bbox_extra_artists=(lgd,), bbox_inches='tight') #Make room for legend

    #Close matplotlib environment
    plt.close()

    return


## Make plot of total phase distribution
def plot_total_phase(distr, outfile):

    #Define figure and axes
    sns.set_style(style="whitegrid")
    sns.set_palette("terrain")
    fig, ax = plt.subplots(1, 1, figsize=(36,32))

    #Parse data into arrays
    x = [0, 1, 2]
    y = [distr[k] for k in sorted(distr.keys())]

    #Set exponent base of y ticks
    majorFormatter = FixedOrderFormatter(6)
    ax.yaxis.set_major_formatter(majorFormatter)
    ax.yaxis.offsetText.set_fontsize(36)

    #Make plot
    sns.barplot(x, y, ax=ax, edgecolor='none')

    #Axis labels
    plt.xlabel('Phase', fontsize=38)
    plt.ylabel('Counts', fontsize=38)
    ax.tick_params(labelsize=34)

    #Remove box lines around plot
    sns.despine()

    #Face color
    try:
        ax.set_facecolor("#f2f2f2")
    except:
        ax.set_axis_bgcolor("#f2f2f2")

    #Finish plot
    plt.tight_layout()

    #Save output
    fig.savefig(outfile)

    return


class FixedOrderFormatter(ScalarFormatter):
    """Formats axis ticks using scientific notation with a constant order of
    magnitude"""
    def __init__(self, order_of_mag=0, useOffset=True, useMathText=False):
        self._order_of_mag = order_of_mag
        ScalarFormatter.__init__(self, useOffset=useOffset,
                                 useMathText=useMathText)
    def _set_orderOfMagnitude(self, range):
        """Over-riding this to avoid having orderOfMagnitude reset elsewhere"""
        self.orderOfMagnitude = self._order_of_mag

## Make plot of RPF against phase as a grouped 2D bar chart
def plot_rpf_phase_grouped2D(phase_distr, outfile):
    # Parse data into Pandas data frame
    df = pd.DataFrame.from_dict(phase_distr, orient="index")
    df = df[['0', '1', '2']]  # Reorder columns
    df = df.stack()  # Put all phases as separated observations
    df = df.reset_index()
    df.columns = ['RPF length', 'Phase', 'Count']

    # Set figure
    sns.set_style(style="whitegrid")
    sns.set_palette("terrain")
    fig, ax = plt.subplots(1, 1)

    # Background color
    try:
        ax.set_facecolor("#f2f2f2")
    except:
        ax.set_axis_bgcolor("#f2f2f2")

    # Plot
    sns.factorplot(x='RPF length', y='Count', hue='Phase', data=df, ax=ax, kind="bar")
    sns.despine()

    # Y label
    ax.set_ylabel("Count")

    # Legend position
    ax.legend(title="Phase", loc='center right', bbox_to_anchor=(1.12, 0.5), ncol=1)

    # Finish plot
    plt.tight_layout()

    # Save figure
    fig.savefig(outfile)

    return

## Make plot of RPF against phase with mayavi
def plot_rpf_phase_mayavi(phase_distr, outfile):
    from mayavi import mlab

    # Parameters
    lensoffset = 0.5
    spread_factor = 2

    # Init data structures
    s = []
    x = [0, 1, 2]

    # Parse
    for rpf in sorted(phase_distr.keys()):
        data_row = []
        for phase in x:
            data_row.append(phase_distr[rpf][str(phase)])
        s.append(data_row)
    s = np.array(s)
    y = sorted(phase_distr.keys(), reverse=True)

    # Z data scaling parameters
    mean_s = np.mean(s)
    order_s = math.floor(math.log10(mean_s))
    max_s_axis = int(math.ceil(np.max(s) / (10 ** order_s)))

    # Bottom mesh
    xpos, ypos = np.meshgrid(x, y)

    # Figure
    mlab.figure(size=(400000, 320000), bgcolor=(1, 1, 1))

    # Colors
    # Color map
    cs = []
    for i in range(0, 220, int(220 / len(y)) + 1):
        clr = cm.terrain(i)
        cs.append(clr[0:3])  # Save as RGB, not RGBA

    # Bars
    s_plot = s / (10 ** order_s)
    for i in range(0, len(y)):
        RPF = ypos[i]
        phase = xpos[i]
        data_z = s_plot[i]
        color = cs[i]
        mlab.barchart(RPF * spread_factor, phase * spread_factor, data_z, color=color, line_width=10)

    # Axes
    yy = np.arange((min(x) - 0.25) * spread_factor, (max(x) + lensoffset) * spread_factor, 0.01)
    yx = np.repeat((max(y) + lensoffset) * spread_factor, len(yy))
    yz = np.repeat(0, len(yy))

    xx = np.arange((min(y) - 0.25) * spread_factor, (max(y) + lensoffset) * spread_factor, 0.01)
    xy = np.repeat((max(x) + lensoffset) * spread_factor, len(xx))
    xz = np.repeat(0, len(xx))

    z_ax_max = math.ceil(np.max(s) / (10 ** order_s))
    zz = np.arange(0, z_ax_max, 0.01)
    zx = np.repeat((min(y) - 0.25) * spread_factor, len(zz))
    zy = np.repeat((max(x) + lensoffset) * spread_factor, len(zz))

    mlab.plot3d(yx, yy, yz, line_width=0.01, tube_radius=0.01, color=(0, 0, 0))
    mlab.plot3d(zx, zy, zz, line_width=0.01, tube_radius=0.01, color=(0, 0, 0))
    mlab.plot3d(xx, xy, xz, line_width=0.01, tube_radius=0.01, color=(0, 0, 0))

    # Axes labels
    mlab.text3d((np.mean(y) + 1) * spread_factor, (max(x) + lensoffset + 1.5) * spread_factor, 0, 'RPF',
                color=(0, 0, 0), scale=0.6, orient_to_camera=False, orientation=[0, 0, 180])
    mlab.text3d((max(y) + lensoffset + 1.25) * spread_factor, (np.mean(x) - 1) * spread_factor, 0, 'Phase',
                color=(0, 0, 0), scale=0.6, orient_to_camera=False, orientation=[0, 0, 90])
    mlab.text3d((min(y) - 0.5) * spread_factor, (max(x) + 1.75) * spread_factor, np.max(s_plot - 4) / 2,
                'Counts (10^{0:g})'.format(order_s), color=(0, 0, 0), scale=0.6, orient_to_camera=False,
                orientation=[180, 90, 45])

    # Axes ticks
    # Phase
    for i in x:
        mlab.text3d((max(y) + lensoffset + 0.5) * spread_factor, i * spread_factor, 0, str(i), color=(0, 0, 0),
                    scale=0.5)
        xx_tick = np.arange((min(y) - 0.25) * spread_factor, (max(y) + lensoffset) * spread_factor + 0.2, 0.01)
        xy_tick = np.repeat(i * spread_factor, len(xx_tick))
        xz_tick = np.repeat(0, len(xx_tick))
        mlab.plot3d(xx_tick, xy_tick, xz_tick, line_width=0.01, tube_radius=0.01, color=(0.7, 0.7, 0.7))
    # RPF
    for j in range(0, len(y)):
        mlab.text3d((j + min(y) + 0.35) * spread_factor, (max(x) + lensoffset + 0.5) * spread_factor, 0, str(y[j]),
                    color=(0, 0, 0), scale=0.5)
        yy_tick = np.arange((min(x) - 0.25) * spread_factor, (max(x) + lensoffset) * spread_factor + 0.2, 0.01)
        yx_tick = np.repeat((j + min(y)) * spread_factor, len(yy_tick))
        yz_tick = np.repeat(0, len(yy_tick))
        mlab.plot3d(yx_tick, yy_tick, yz_tick, line_width=0.01, tube_radius=0.01, color=(0.7, 0.7, 0.7))
    # Counts
    for k in range(0, max_s_axis + 1, 2):
        mlab.text3d((min(y) - 0.5) * spread_factor, (max(x) + 0.75) * spread_factor, k, str(k), color=(0, 0, 0),
                    scale=0.5)
        xx_back = np.arange((min(y) - 0.25) * spread_factor, (max(y) + lensoffset) * spread_factor, 0.01)
        xy_back = np.repeat((min(x) - 0.25) * spread_factor, len(xx_back))
        xz_back = np.repeat(k, len(xx_back))
        mlab.plot3d(xx_back, xy_back, xz_back, line_width=0.01, tube_radius=0.01, color=(0.7, 0.7, 0.7))
        yy_side = np.arange((min(x) - 0.25) * spread_factor, (max(x) + lensoffset) * spread_factor, 0.01)
        yx_side = np.repeat((min(y) - 0.25) * spread_factor, len(yy_side))
        yz_side = np.repeat(k, len(yy_side))
        mlab.plot3d(yx_side, yy_side, yz_side, line_width=0.01, tube_radius=0.01, color=(0.7, 0.7, 0.7))

    # Side surfaces
    y_wall, z_wall = np.mgrid[(min(y) - 0.25) * spread_factor:(max(y) + lensoffset) * spread_factor:2j, 0:z_ax_max:2j]
    x_wall = np.zeros_like(y_wall) - lensoffset

    x_wall2, z_wall2 = np.mgrid[(min(x) - 0.25) * spread_factor:(max(x) + lensoffset) * spread_factor:2j, 0:z_ax_max:2j]
    y_wall2 = np.zeros_like(x_wall2) + ((min(y) - 0.25) * spread_factor)

    y_wall3, x_wall3 = np.mgrid[(min(y) - 0.25) * spread_factor:(max(y) + lensoffset) * spread_factor:2j,
                       (min(x) - 0.25) * spread_factor:(max(x) + lensoffset) * spread_factor:2j, ]
    z_wall3 = np.zeros_like(y_wall3)

    mlab.mesh(y_wall, x_wall, z_wall, color=(1, 1, 1))
    mlab.mesh(y_wall2, x_wall2, z_wall2, color=(1, 1, 1))
    mlab.mesh(y_wall3, x_wall3, z_wall3, color=(1, 1, 1))

    # Camera
    mlab.view(azimuth=55, elevation=65, distance='auto')

    # Export figure
    mlab.savefig(outfile)

    return

## Make plot of RPF against phase
def plot_rpf_phase_pyplot3D(phase_distr, outfile):
    from mpl_toolkits.mplot3d import Axes3D

    #Initialize fig and axes
    fig = plt.figure(figsize=(20,13))
    ax = fig.gca(projection='3d')
    axis_rate = 3 / float(len(phase_distr.keys()))
    ax.pbaspect = [axis_rate, 1, 0.5]  # Aspect ratio based on proportions of x and y axis
    '''
    If you want no square axis:
    edit the get_proj function inside site-packages\mpl_toolkits\mplot3d\axes3d.py (and remove .pyc file!):

    try:
        self.localPbAspect=self.pbaspect
        zoom_out = (self.localPbAspect[0]+self.localPbAspect[1]+self.localPbAspect[2])
    except AttributeError:
        self.localPbAspect=[1,1,1]
        zoom_out = 0
    xmin, xmax = self.get_xlim3d() /  self.localPbAspect[0]
    ymin, ymax = self.get_ylim3d() /  self.localPbAspect[1]
    zmin, zmax = self.get_zlim3d() /  self.localPbAspect[2]

    # transform to uniform world coordinates 0-1.0,0-1.0,0-1.0
    worldM = proj3d.world_transformation(xmin, xmax,
                                             ymin, ymax,
                                             zmin, zmax)

    # look into the middle of the new coordinates
    R = np.array([0.5*self.localPbAspect[0], 0.5*self.localPbAspect[1], 0.5*self.localPbAspect[2]])
    xp = R[0] + np.cos(razim) * np.cos(relev) * (self.dist+zoom_out)
    yp = R[1] + np.sin(razim) * np.cos(relev) * (self.dist+zoom_out)
    zp = R[2] + np.sin(relev) * (self.dist+zoom_out)
    E = np.array((xp, yp, zp))


    then add one line to this script to set pbaspect:

    ax = fig.gca(projection='3d')
    ax.pbaspect = [2.0, 0.6, 0.25] #e.g.
    '''

    #Make arrays for x- and y-coordinates
    x = np.asarray([0,1,2])
    y = np.asarray(sorted(phase_distr.keys()))
    xmin = min(x)
    xmax = max(x)
    ymin=min(y)
    ymax=max(y)

    #Make a grid so that x- and y-coordinate for each point are given
    xpos, ypos = np.meshgrid(x[:] - 0.125, y[:] - 0.375)

    #Flatten: transform 2 dimensional array to 1 dimensional
    xpos = xpos.flatten('F')
    ypos = ypos.flatten('F')
    zpos = np.zeros_like(xpos)  # Make for z-coordinate an array of zeroes of the same length

    #Length of the bars
    dx = 0.5 * np.ones_like(xpos)
    dy = dx.copy()

    #Sort phase distr dict based on keys and save as z lengths
    keys_sorted_1 = sorted(phase_distr.keys())
    keys_sorted_2 = sorted(phase_distr[keys_sorted_1[0]].keys())
    dz = []


    for key2 in keys_sorted_2:
        for key1 in keys_sorted_1:
            dz.append(phase_distr[key1][key2])

    #Color map
    cs=[]
    for i in range(0,300,int(300/len(y))+1):
        clr = cm.terrain(i)
        cs.append(clr)

    #Plot
    bars = ax.bar3d(xpos, ypos, zpos, dx, dy, dz, zsort='max', alpha=0.99, color=cs*3)
    bars.set_edgecolor('k')

    #Axis ticks
    ticksxpos = np.arange(xmin, xmax+1, 1)
    ticksxlabs =range(xmin,xmax+1,1)
    plt.xticks(ticksxpos,ticksxlabs)
    ax.tick_params(axis='x', which='major', labelsize=20, pad=15)

    ticksypos = np.arange(ymin, ymax+1,1)
    ticksylabs=range(ymin,ymax+1,1)
    plt.yticks(ticksypos,ticksylabs)
    ax.tick_params(axis='y', which='major', labelsize=20, pad=15)

    ax.tick_params(axis='z', which='major', labelsize=20, pad=40)

    #Axis labels
    plt.xlabel('Phase', labelpad=60, fontsize=30)
    plt.ylabel('RPF length', labelpad=100, fontsize=30)
    ax.set_zlabel('counts', labelpad=90, fontsize=30)

    #Camera point angle
    ax.view_init(azim=-30)

    #Place axes in the middle of the fig environment
    ax.set_position([-0.35, -0.15, 1.5, 1.5])

    #Save as output
    fig.savefig(outfile)

    plt.close()

    return

## Get plot data out of results DB
def get_plot_data(db, treated):

    #Parse
    fastq=""
    if treated == "untreated":
        fastq = "fastq1"
    elif treated == "treated":
        fastq = "fastq2"
    else:
        print "ERROR: treated option cannot be parsed!"
        sys.exit()

    #Init
    phase_distr = defaultdict(lambda: defaultdict())
    total = defaultdict()
    total['0'] = 0
    total['1'] = 0
    total['2'] = 0
    triplet_data = defaultdict(lambda: defaultdict())

    try:
        conn = sqlite3.connect(db)
    except:
        print "ERROR: could not connect to "+db
        sys.exit()

    with conn:
        cur = conn.cursor()

        #Query rpf phase
        query = "SELECT * FROM rpf_phase_"+treated+";"
        cur.execute(query)
        output = cur.fetchall()

        for i in range(0, len(output)):
            #Parse
            rpf = output[i][0]
            phase0count = output[i][1]
            phase1count = output[i][2]
            phase2count = output[i][3]
            #RPF splitted phase distribution
            phase_distr[rpf]['0'] = phase0count
            phase_distr[rpf]['1'] = phase1count
            phase_distr[rpf]['2'] = phase2count
            #Total phase distribution
            total['0'] = total['0'] + phase0count
            total['1'] = total['1'] + phase1count
            total['2'] = total['2'] + phase2count

        #Get triplet data
        query = "SELECT triplet, phase, count FROM phase_triplet_"+treated+";"
        cur.execute(query)
        output_triplet = cur.fetchall()

        for i in range(0, len(output_triplet)):
            # Parse
            triplet = output_triplet[i][0]
            phase = output_triplet[i][1]
            count = output_triplet[i][2]
            triplet_data[triplet][phase] = count

        #Get arguments table variables
        query = "SELECT value FROM arguments WHERE variable=\"run_name\";"
        cur.execute(query)
        run_name = cur.fetchone()[0]
        query = "SELECT value FROM arguments WHERE variable=\"mapper\";"
        cur.execute(query)
        mapper = cur.fetchone()[0]
        query = "SELECT value FROM arguments WHERE variable=\"unique\";"
        cur.execute(query)
        unique = cur.fetchone()[0]
        query = "SELECT value FROM arguments WHERE variable=\"ensembl_version\";"
        cur.execute(query)
        ensembl_version = cur.fetchone()[0]
        query = "SELECT value FROM arguments WHERE variable=\"species\";"
        cur.execute(query)
        species = cur.fetchone()[0]

        #Get total mapped genomic sequences
        sample_name = run_name+"_"+mapper+"_"+unique+"_"+ensembl_version+"."+fastq
        gene_distr_name_prefix = run_name+"_ens"+str(ensembl_version)+"_"+mapper
        query = "SELECT mapped_T FROM statistics WHERE sample=\""+sample_name+"\" AND type=\"genomic\";"
        cur.execute(query)
        tot_maps = cur.fetchone()[0]

    return phase_distr, total, triplet_data, run_name, tot_maps, gene_distr_name_prefix, ensembl_version, species

def format_thousands(x):
    return '{:,}'.format(int(x)).replace(",", " ")


def get_AA(item):

    codontable = get_codontable()
    AA = codontable[item]

    return AA


def get_species(result_db):

    try:
        conn = sqlite3.connect(result_db)
    except:
        print "Could not connect to "+result_db
        sys.exit()

    with conn:
        cur = conn.cursor()

        query = "SELECT value FROM arguments WHERE variable='species';"
        cur.execute(query)
        species = cur.fetchone()[0]

    return species

def get_codontable():

    codontable = {
        'ATA': 'Ile', 'ATC': 'Ile', 'ATT': 'Ile', 'ATG': 'Met',
        'ACA': 'Thr', 'ACC': 'Thr', 'ACG': 'Thr', 'ACT': 'Thr',
        'AAC': 'Asn', 'AAT': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
        'AGC': 'Ser', 'AGT': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
        'CTA': 'Leu', 'CTC': 'Leu', 'CTG': 'Leu', 'CTT': 'Leu',
        'CCA': 'Pro', 'CCC': 'Pro', 'CCG': 'Pro', 'CCT': 'Pro',
        'CAC': 'His', 'CAT': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
        'CGA': 'Arg', 'CGC': 'Arg', 'CGG': 'Arg', 'CGT': 'Arg',
        'GTA': 'Val', 'GTC': 'Val', 'GTG': 'Val', 'GTT': 'Val',
        'GCA': 'Ala', 'GCC': 'Ala', 'GCG': 'Ala', 'GCT': 'Ala',
        'GAC': 'Asp', 'GAT': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
        'GGA': 'Gly', 'GGC': 'Gly', 'GGG': 'Gly', 'GGT': 'Gly',
        'TCA': 'Ser', 'TCC': 'Ser', 'TCG': 'Ser', 'TCT': 'Ser',
        'TTC': 'Phe', 'TTT': 'Phe', 'TTA': 'Leu', 'TTG': 'Leu',
        'TAC': 'Tyr', 'TAT': 'Tyr', 'TAA': 'STOP', 'TAG': 'STOP',
        'TGC': 'Cys', 'TGT': 'Cys', 'TGA': 'STOP', 'TGG': 'Trp',
    }

    return codontable


## Data Dumper for recursively printing nested dictionaries and defaultDicts ##
def print_dict(dictionary, indent='', braces=0):
    """
    :param dictionary: dict or defaultdict to be printed
    :param indent: indentation
    :param braces:
    :return: void
    """

    for key, value in dictionary.iteritems():
        if isinstance(value, dict):
            print '%s%s%s%s' % (indent, braces * '[', key, braces * ']')
            print_dict(value, indent + '  ', braces)
        else:
            print indent + '%s = %s' % (key, value)

#######Set Main##################
if __name__ == "__main__":
    try:
        main()
    except Exception, e:
        traceback.print_exc()
#################################