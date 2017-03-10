import traceback
import getopt
from collections import defaultdict
import sqlite3
import os
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
    -o | --outfolder                        The output folder for plots
                                                (default mappingQC_(un)treated)
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
        myopts, args = getopt.getopt(sys.argv[1:], "w:r:s:t:o:h:z:p:i:e:", ["work_dir=", "result_db=", "input_samfile=", "treated=", "outfile=", "outhtml=", "outzip=", "plastid_option=", "plastid_img=" ,"ensembl_db="])
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
    if outfolder == '':
        outfolder = workdir+"/mappingQC_"+treated
    #Trim off last backslash
    m = re.search("^.*/(.*?)/$", outfolder)
    if m:
        outfolder = m.group(1)
    if outhtml == '':
        outhtml = workdir+"/mappingQC_out_"+treated+".html"
    if outzip == '':
        outzip = workdir+"/mappingQC_"+treated+".zip"
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

    ########
    # MAIN #
    ########

    #Make plot  directory
    if not os.path.exists(outfolder):
        os.system("mkdir "+outfolder)

    #Download biobix image
    os.system("wget --quiet \"http://www.nxtgnt.ugent.be/images/img/logos/BIOBIX_logo.png\"")
    os.system("mv BIOBIX_logo.png "+outfolder)

    #Get plot data out of results DB
    phase_distr, total_phase_distr, triplet_distr, run_name, totmaps, prefix_gene_distr, ensembl_version, species = get_plot_data(result_db, treated)

    #Make total phase distribution plot
    outfile = outfolder+"/tot_phase.png"
    plot_total_phase(total_phase_distr, outfile)

    #Make RPF-phase distribution plot
    outfile = outfolder+"/rpf_phase.png"
    plot_rpf_phase(phase_distr, outfile)

    #Make phase position distribution
    phase_position_distr(tmpfolder, outfolder, treated)

    #Make triplet identity plots
    triplet_plots(triplet_distr, outfolder)

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
    write_out_html(outhtml, outfolder, samfile, run_name, totmaps, plastid_option, offsets_file, offset_img, prefix_gene_distr, ensembl_version, species, ens_db, treated)
    #Make copy of html so that is certainly in the output folder as well
    os.system("cp "+outhtml+" "+outfolder)

    ##Archive and collect output
    #Make output archive
    output_arch = "MappingQC_archive/"
    os.system("mkdir " + output_arch)
    #Bring output html and output images folder to archive
    os.system("cp -r " + outhtml + " " + output_arch)
    os.system("cp -r " + outfolder + " " + output_arch)
    #zip output archive
    tmpZip = workdir+"tmp.zip"
    os.system("zip -r -q "+tmpZip+" "+output_arch)
    os.system("mv "+tmpZip+" "+outzip)

    return

############
### SUBS ###
############


## Write output html file
def write_out_html(outfile, output_folder, samfile, run_name, totmaps, plastid, offsets_file, offsets_img, prefix_gene_distr, ensembl_version, species, ens_db, treated):

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

    #Structure of html file
    html_string = """<!DOCTYPE html>
<html>
<head>
   <title>Mapping QC Report """+run_name+"""</title>
   <meta charset="utf-8"></meta>
   <meta name="description" content="Overview HTML of all mappingQC results"></meta>
   <link href="https://fonts.googleapis.com/css?family=Indie+Flower" rel="stylesheet">
   <style>
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
</head>

<body>
    <div id="header">
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
            <img src=\""""+prefix_gene_distr+"""rankedgenes.png\" alt="Ranked genes" id="ranked_genes">
            </div>
        </p>
        <p>
            <div class="img">
            <img src=\""""+prefix_gene_distr+"""cumulative.png\" alt="Cumulative genes" id="cumulative">
            </div>
        </p>
        <p>
            <div class="img">
            <img src=\""""+prefix_gene_distr+"""density.png\" alt="Genes density" id="genes_density">
            </div>
        </p>

        <span class="anchor" id="section4"></span>
        <h2 id="metagenic_classification">Metagenic classification</h2>
        <p>
            <div class="img">
            <img src=\""""+prefix_gene_distr+"""annotation_coding.png\" alt="Metagenic classification coding" id="annotation_coding">
            </div>
        </p>
        <p>
            <div class="img">
            <img src=\""""+prefix_gene_distr+"""annotation_noncoding.png\" alt="Noncoding classification" id="annotation_noncoding">
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


## Plot triplet identity data
def triplet_plots(data, outputfolder):
    outfile = outputfolder+"/triplet_id.png"

    plt.style.use('fivethirtyeight') #Stylesheet
    fig = plt.figure(figsize=(36, 32))
    grid = GridSpec(8, 9) #Construct grid for subplots
    grid_i = -1 #Grid coordinates
    grid_j = 0
    for triplet in sorted(data.keys(), key=get_AA):
        grid_i += 1
        if grid_i == 8:
            grid_i = 0
            grid_j += 1
        ax = plt.subplot(grid[grid_j,grid_i]) #Define subplot axes element in grid
        df = pd.DataFrame.from_dict(data[triplet], orient="index") #By index for right orientation
        df = df.sort_index(axis=0)
        labels_list = df[0].values.tolist()
        labels_list = map(format_thousands, labels_list)
        df.plot(kind="pie", subplots="True", autopct='%.1f%%', ax=ax, legend=None, labels=labels_list)
        #Individual legends off, subplots="True" for evading selecting y column error, autopercentage
        ax.set(ylabel='') #Do not plot y column label
        if triplet=="ATG":
            title_color=hex2color('#00ff00')
        elif triplet=="TGA" or triplet=="TAG" or triplet=="TAA":
            title_color=hex2color("#ff0000")
        else:
            title_color="k"
        ax.set_title(triplet+": "+get_AA(triplet), {'fontsize': 36}, color=title_color) #Put triplet in title
    handles, labels = ax.get_legend_handles_labels() #Get legend information of last ax object
    leg = fig.legend(handles, ["Phase 0", "Phase 1", "Phase 2"], bbox_to_anchor=(1, 0.53), fontsize=36)#Define legend
    leg.get_frame().set_edgecolor('b')
    plt.tight_layout() #Prevent overlapping elements
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
    bins = np.linspace(0,1,21)

    #Plot data
    fig, ax = plt.subplots(1, 1)
    ax.hist([data0,data1,data2], bins, label=["Phase 0", "Phase 1", "Phase 2"])
    ax.set_facecolor("#f2f2f2")
    lgd = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    #Axis info
    plt.ylabel("Counts")
    plt.xlabel("Relative position in sequence")
    plt.xlim([0, 1])

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
    fig, ax = plt.subplots(1, 1)

    #Parse data into arrays
    x = [0, 1, 2]
    y = [distr[k] for k in sorted(distr.keys())]

    #Set exponent base of y ticks
    majorFormatter = FixedOrderFormatter(6)
    ax.yaxis.set_major_formatter(majorFormatter)

    #Make plot
    sns.barplot(x, y, ax=ax)

    #Axis labels
    plt.xlabel('Phase')
    plt.ylabel('Counts')

    #Remove box lines around plot
    sns.despine()

    #Face color
    ax.set_facecolor("#f2f2f2")

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


## Make plot of RPF against phase
def plot_rpf_phase(phase_distr, outfile):

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
    xpos, ypos = np.meshgrid(x[:] - 0.125, y[:] -0.375)

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