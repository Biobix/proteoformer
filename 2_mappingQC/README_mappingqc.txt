##################################################################
#								 #
# MappingQC (PROTEOFORMER version): some additional information  #
#								 #
##################################################################


EXAMPLE

Main script is mappingQC.pl and you can run it as in following example:
perl ./mappingQC.pl --samfile yoursamfile.sam --treated untreated --cores 20 --result_db yourproteoformerresults.db --ens_db ENS_mmu_82.db --unique N --offset plastid --offset_img youruntreatedoffsetsimg.png --tool_dir mqc_tools

A Galaxy version of MappingQC and the whole PROTEOFORMER version is also available.

There is also a stand-alone version of mappingQC available at https://github.com/Biobix/mQC


INPUT PARAMETERS

* work_dir: working directory to run the scripts in (default: current working directory)
* samfile: path to the SAM file that comes out of the mapping script of PROTEOFORMER (mandatory)
* treated: whether you analyse the alignments of the treated (LTM/HARR/PMY) or the untreated (CHX/freezing/no treatment) sample
		Possible options: treated, untreated (default untreated)
* cores: the amount of cores to run the script on (integer, default: 5)
* result_db: path to the SQLite database where PROTEOFORMER mapping has stored all mapping results (mandatory)
* tmp: temporary folder for storing temporary files of mappingQC (default: work_dir/tmp)
* unique: whether to use only the unique alignments.
		Possible option: Y, N (default Y)
	If you mapped uniquely in the previous step of PROTEOFORMER, only option Y is possible. If you mapped non-uniquely, Y or N available.
* ens_db: path to the Ensembl SQLite database with annotation info. Use the ENS_db.py script to build this (this script is located in the Information folder). (mandatory)
* offset: the offset determination method.
		Possible options:		
			- plastid: use the Plastid offsets (Dunn et al. 2016) which were calculated during the mapping step of PROTEOFORMER
				The image from the offset calculation should be given in the —offset_img argument
			- standard: use the standard offsets from the paper of Ingolia et al. (2012) (default)
			- from_file: use offsets from an input file
				The offsets input file should be given in the —offset_file argument
* output_folder: the folder to store the output files (default: work_dir/mappingQC_output)
* tool_dir: folder with necessary additional mappingQC tools. More information below in the ‘dependencies’ section. (default: work_dir/mqc_tools)
* plotrpftool: the module that will be used for plotting the RPF-phase figure
		Possible options:
                	- grouped2D: use Seaborn to plot a grouped 2D bar chart (default)
                       	- pyplot3D: use mplot3d to plot a 3D bar chart
				This tool can suffer sometimes from Escher effects, as it tries to plot a 3D plot with the 2D software of pyplot and matplotlib.
                     	- mayavi: use the mayavi package to plot a 3D bar chart
				This tool only works on local systems with graphical cards.
* html: custom name for the output HTML file (default: work_dir/mappingqc_out.html)
* zip: custom name for output ZIP file (default: work_dir/mappingQC_(un)treated.zip)



OUTPUT

MappingQC makes a folder with different plots. The main output file in this folder however is an HTML file that combines all these figures in a nice overview. To transfer your results to another system or location, mappingQC compresses the results folder also to a zip file, so that all results can be unpacked together on the target location.
Following figures are included in the folder and in the overview HTML file:
* A table with practical information about the mappingQC analysis (analysis time, input parameters, size of the SAM file,..)
* Plastid offset information (if applicable)
* Gene distributions: counts summed over genes and shown in 3 different plots.
* Metagenic plots: counts summed over different annotations. For the non-coding supergroup, a detail plot is also given.
* Total phase distribution: shows you how the counts are divided over the three reading frames in the canonical translation products of canonical protein-coding transcripts.
* RPF-phase distribution: shows you how the counts are divided over the three reading frames, but also separated over the different RPF lengths. This is also based on the canonical translation products of canonical protein-coding transcripts.
* Phase - relative position distribution: shows you how the counts are divided over the three reading frames on a metagenic scale over all canonical reading frames. Relative position 0 means metagenically at the beginning of the sequences; relative position 1 means at the metagenic end of the sequences.
* Triplet identity plots: shows you the reading frame distribution separated for all possible codons. This can pick up if there are any reading frame distribution biases for specific codons. The resulting amino acid or Start/Stop signal of each codon is given as well.



DEPENDENCIES

As you can see in the command line, mappingQC relies on a tool directory with some additional tools. These include:
* metagenic_piecharts.R				An R tool to plot the metagenic piecharts in R
* quality_plots.R				An R tool to plot the gene distribution quality plots in R
* mappingQC.py					A python (Python2) script to plot all the other plots and assemble all the output in an HTML overview file.

MappingQC, like the rest of the PROTEOFORMER pipeline, relies also on SQLite and the sqlite3 command line tool for storing results in database structures.

MappingQC relies on following Perl modules which have to be installed on your system:
* DBI
* Getopt::Long
* Parallel::ForkManager
* CWD
* Data::Dumper (for debugging purposes)

Furthermore, mappingQC relies on following Python2 modules which have to be installed on your system:
* getopt
* defaultdict (collections)
* sqlite3
* pandas
* numpy
* matplotlib (including pyplot, colors, cm, gridspec, ticker and mplot3d)
* seaborn

!! For the 3D plot (counts as a function of phase and RPF length), you have to make an adaptation in the Python2 libraries of mplot3d. The default axes3d.py script (that will be installed if you download and install mplot3d, the one that your python2 actually uses!) needs to be replaced by the axes3d.py script you can find at https://github.com/Biobix/proteoformer/tree/master/MappingQC/mqc_tools/site-packages. You also need to delete the axes3d.pyc script!
mplot3d is not able to plot 3D barcharts in non-cubic environments and this adapted script will solve this issue.


MORE INFORMATION

For more information about mappingQC: contact Steven.Verbruggen@UGent.be or Gerben.Menschaert@UGent.be