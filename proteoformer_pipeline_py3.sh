#!/bin/bash


echo ' ________  ________  ________  _________  _______   ________  ________ ________  ________  _____ ______   _______   ________     '
echo '|\   __  \|\   __  \|\   __  \|\___   ___\\  ___ \ |\   __  \|\  _____\\   __  \|\   __  \|\   _ \  _   \|\  ___ \ |\   __  \    '
echo '\ \  \|\  \ \  \|\  \ \  \|\  \|___ \  \_\ \   __/|\ \  \|\  \ \  \__/\ \  \|\  \ \  \|\  \ \  \\\__\ \  \ \   __/|\ \  \|\  \   '
echo ' \ \   ____\ \   _  _\ \  \\\  \   \ \  \ \ \  \_|/_\ \  \\\  \ \   __\\ \  \\\  \ \   _  _\ \  \\|__| \  \ \  \_|/_\ \   _  _\  '
echo '  \ \  \___|\ \  \\  \\ \  \\\  \   \ \  \ \ \  \_|\ \ \  \\\  \ \  \_| \ \  \\\  \ \  \\  \\ \  \    \ \  \ \  \_|\ \ \  \\  \| '
echo '   \ \__\    \ \__\\ _\\ \_______\   \ \__\ \ \_______\ \_______\ \__\   \ \_______\ \__\\ _\\ \__\    \ \__\ \_______\ \__\\ _\ '
echo '    \|__|     \|__|\|__|\|_______|    \|__|  \|_______|\|_______|\|__|    \|_______|\|__|\|__|\|__|     \|__|\|_______|\|__|\|__|'
echo -e "\n\n"                                                                                                                                 
                                                                                                                                 
                                                                                                                                 


## GLOBAL VARIABLES

declare -A datasets=(["OHMX202X0XXX_00X"]="/path/to/fastq/gz/file1.fastq.gz"
                     ["OHMX202X0XXX_00X"]="/path/to/fastq/gz/file2.fastq.gz")

##Input arguments##
echo "INPUT ARUGMENTS:"
ENSEMBL_ANNOT="100"
SPECIES="human"
SPECIES_SHORT="hsa"
SUITE="standard"  											# Mostly standard or plastid
PRICE="Y"
CORES="20"
UNIQUEMAPPING="Y"
CLIPPER="trimmomatic" 									# Mostly trimmomatic or fastx
ADAPTORSEQ="TGGAATTCTCGGGTGCCAAGG"
PHIX="Y"												# Y or N
RRNA="Y"												# Y or N
SNORNA="Y"												# Y or N
TRNA="Y"												# Y or N
SAVE_UNMAPPED="N"                                       # Y or N
TRCOORD="Y"                                             # Y or N
TESTRUN="N"                                             # Y or N
IGENOMESROOT="/data/igenomes/"
ENSEMBLDICT="/share/user/Ensembl/"
ENSEMBLDB="ENS_${SPECIES_SHORT}_${ENSEMBL_ANNOT}.db"
COMPLOGO='ohmx'											# ohmx or biobix
SCRIPTDIR="/home/user/scripts/"
BASEDIR="/data/user/project_folder/folder/"
echo "basedir = $BASEDIR"
echo -e "scriptdir = $SCRIPTDIR\n"

echo "Ensembl annotation = $ENSEMBL_ANNOT"
echo "Species = $SPECIES"
echo "Offset calling method = $SUITE"
echo "PRICE-adapted mapping files = $PRICE"
echo "Cores = $CORES"
echo "Unique mapping = $UNIQUEMAPPING"
echo "Adaptor sequence = $ADAPTORSEQ"
echo "Clipper = $CLIPPER"
echo "Phix filtering = $PHIX"
echo "rRNA filtering = $RRNA"
echo "tRNA filtering = $TRNA"
echo "sn(o)RNA filtering = $SNORNA"
echo "Test run = $TESTRUN"
echo "iGenomes root folder = $IGENOMESROOT"
echo "Ensembl DB location = ${ENSEMBLDICT}${ENSEMBLDB}" 
echo -e "\n\n"

#Activate conda for in-shell usage
eval "$(conda shell.bash hook)"

##Activate PROTEOFORMER Conda Environment##
conda activate proteoformer_general

##Reference information##
#echo -e "Download reference info\n\n"
#python $SCRIPTDIR/proteoformer/Additional_tools/ENS_db.py -v $ENSEMBL_ANNOT -s $SPECIES
#chmod 755 $ENSEMBLDB
#mv $ENSEMBLDB $ENSEMBLDICT
#python $SCRIPTDIR/proteoformer/Additional_tools/get_igenomes.py -v $ENSEMBL_ANNOT -s $SPECIES -d $IGENOMESROOT -c 15


##Internal dict structure##
mkdir $BASEDIR/fastqc_raw
mkdir $BASEDIR/fastqc_mapped
mkdir $BASEDIR/statistics
mkdir $BASEDIR/mqc
mkdir $BASEDIR/mqc_suppl

# Cat files for statistics
TOTAL_STATS_FILE=${BASEDIR}/statistics/total_stats.csv
TEMP_TOTAL_STATS=${BASEDIR}/statistics/tmp_total_stats.csv
rm -rf ${TOTAL_STATS_FILE}
touch ${TOTAL_STATS_FILE}

##Start loop
for i in "${!datasets[@]}"
do
  :

ID=$i
FILE=${datasets[$i]}
echo -e "START NEW LOOP:\n\t${ID}\n\t${FILE}\n"

##gunzip if necessary
if [[ "$FILE" =~ .*\.gz$ ]]
then
    pigz -p $CORES -d $FILE
    #gunzip $FILE
    UNZIPFILE="${FILE%.*}"
    echo "File $FILE unzipped"
else
	UNZIPFILE=$FILE
fi

echo -e "Unzipped file: $UNZIPFILE \n"

mkdir $ID
cd $ID

echo -e "1) FastQC raw file $ID \n"
fastqc $UNZIPFILE -o $BASEDIR/fastqc_raw -t $CORES
echo -e "FastQC raw file $ID done \n"



echo -e "2) Mapping $ID \n"
##Mapping specific arguments
READTYPE="ribo_untr"
OFFSET_FILE_UNTR="path_to_file"
OFFSET_FILE_TR="path_to_file"
OFFSET_IMG_UNTR=$BASEDIR/$ID/plastid/${ID}_untreated_p_offsets.png
OFFSET_IMG_TR=$BASEDIR/$ID/plastid/${ID}_treated_p_offsets.png

## Start with basic mapping
perl $SCRIPTDIR/proteoformer/1_mapping_py3/mapping.pl --inputfile1 $UNZIPFILE --readtype $READTYPE --name $ID --species $SPECIES --ensembl $ENSEMBL_ANNOT --cores $CORES --unique $UNIQUEMAPPING --igenomes_root $IGENOMESROOT --clipper $CLIPPER --adaptor $ADAPTORSEQ --phix $PHIX --rRNA $RRNA --snRNA $SNORNA --tRNA $TRNA --save_unmapped $SAVE_UNMAPPED --rpf_split N --tr_coord $TRCOORD --price_files $PRICE --suite $SUITE > $BASEDIR/$ID/mapping_$ID.txt 2>&1

##Then, go over the different suite options
if [[ $SUITE = "standard" ]] || [[ $SUITE = "cst_5prime" ]] || [[ $SUITE = "cst_3prime" ]]
then
    echo -e "\n\n\n\n\n\n\n\t\tS U I T E:     GO THROUGH WITH MAPPING PARSING\n\n\n"
    perl $SCRIPTDIR/proteoformer/1_mapping_py3/mapping_parsing.pl --out_sqlite $BASEDIR/$ID/SQLite/results.db --offset $SUITE > $BASEDIR/$ID/mapping_parsing_$ID.txt 2>&1
elif [[ $SUITE = "from_file" ]]
then
    if [[ $READTYPE = "ribo_untr" ]]
    then
        echo -e "\n\n\n\n\n\n\n\t\tS U I T E:     GO THROUGH WITH MAPPING PARSING\n\n\n"
        perl $SCRIPTDIR/proteoformer/1_mapping_py3/mapping_parsing.pl --out_sqlite $BASEDIR/$ID/SQLite/results.db --offset $SUITE --offset_file_untr $OFFSET_FILE_UNTR > $BASEDIR/$ID/mapping_parsing_$ID.txt 2>&1
    elif [[ $READTYPE = "ribo" ]]
    then
        echo -e "\n\n\n\n\n\n\n\t\tS U I T E:     GO THROUGH WITH MAPPING PARSING (treated+untreated)\n\n\n"
        perl $SCRIPTDIR/proteoformer/1_mapping_py3/mapping_parsing.pl --out_sqlite $BASEDIR/$ID/SQLite/results.db --offset $SUITE --offset_file_untr $OFFSET_FILE_UNTR --offset_file_tr $OFFSET_FILE_TR > $BASEDIR/$ID/mapping_parsing_$ID.txt 2>&1
    else
        echo -e "from_file suite can only occur with the ribo and ribo_untr readtypes!"
    fi
elif [[ $SUITE = "plastid" ]]
then
    conda deactivate
    conda activate proteoformer_plastid
    echo -e "\n\n\n\n\n\n\n\t\tS U I T E:     GO THROUGH WITH PLASTID (UNTREATED)\n\n\n"
    perl $SCRIPTDIR/proteoformer/1_mapping_py3/mapping_plastid.pl --out_sqlite $BASEDIR/$ID/SQLite/results.db --treated untreated  --offset_img $OFFSET_IMG_UNTR > $BASEDIR/$ID/mapping_plastid_$ID.txt 2>&1
    if [[ $READTYPE = "ribo" ]]
    then
        echo -e "\n\n\n\n\n\n\n\t\tS U I T E:     GO THROUGH WITH PLASTID (TREATED)\n\n\n"
        perl $SCRIPTDIR/proteoformer/1_mapping_py3/mapping_plastid.pl --out_sqlite $BASEDIR/$ID/SQLite/results.db --treated treated --offset_img $OFFSET_IMG_TR > $BASEDIR/$ID/mapping_plastid_tr_$ID.txt 2>&1
    fi
    conda deactivate
    conda activate proteoformer_general
    echo -e "\n\n\n\n\n\n\n\t\tS U I T E:     GO THROUGH WITH MAPPING PARSING\n\n\n"
    perl $SCRIPTDIR/proteoformer/1_mapping_py3/mapping_parsing.pl --out_sqlite $BASEDIR/$ID/SQLite/results.db --offset $SUITE > $BASEDIR/$ID/mapping_parsing_$ID.txt 2>&1
else
    echo "Suite option is not in the list of possible options!"
fi

#Copy stats for multiQC
cp STAR/fastq1/Log.final.out $BASEDIR/fastqc_mapped/${ID}_fastqc.Log.final.out
echo -e "Mapping done for $ID \n"

ln -s  $BASEDIR/$ID/STAR/fastq1/untreat.bam $BASEDIR/$ID/$ID.bam



echo -e "3) FastQC mapped file $ID \n"
fastqc $BASEDIR/$ID/$ID.bam -o $BASEDIR/fastqc_mapped -t $CORES
echo -e "FastQC on the mapped file done for $ID\n"

sqlite3 SQLite/results.db <<END_SQL
.timeout 2000
.header on
.mode csv
.output ${BASEDIR}/statistics/${ID}_stats.csv
select * from statistics order by total desc;
END_SQL

echo -e "Statistics written for $ID \n"

mkdir $BASEDIR/mqc/$ID
mkdir $BASEDIR/mqc_suppl/$ID
echo -e "4) mQC $ID \n"
perl $SCRIPTDIR/proteoformer/2_mappingQC_py3/mappingQC.pl --samfile $BASEDIR/$ID/STAR/fastq1/untreat.sam --treated untreated --testrun $TESTRUN --cores $CORES --result_db $BASEDIR/$ID/SQLite/results.db --unique $UNIQUEMAPPING --ens_db $ENSEMBLDICT/$ENSEMBLDB --offset $SUITE  --offset_img $OFFSET_IMG_UNTR --tool_dir $SCRIPTDIR/proteoformer/2_mappingQC_py3/mqc_tools/ --plotrpftool pyplot3D --output_folder $BASEDIR/mqc/$ID --suppl_out_folder $BASEDIR/mqc_suppl/$ID --html $BASEDIR/mqc/mqc_$ID.html --zip $BASEDIR/mqc/mqc_$ID.zip --tmp $BASEDIR/$ID/tmp/ --comp_logo $COMPLOGO > $BASEDIR/$ID/mQC_$ID.txt 2>&1
rm -rf $BASEDIR/$ID/tmp/mappingqc_untreated
echo -e "mQC performed for $ID \n"
echo -e "\n"

rm -rf  fastq/fastq1_nophix.fq fastq/nophix/  fastq/fastq1_clip* fastq/Unmapped.out.mate1  fastq/Log.* fastq/SJ.out.tab  fastq/Aligned.out.sam fastq/fastq1_norrna.fq fastq/fastq1_norrna_nosnrna.fq fastq/fastq1_norrna_nosnrna_notrna.fq

# echo -e "5) Transcript calling $ID \n"
# perl $SCRIPTDIR/proteoformer/3_tr_calling/Rule-based/ribo_translation.pl --in_sqlite $BASEDIR/$ID/SQLite/results.db --out_sqlite $BASEDIR/$ID/SQLite/results.db --ens_db $ENSEMBLDICT/$ENSEMBLDB > $BASEDIR/$ID/tr_translation_$ID.txt 2>&1
# echo -e "Transcript calling performed on $ID \n"


# ##To run price, another ENV needs to be loaded.
# conda deactivate
# conda activate price

# echo -e "6) PRICE $ID \n"
# python $SCRIPTDIR/proteoformer/4_ORF_calling/using_PRICE/PRICE.py -r $BASEDIR/$ID/SQLite/results.db > $BASEDIR/$ID/PRICE_orf_$ID.txt 2>&1
# echo -e "PRICE ORF calling performed on $ID \n"

# conda deactivate
# conda activate proteoformer


# ##To run spectre, another ENV needs to be loaded.
# conda deactivate
# conda activate spectre

# echo -e "7) SPECTRE $ID \n"
# python $SCRIPTDIR/proteoformer/4_ORF_calling/using_SPECtre/SPECtre.py -r $BASEDIR/$ID/SQLite/results.db -o 28:12,29:12,30:12 -c 60 -x 3 > $BASEDIR/$ID/SPECtre_orf_$ID.txt 2>&1
# echo -e "SPECTRE ORF calling performed on $ID \n"

# conda deactivate
# conda activate proteoformer

##Go back to the BaseDir
cd ..

pigz -p $CORES $UNZIPFILE
#gzip $UNZIPFILE
echo "File $UNZIPFILE zipped"
echo -e "\n\n"

done
#End loop

#MultiQC is in base environment
conda deactivate
conda activate proteoformer_multiqc

cd $BASEDIR/fastqc_raw
multiqc -i "Raw Read Data" .
cd $BASEDIR

cd $BASEDIR/fastqc_mapped
multiqc -i "Mapped Read Data" .
cd $BASEDIR

echo -e "MultiQC on all raw/mapped data performed\n"

conda deactivate

##Cat all stats files
#Do this in a sorted ID manner
readarray -t sorted_datasets < <(for a in "${!datasets[@]}"; do echo "$a"; done | sort)
for i in "${sorted_datasets[@]}"
do
    ID=$i
    #Delete the header lines
    tail -n +2 "${BASEDIR}/statistics/${ID}_stats.csv" > ${BASEDIR}/statistics/tmp_stats.csv
    cat ${TOTAL_STATS_FILE} ${BASEDIR}/statistics/tmp_stats.csv >> ${TEMP_TOTAL_STATS}
    rm -rf ${BASEDIR}/statistics/tmp_stats.csv
    mv ${TEMP_TOTAL_STATS} ${TOTAL_STATS_FILE}
done

##Deactivate conda environment
#conda deactivate

