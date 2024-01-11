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

declare -A datasets=(["OHMX20200619_005"]="/data/gerbenm/20200619_Novo_GST_Frame/Novo/20191030_Ribo1_MHU-L-041709_S177_L008_R1_001.fastq.gz"
                     ["OHMX20200619_006"]="/data/gerbenm/20200619_Novo_GST_Frame/Novo/20191030_Ribo2_MHU-L-061019_S178_L008_R1_001.fastq.gz"
                     ["OHMX20200619_007"]="/data/gerbenm/20200619_Novo_GST_Frame/Novo/20191030_Ribo3_MHU-L-091307_S179_L008_R1_001.fastq.gz"
                     ["OHMX20200619_001"]="/data/gerbenm/20200619_Novo_GST_Frame/Novo/a20190013_001_FKDL202574961-1a-AK2178-AK1031_HK7NYDRXX_L1.fq.gz"
                     ["OHMX20200619_004"]="/data/gerbenm/20200619_Novo_GST_Frame/Novo/a20190013_004_FKDL202574961-1a-AK2180-AK10751_HK7NYDRXX_L1.fq.gz"
                     ["OHMX20200619_002"]="/data/gerbenm/20200619_Novo_GST_Frame/Novo/a20190013_002_FKDL202574961-1a-AK870-AK10750_HK7NYDRXX_L1.fq.gz"
                     ["OHMX20200619_003"]="/data/gerbenm/20200619_Novo_GST_Frame/Novo/a20190013_003_FKDL202574961-1a-AK1958-AK2941_HK7NYDRXX_L1.fq.gz")

##Input arguments##
echo "INPUT ARUGMENTS:"
ENSEMBL_ANNOT="100"
SPECIES="human"
SPECIES_SHORT="hsa"
ORF="plastid"  											# Mostly standard or plastid
PRICE="Y"
CORES="20"
UNIQUEMAPPING="Y"
CLIPPER="trimmomatic" 									# Mostly trimmomatic or fastx
ADAPTORSEQ="TGGAATTCTCGGGTGCCAAGG"
PHIX="Y"												# Y or N
RRNA="Y"												# Y or N
SNORNA="Y"												# Y or N
TRNA="Y"												# Y or N
TRCOORD="Y"                                             # Y or N
TESTRUN="N"                                             # Y or N
IGENOMESROOT="/data/igenomes/"
ENSEMBLDICT="/share/steven/Ensembl/"
ENSEMBLDB="ENS_${SPECIES_SHORT}_${ENSEMBL_ANNOT}.db"
COMPLOGO='ohmx'											# ohmx or biobix
SCRIPTDIR="/home/gerben/scripts/"
BASEDIR="/data/gerbenm/20200619_Novo_GST_Frame/Novo"
echo "basedir = $BASEDIR"
echo -e "scriptdir = $SCRIPTDIR\n"

echo "Ensembl annotation = $ENSEMBL_ANNOT"
echo "Species = $SPECIES"
echo "Offset calling method = $ORF"
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
#conda activate proteoformer

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
perl $SCRIPTDIR/proteoformer/1_mapping/mapping.pl --inputfile1 $UNZIPFILE --readtype ribo_untr --name $ID --species $SPECIES --ensembl $ENSEMBL_ANNOT --cores $CORES --unique $UNIQUEMAPPING --igenomes_root $IGENOMESROOT --clipper $CLIPPER --adaptor $ADAPTORSEQ --phix $PHIX --rRNA $RRNA --snRNA $SNORNA --tRNA $TRNA --rpf_split N --tr_coord $TRCOORD --price_files $PRICE --suite $ORF --suite_tools_loc $SCRIPTDIR/proteoformer/1_mapping/ > $BASEDIR/$ID/mapping_$ID.txt 2>&1
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
perl $SCRIPTDIR/proteoformer/2_mappingQC_py3/mappingQC.pl --samfile $BASEDIR/$ID/STAR/fastq1/untreat.sam --treated untreated --testrun $TESTRUN --cores $CORES --result_db $BASEDIR/$ID/SQLite/results.db --unique $UNIQUEMAPPING --ens_db $ENSEMBLDICT/$ENSEMBLDB --offset $ORF  --offset_img $BASEDIR/$ID/plastid/${ID}_untreated_p_offsets.png --tool_dir $SCRIPTDIR/proteoformer/2_mappingQC_py3/mqc_tools/ --plotrpftool pyplot3D --output_folder $BASEDIR/mqc/$ID --suppl_out_folder $BASEDIR/mqc_suppl/$ID --html $BASEDIR/mqc/mqc_$ID.html --zip $BASEDIR/mqc/mqc_$ID.zip --tmp $BASEDIR/$ID/tmp/ --comp_logo $COMPLOGO > $BASEDIR/$ID/mQC_$ID.txt 2>&1
rm -rf $BASEDIR/$ID/tmp/mappingqc_untreated
echo -e "mQC performed for $ID \n"
echo -e "\n"

rm -rf  fastq/fastq1_nophix.fq fastq/nophix/  fastq/fastq1_clip* fastq/Unmapped.out.mate1  fastq/Log.* fastq/SJ.out.tab  fastq/Aligned.out.sam fastq/fastq1_norrna.fq fastq/fastq1_norrna_nosnrna.fq fastq/fastq1_norrna_nosnrna_notrna.fq

echo -e "5) Transcript calling $ID \n"
perl $SCRIPTDIR/proteoformer/3_tr_calling/Rule-based/ribo_translation.pl --in_sqlite $BASEDIR/$ID/SQLite/results.db --out_sqlite $BASEDIR/$ID/SQLite/results.db --ens_db $ENSEMBLDICT/$ENSEMBLDB > $BASEDIR/$ID/tr_translation_$ID.txt 2>&1
echo -e "Transcript calling performed on $ID \n"


##To run price, another ENV needs to be loaded.
conda deactivate
conda activate price

echo -e "6) PRICE $ID \n"
python $SCRIPTDIR/proteoformer/4_ORF_calling/using_PRICE/PRICE.py -r $BASEDIR/$ID/SQLite/results.db > $BASEDIR/$ID/PRICE_orf_$ID.txt 2>&1
echo -e "PRICE ORF calling performed on $ID \n"

conda deactivate
conda activate proteoformer


##To run spectre, another ENV needs to be loaded.
conda deactivate
conda activate spectre

echo -e "7) SPECTRE $ID \n"
python $SCRIPTDIR/proteoformer/4_ORF_calling/using_SPECtre/SPECtre.py -r $BASEDIR/$ID/SQLite/results.db -o 28:12,29:12,30:12 -c 60 -x 3 > $BASEDIR/$ID/SPECtre_orf_$ID.txt 2>&1
echo -e "SPECTRE ORF calling performed on $ID \n"

conda deactivate
source activate proteoformer

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

cd $BASEDIR/fastqc_raw
multiqc -i "Raw Read Data" .
cd $BASEDIR

cd $BASEDIR/fastqc_mapped
multiqc -i "Mapped Read Data" .
cd $BASEDIR

echo -e "MultiQC on all raw/mapped data performed\n"

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

