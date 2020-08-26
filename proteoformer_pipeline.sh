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
declare -A datasets=( ["20200001_001"]="/data2/steven/OHMX/OHMX20200001_bimini/raw/RIBO_X204SC20061628-Z01-F001/raw_data/a20200001_001/a20200001_001_FKDL202589086-1a-AK1912-AK1031_HLTJHDRXX_L1.fq.gz" ["20200001_002"]="/data2/steven/OHMX/OHMX20200001_bimini/raw/RIBO_X204SC20061628-Z01-F001/raw_data/a20200001_002/a20200001_002_FKDL202589086-1a-AK8599-AK10750_HLTJHDRXX_L1.fq.gz" ["20200001_003"]="/data2/steven/OHMX/OHMX20200001_bimini/raw/RIBO_X204SC20061628-Z01-F001/raw_data/a20200001_003/a20200001_003_FKDL202589086-1a-AK393-AK2941_HLTJHDRXX_L1.fq.gz" ["20200001_004"]="/data2/steven/OHMX/OHMX20200001_bimini/raw/RIBO_X204SC20061628-Z01-F001/raw_data/a20200001_004/a20200001_004_FKDL202589086-1a-AK8600-AK10751_HLTJHDRXX_L1.fq.gz")

##Input arguments##
echo "INPUT ARUGMENTS:"
ENSEMBL_ANNOT="100"
SPECIES="human"
SPECIES_SHORT="hsa"
ORF="plastid"
PRICE="Y"
CORES="20"
UNIQUEMAPPING="Y"
ADAPTORSEQ="TGGAATTCTCGGGTGCCAAGG"
IGENOMESROOT="/data/igenomes/"
ENSEMBLDICT="/share/steven/Ensembl/"
ENSEMBLDB="ENS_${SPECIES_SHORT}_${ENSEMBL_ANNOT}.db"

SCRIPTDIR="/home/steven/"
BASEDIR="/data2/steven/OHMX/OHMX20200001_bimini/"
echo "basedir = $BASEDIR"
echo -e "scriptdir = $SCRIPTDIR\n"

echo "Ensembl annotation = $ENSEMBL_ANNOT"
echo "Species = $SPECIES"
echo "Offset calling method = $ORF"
echo "PRICE-adapted mapping files = $PRICE"
echo "Cores = $CORES"
echo "Unique mapping = $UNIQUEMAPPING"
echo "Adaptor sequence = $ADAPTORSEQ"
echo "iGenomes root folder = $IGENOMESROOT"
echo "Ensembl DB location = ${ENSEMBLDICT}${ENSEMBLDB}" 
echo -e "\n\n"

##Activate PROTEOFORMER Conda Environment##
source activate proteoformer

##Reference information##
echo -e "Download reference info\n\n"
python $SCRIPTDIR/proteoformer/Additional_tools/ENS_db.py -v $ENSEMBL_ANNOT -s $SPECIES
chmod 755 $ENSEMBLDB
mv $ENSEMBLDB $ENSEMBLDICT
python $SCRIPTDIR/proteoformer/Additional_tools/get_igenomes.py -v $ENSEMBL_ANNOT -s $SPECIES -d $IGENOMESROOT -r -c 15


##Internal dict structure##
mkdir $BASEDIR/fastqc_raw
mkdir $BASEDIR/fastqc_mapped
mkdir $BASEDIR/statistics
mkdir $BASEDIR/mqc



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
    gunzip $FILE
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
perl $SCRIPTDIR/proteoformer/1_mapping/mapping.pl --inputfile1 $UNZIPFILE --readtype ribo_untr --name $ID --species $SPECIES --ensembl $ENSEMBL_ANNOT --cores $CORES --unique $UNIQUEMAPPING --igenomes_root $IGENOMESROOT --clipper fastx --adaptor $ADAPTORSEQ --phix Y --rRNA Y --snRNA Y --tRNA Y --rpf_split N --price_files $PRICE --suite $ORF --suite_tools_loc $SCRIPTDIR/proteoformer/1_mapping/ > $BASEDIR/$ID/mapping_$ID.txt 2>&1
echo -e "Mapping done for $ID \n"

ln -s  $BASEDIR/$ID/STAR/fastq1/untreat.bam $BASEDIR/$ID/$ID.bam

echo -e "3) FastQC mapped file $ID \n"
fastqc $BASEDIR/$ID/$ID.bam -o $BASEDIR/fastqc_mapped -t 20
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
echo -e "4) mQC $ID \n"
perl $SCRIPTDIR/proteoformer/2_mappingQC/mappingQC.pl --samfile $BASEDIR/$ID/STAR/fastq1/untreat.sam --treated untreated --cores $CORES --result_db $BASEDIR/$ID/SQLite/results.db --unique $UNIQUEMAPPING --ens_db $ENSEMBLDICT/$ENSEMBLDB --offset $ORF  --offset_img $BASEDIR/$ID/plastid/${ID}_untreated_p_offsets.png --tool_dir $SCRIPTDIR/proteoformer/2_mappingQC/mqc_tools/ --plotrpftool pyplot3D --output_folder $BASEDIR/mqc/$ID --html $BASEDIR/mqc/mqc_$ID.html --zip $BASEDIR/mqc/mqc_$ID.zip --tmp $BASEDIR/$ID/tmp/ > $BASEDIR/$ID/mQC_$ID.txt 2>&1
rm -rf $BASEDIR/$ID/tmp/mappingqc_untreated
echo -e "mQC performed for $ID \n"
echo -e "\n"

rm -rf  fastq/fastq1_nophix.fq fastq/nophix/  fastq/fastq1_clip* fastq/Unmapped.out.mate1  fastq/Log.* fastq/SJ.out.tab  fastq/Aligned.out.sam fastq/fastq1_norrna.fq fastq/fastq1_norrna_nosnrna.fq fastq/fastq1_norrna_nosnrna_notrna.fq

echo -e "5) Transcript calling $ID \n"
perl $SCRIPTDIR/proteoformer/3_tr_calling/Rule-based/ribo_translation.pl --in_sqlite $BASEDIR/$ID/SQLite/results.db --out_sqlite $BASEDIR/$ID/SQLite/results.db --ens_db $ENSEMBLDICT/$ENSEMBLDB > $BASEDIR/$ID/tr_translation_$ID.txt 2>&1
echo -e "Transcript calling performed on $ID \n"


##To run price, another ENV needs to be loaded.
source deactivate
source acitvate price

echo -e "6) PRICE $ID \n"
python $SCRIPTDIR/proteoformer/4_ORF_calling/using_PRICE/PRICE.py -r $BASEDIR/$ID/SQLite/results.db > $BASEDIR/$ID/PRICE_orf_$ID.txt 2>&1
echo -e "PRICE ORF calling performed on $ID \n"

source deactivate
source activate proteoformer

##Go back to the BaseDir
cd ..

gzip $UNZIPFILE
echo "File $UNZIPFILE zipped"
echo -e "\n\n"

done
#End loop

cd $BASEDIR/fastqc_raw
multiqc -i "Raw Read Data" .
cd $BASEDIR

cd $BASEDIR/fastqc_mapped
multiqc -i "Mapped Read Data" .
cd $BASEDIR

echo -e "MultiQC on all raw/mapped data performed\n"

##Deactivate conda environment
source deactivate

