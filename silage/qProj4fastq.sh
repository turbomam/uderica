# . /home/mark/gitstage/uderica/silage/qProj4fastq.sh

source /home/mark/miniconda3/bin/activate qiime1

export IN2CSV_PATH="/usr/local/bin/in2csv"
export FASTQC_PATH="/home/mark/gitstage/uderica/silage/FastQC/fastqc"
export MULTIQC_PATH="/usr/local/bin/multiqc"

export MULTIQC_CONFIG="/home/mark/gitstage/uderica/silage/multiqc_config.yaml"

export THREAD_COUNT="7"

export rootFold="/home/mark/gitstage/uderica/silage"

export projName=`echo $projFold | sed -r 's/\s+/_/g' | sed 's/^\///'`

# this is WAY easier and more robust if there are no space in the file and directory names

export projFold="/15_Hil_CS"
export bactFold="/Bacteria"
export bactDataFold="/16S"
export bactKeyFile="15 Hil CS bacteria file path.xlsx"
export fungiFold="/Fungi"
export fungiDataFold="/ITS"
export fungiKeyFile="15 Hil CS yeast file path.xlsx"

export sampleIdCol="5"
export R1Col="6"
export R2Col="7"

export fastqFileSuffix=".fastq.gz"

export FLASH_PATH="/home/mark/gitstage/uderica/silage/FLASH-1.2.11/flash"
export MIN_FLASH_OVERLAP=50
export MAX_FLASH_OVERLAP=250

# export QIIME_MAPPING_FILE="${PROJ_GENE_HOME}/map.txt"
# haven't expliclitly checked the mapping file yet

export REF_SEQ_FILE="/home/mark/gitstage/uderica/silage/gg_13_8_otus/rep_set/97_otus.fasta"

export BACT_ALIAS="${rootFold}${projFold}${bactFold}${bactDataFold}"

export CORE_DIV_MIN_DEPTH=10000




${IN2CSV_PATH} --no-inference "${rootFold}${projFold}${bactFold}/${bactKeyFile}" > ${projName}.csv

sed -e 's/,/\t/g' ${projName}.csv > ${projName}.txt

awk 'BEGIN {FS="\t"}; NR>1 {print $5 "\t" $6 "\t" $7 }' ${projName}.txt > ${projName}_min.txt

while read p; 
	do 
	arr=($p) ; 
        ls -l "${BACT_ALIAS}/${arr[1]}${fastqFileSuffix}"
	${FASTQC_PATH} -t ${THREAD_COUNT} "${BACT_ALIAS}/${arr[1]}${fastqFileSuffix}" "${BACT_ALIAS}/${arr[2]}${fastqFileSuffix}"
done < 15_Hil_CS_min.txt

${MULTIQC_PATH} --config ${MULTIQC_CONFIG} "${BACT_ALIAS}" -o "${BACT_ALIAS}"

rm -rf "${BACT_ALIAS}/extendedFrags/"
mkdir -p "${BACT_ALIAS}/extendedFrags/"

  while read p; do
	arr=($p) ; 
    echo `date` "flashing ${BACT_ALIAS}/${arr[1]}${fastqFileSuffix} and ${BACT_ALIAS}/${arr[2]}${fastqFileSuffix}"

${FLASH_PATH} "${BACT_ALIAS}/${arr[1]}${fastqFileSuffix}" "${BACT_ALIAS}/${arr[2]}${fastqFileSuffix}" -m ${MIN_FLASH_OVERLAP} -M ${MAX_FLASH_OVERLAP} -d "${BACT_ALIAS}/extendedFrags/"  -o ${arr[0]}

    done < 15_Hil_CS_min.txt

# in2csv may cast ints to floats, which may mean that the sample names end in .0
${FASTQC_PATH} -t ${THREAD_COUNT}  "${BACT_ALIAS}/extendedFrags"/*.extendedFrags.fastq

${MULTIQC_PATH} --config ${MULTIQC_CONFIG} "${BACT_ALIAS}/extendedFrags/" -o "${BACT_ALIAS}/extendedFrags/"

# rm -rf *fastqc* ; rm -rf *extended* ; rm -rf *hist* ; rm -rf *Combined* ; rm -rf *multi* ; rm -rf pick ; rm -rf split



rm -rf "${BACT_ALIAS}"/split/
mkdir -p "${BACT_ALIAS}"/split/

    while read p; 
	do 
	arr=($p) 
SAMPLE_RAW=${arr[0]}
SAMPLE_INT=`echo $SAMPLE_RAW | sed -r 's/\..*$//'`

echo "splitting $SAMPLE_INT"

       split_libraries_fastq.py -v \
              -q 29 \
              -m "${rootFold}${projFold}${bactFold}"/map.txt \
              -i "${BACT_ALIAS}"/extendedFrags/${SAMPLE_RAW}.extendedFrags.fastq  \
              --sample_ids ${SAMPLE_INT} \
              -o "${BACT_ALIAS}"/split/${SAMPLE_INT} \
              --barcode_type 'not-barcoded'

 
    done < 15_Hil_CS_min.txt

rm -rf "${BACT_ALIAS}"/pick
mkdir "${BACT_ALIAS}"/pick


    sleep 15

    find "${BACT_ALIAS}/split/" -name seqs.fna -exec cat {} \; > "${BACT_ALIAS}"/pick/flash_pick_cat.fna

    sleep 15


    echo `date` "start OTU picking" | tee -a ${LOG_FILE}

    # -a = parllel
    # -f = force overwrite
    pick_open_reference_otus.py \
        -a \
        -f \
        -i "${BACT_ALIAS}"/pick/flash_pick_cat.fna \
        -m uclust \
        -o "${BACT_ALIAS}"/pick/ \
        -O ${THREAD_COUNT} \
        -p ${rootFold}/bacteria_otu_pick_param.txt \
        -r ${REF_SEQ_FILE}

    echo `date` "OTU picking complete" | tee -a ${LOG_FILE}


biom summarize-table -i "${BACT_ALIAS}"/pick/otu_table_mc2_w_tax_no_pynast_failures.biom

filter_otus_from_otu_table.py \
-n 3 \
-i "${BACT_ALIAS}"/pick/otu_table_mc2_w_tax_no_pynast_failures.biom \
-o "${BACT_ALIAS}"/pick/otu_table_3plus.biom

rm -rf "${BACT_ALIAS}"/corediv

# get mean counts and sd from biom_table_summary.txt, set min to mean - ( 2 x sd )

core_diversity_analyses.py \
    -a \
    -c Treatment \
    -e ${CORE_DIV_MIN_DEPTH} \
    -i "${BACT_ALIAS}"/pick/otu_table_3plus.biom \
    -m "${rootFold}${projFold}${bactFold}"/map.txt \
    -o "${BACT_ALIAS}"/corediv \
    -O ${THREAD_COUNT} \
    -t "${BACT_ALIAS}"/pick/rep_set.tre

