# . /home/mark/gitstage/uderica/silage/qProj4fastq.sh

source /home/mark/miniconda3/bin/activate qiime1

# source the project's envvars script
# reuse getopts approach from previous "handsoff" script
source $1

${IN2CSV_PATH} --no-inference "${rootFold}${projFold}${domainFold}/${filePathFile}" > "${rootFold}${projFold}${domainFold}/${projName}".csv

sed -e 's/,/\t/g' "${rootFold}${projFold}${domainFold}/${projName}.csv" > "${rootFold}${projFold}${domainFold}/${projName}.txt"

awk 'BEGIN {FS="\t"}; NR>1 {print $5 "\t" $6 "\t" $7 }' "${rootFold}${projFold}${domainFold}/${projName}.txt" > "${rootFold}${projFold}${domainFold}/${projName}_min.txt"

while read p; 
	do 
	arr=($p) ; 

ls -l "${BACT_ALIAS}/${arr[1]}${fastqFileSuffix}"
ls -l "${BACT_ALIAS}/${arr[2]}${fastqFileSuffix}"

	${FASTQC_PATH} -t ${THREAD_COUNT} "${BACT_ALIAS}/${arr[1]}${fastqFileSuffix}" "${BACT_ALIAS}/${arr[2]}${fastqFileSuffix}"
done < "${rootFold}${projFold}${domainFold}/${projName}_min.txt"

${MULTIQC_PATH} --config ${MULTIQC_CONFIG} "${BACT_ALIAS}" -o "${BACT_ALIAS}"

rm -rf "${BACT_ALIAS}/extendedFrags/"
mkdir -p "${BACT_ALIAS}/extendedFrags/"

  while read p; do
	arr=($p) ; 
    echo `date` "flashing ${BACT_ALIAS}/${arr[1]}${fastqFileSuffix} and ${BACT_ALIAS}/${arr[2]}${fastqFileSuffix}"

${FLASH_PATH} "${BACT_ALIAS}/${arr[1]}${fastqFileSuffix}" "${BACT_ALIAS}/${arr[2]}${fastqFileSuffix}" -m ${MIN_FLASH_OVERLAP} -M ${MAX_FLASH_OVERLAP} -d "${BACT_ALIAS}/extendedFrags/"  -o ${arr[0]}

    done < "${rootFold}${projFold}${domainFold}/${projName}_min.txt"

# in2csv may cast ints to floats, which may mean that the sample names end in .0
${FASTQC_PATH} -t ${THREAD_COUNT}  "${BACT_ALIAS}/extendedFrags"/*.extendedFrags.fastq

${MULTIQC_PATH} --config ${MULTIQC_CONFIG} "${BACT_ALIAS}/extendedFrags/" -o "${BACT_ALIAS}/extendedFrags/"

# rm -rf *fastqc* ; rm -rf *extended* ; rm -rf *hist* ; rm -rf *Combined* ; rm -rf *multi* ; rm -rf pick ; rm -rf split ; rm -rf corediv

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
              -m "${rootFold}${projFold}${domainFold}"/map.txt \
              -i "${BACT_ALIAS}"/extendedFrags/${SAMPLE_RAW}.extendedFrags.fastq  \
              --sample_ids ${SAMPLE_INT} \
              -o "${BACT_ALIAS}"/split/${SAMPLE_INT} \
              --barcode_type 'not-barcoded'

 
    done < "${rootFold}${projFold}${domainFold}/${projName}_min.txt"

rm -rf "${BACT_ALIAS}"/pick
mkdir "${BACT_ALIAS}"/pick


    sleep 15

    find "${BACT_ALIAS}/split/" -name seqs.fna -exec cat {} \; > "${BACT_ALIAS}"/pick/flash_pick_cat.fna

    sleep 15


echo `date` "start OTU picking"


# elif fungi...
if [ $lifeDomain == "bacteria" ]
then
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
elif [ $lifeDomain == "fungi" ]
then
    pick_open_reference_otus.py \
        -a \
        -f \
        -i "${BACT_ALIAS}"/pick/flash_pick_cat.fna \
        -m uclust \
        -o "${BACT_ALIAS}"/pick/ \
        -O ${THREAD_COUNT} \
        -p ${rootFold}/fungi_otu_pick_param.txt \
        -r ${REF_SEQ_FILE} \
        --suppress_align_and_tree
fi

    echo `date` "OTU picking complete"


biom summarize-table -i "${BACT_ALIAS}"/pick/otu_table_mc2_w_tax_no_pynast_failures.biom

filter_taxa_from_otu_table.py -i "${BACT_ALIAS}"/pick/otu_table_mc2_w_tax_no_pynast_failures.biom -o "${BACT_ALIAS}"/pick/picked_no_chloro.biom -n c__Chloroplast

biom summarize-table -i "${BACT_ALIAS}"/pick/picked_no_chloro.biom

filter_otus_from_otu_table.py \
-n 3 \
-i "${BACT_ALIAS}"/pick/picked_no_chloro.biom \
-o "${BACT_ALIAS}"/pick/otu_table_3plus.biom

rm -rf "${BACT_ALIAS}"/corediv

# get mean counts and sd from biom_table_summary.txt, set min to mean - ( 2 x sd )

# Min: 27347.0
# Max: 99927.0
# Median: 54502.000
# Mean: 55241.407
# Std. dev.: 16687.858




# elif fungi...
if [ $lifeDomain == "bacteria" ]
then
core_diversity_analyses.py \
    -a \
    -c ${CORE_DIV_EXP_FACTOR} \
    -e ${CORE_DIV_MIN_DEPTH} \
    -i "${BACT_ALIAS}"/pick/otu_table_3plus.biom \
    -m "${rootFold}${projFold}${domainFold}"/map.txt \
    -o "${BACT_ALIAS}"/corediv \
    -O ${THREAD_COUNT} \
    -t "${BACT_ALIAS}"/pick/rep_set.tre
elif [ $lifeDomain == "fungi" ]
then
core_diversity_analyses.py \
    -a \
    -c ${CORE_DIV_EXP_FACTOR} \
    -e ${CORE_DIV_MIN_DEPTH} \
    -i "${BACT_ALIAS}"/pick/otu_table_3plus.biom \
    -m "${rootFold}${projFold}${domainFold}"/map.txt \
    -o "${BACT_ALIAS}"/corediv \
    -O ${THREAD_COUNT} \
    --nonphylogenetic_diversity
fi

