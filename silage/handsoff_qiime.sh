#!/bin/bash

# TODO add usage statement (on -h or if problems with getopts?)

# reset this POSIX env var in case getopts has been used previously in the shell
OPTIND=1

# initialize options-related variables if necessary

function show_usage_help()
{
    echo "usage: $0 -c <configuration file>"
}

while getopts "h?c:" opt; do
    case "$opt" in
    h|\?)
        show_usage_help
        exit 0
        ;;
    c)  CONFIG_FILE=$OPTARG
        ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

echo "CONFIG_FILE=$CONFIG_FILE; unprocessed: $@"

###

# all R1 and R2 files for a given gene 
# from a given project/study should be in one folder
# with nothing else

# if data files aren't organized above and project/sample metadata isn't complete
# try making a "fasta" file out of the filenames
# list2fasta4mafft.py
# and run through mafft in full text mode to better understand setup

# the qiime miniconda environment should already have been activated, if relevant

# read project/gene specific information
source ${CONFIG_FILE}

# start with a fresh log file?
# harmless error if file doesn't exist yet?
if [ $FRESH_LOG == "TRUE" ]
  then
    rm ${LOG_FILE}
fi

echo "configuration from $CONFIG_FILE" >> ${LOG_FILE}

cat ${CONFIG_FILE} >> ${LOG_FILE}


# the mapping file should already have been validated
# I can create sequential barcodes if necessary
# list2fasta4mafft.py
#  merge the mapping file with the experimental metadata?  see expectations below
validate_mapping_file.py -m ${QIIME_MAPPING_FILE} -o ${QIIME_MAPPING_FILE}_mapcheck


# initial fastqc
if [ $INITIAL_FASTQC == "TRUE" ]
  then
    while read p; do
      # make this a function?
      # pass output as env vars
      CURRENT_PAIR=${PROJ_GENE_HOME}/${DATA_DIR}/${PAIR_PREFIX}${p}${PAIR_SUFFIX} ;
      CURRENT_R1=${CURRENT_PAIR}${R1}${FILE_EXT} ;
      CURRENT_R2=${CURRENT_PAIR}${R2}${FILE_EXT} ;
      ls -l ${CURRENT_R1} ;
      ls -l ${CURRENT_R2} ;
      echo `date` "fastqc'ing ${CURRENT_PAIR}" | tee -a ${LOG_FILE}
      # omitting per-process logging
      # sending jobs to background (poor man's parallelization)
      # tune java heap memory?
      ${FASTQC_PATH} ${CURRENT_R1} &
      ${FASTQC_PATH} ${CURRENT_R2} &
#    done < <(awk 'BEGIN {FS="\t"}; NR>2 {print $5}' ${SAMPLE_METADATA_FILE})
    done < <(awk 'BEGIN {FS="\t"}; NR>2 {print $1}' ${QIIME_MAPPING_FILE})
    wait

    # harmless error if folder already exists?
    mkdir ${PROJ_GENE_HOME}/${DATA_DIR}/fastqc

    mv ${PROJ_GENE_HOME}/${DATA_DIR}/*fastqc.* ${PROJ_GENE_HOME}/${DATA_DIR}/fastqc
fi


# multiqc
if [ $INITIAL_MULTIQC == "TRUE" ]
  then
    MULTIQC_COMMAND="${MULTIQC_PATH} ${PROJ_GENE_HOME}/${DATA_DIR}/fastqc -o ${PROJ_GENE_HOME}/${DATA_DIR}/multiqc"
    echo ${MULTIQC_COMMAND}
    ${MULTIQC_COMMAND}
    # parse/analyze multiqc? take action if poor quality, like aborting remainder of script?
fi


# flash... already multithreaded?
# extend overlapping reads (since read len > amplicon len / 2)
if [ $DO_FLASH == "TRUE" ]
  then
  mkdir ${FLASHED_DIR}
  while read p; do
    # make this a function?
    CURRENT_PAIR=${PROJ_GENE_HOME}/${DATA_DIR}/${PAIR_PREFIX}${p}${PAIR_SUFFIX} ;
    CURRENT_R1=${CURRENT_PAIR}${R1}${FILE_EXT} ;
    CURRENT_R2=${CURRENT_PAIR}${R2}${FILE_EXT} ;
    ls -l ${CURRENT_R1} ;
    ls -l ${CURRENT_R2} ;

    echo `date` "flashing ${CURRENT_PAIR}" | tee -a ${LOG_FILE}
    FLASH_COMMAND="${FLASH_PATH} ${CURRENT_R1} ${CURRENT_R2} -M ${FLASH_OVERLAP} -d ${FLASHED_DIR} -o sample_${p}"
    echo ${FLASH_COMMAND}
    ${FLASH_COMMAND} | tee -a ${LOG_FILE}

  done < <(awk 'BEGIN {FS="\t"}; NR>2 {print $5}' ${SAMPLE_METADATA_FILE})

  mkdir ${EXTENDED_DIR}

  # segregate flash-extended files
  # what followup to perform on unextended reads?  what were they?  something like BLAST?
  mv ${FLASHED_DIR}/sample*.extendedFrags.fastq ${EXTENDED_DIR}
fi


###  add flags & if statements to determine which steps to execute
# post-flash fastqc
# not consistenet here in terms of logging etc.
if [ $POST_FLASH_BOTHQC == "TRUE" ]
  then
    for f in ${EXTENDED_DIR}/*.fastq
    do
      echo `date` "fastqc'ing $f" | tee -a ${LOG_FILE}
      echo "${FASTQC_PATH} $f" 
      ${FASTQC_PATH} $f &
    done

    wait


    # post-flash multiqc
    # harmless error if folder already exists?
    mkdir ${EXTENDED_DIR}/fastqc

    mv ${EXTENDED_DIR}/*fastqc.* ${EXTENDED_DIR}/fastqc

    # multiqc

    MULTIQC_COMMAND="${MULTIQC_PATH} ${EXTENDED_DIR}/fastqc -o ${EXTENDED_DIR}/multiqc"
    echo ${MULTIQC_COMMAND}

    ${MULTIQC_COMMAND}
fi


###

# QIIME SPLIT (& CLEANUP ?)
if [ $DO_QIIME_SPLIT == "TRUE" ]
  then
    while read p; do
      # make this a function?
      # pass output as env vars

      echo `date` "'splitting' ${p}" | tee -a ${LOG_FILE}

      # expectation: there should be a one-to-one match between sample IDs in the QIIME mapping file
      #  and ${SAMPLE_METADATA_FILE}
      #  haven't doen anything to test that yet
 
      SPLIT_CMD="split_libraries_fastq.py  -m ${QIIME_MAPPING_FILE} -i ${EXTENDED_DIR}/sample_${p}.extendedFrags.fastq --sample_ids ${p} -o ${EXTENDED_DIR}/sample_${p}/ --barcode_type 'not-barcoded' "
      echo ${SPLIT_CMD}
      # ${SPLIT_CMD} 
      # doesn't like being run as an environment variable ?!

split_libraries_fastq.py -m ${QIIME_MAPPING_FILE} -i ${EXTENDED_DIR}/sample_${p}.extendedFrags.fastq --sample_ids ${p} -o ${EXTENDED_DIR}/sample_${p}/ --barcode_type 'not-barcoded' &

      # omit per-process logging?  doesn't generate any stdout log... saves to file?
      # send jobs to background (poor man's & parallelization)
      # or set number of threads?  by default appears to be single-threaded
      # java heap memory prob irrelevant?
  
    done < <(awk 'BEGIN {FS="\t"}; NR>2 {print $5}' ${SAMPLE_METADATA_FILE})
  
    wait

fi


if [ $DO_QIIME_PICK == "TRUE" ]
  then
    find ${EXTENDED_DIR} -name seqs.fna -exec cat {} \; > ${EXTENDED_DIR}/flash_trim_cat.fna

    echo `date` "start OTU picking" | tee -a ${LOG_FILE}

    # pick_open_reference_otus.py -i ./14All.fna -o./14All_otus -r /Users/amybiddle/Downloads/gg_13_8_otus/rep_set/97_otus.fasta -m uclust -p ./otu_pick_para.txt
    pick_open_reference_otus.py -i ${EXTENDED_DIR}/flash_trim_cat.fna -o ${PICKED_DIR} -r ${REF_SEQ_FILE} -m uclust -p ${OTU_PICK_PARAM_FILE} -a -O ${THREAD_CT} -f

    echo `date` "OTU picking complete" | tee -a ${LOG_FILE}

    # using force, make that optional
    # running 6 threads in parallel if possible, make that optyional
    # took ~ 9 minutes for the understudy 16s
fi


## skip this?
# map_reads_to_reference.py -i allfile.fna -o all_otu -r 99_otus.fasta -t 99_otu_taxonomy.txt -m usearch


# fast
biom summarize-table -i ${PICKED_DIR}/otu_table_mc2_w_tax_no_pynast_failures.biom | tee -a ${LOG_FILE}


# make these settings parameters
MIN_CT_FILTERED_FILE=${PICKED_DIR}/otu_table_no_sing_doub.biom
filter_otus_from_otu_table.py -i ${PICKED_DIR}/otu_table_mc2_w_tax_no_pynast_failures.biom -o ${MIN_CT_FILTERED_FILE} -n 3 | tee -a ${LOG_FILE}

biom convert \
-i ${MIN_CT_FILTERED_FILE} \
-o ${PICKED_DIR}/otu_table_json.biom \
--to-json \
--table-type="OTU table"

echo `date` "start core diversity" | tee -a ${LOG_FILE}
# watch out for large negative eigenvalues
core_diversity_analyses.py -i ${MIN_CT_FILTERED_FILE} -m ${QIIME_MAPPING_FILE} -o ${PICKED_DIR}/corediv -e ${CORE_DIV_MIN_DEPTH} -t ${PICKED_DIR}/rep_set.tre  -c treatment -a -O ${THREAD_CT}  | tee -a ${LOG_FILE}
echo `date` "core diversity complete" | tee -a ${LOG_FILE}

# now start phyloseq analysis in R
# just run as a background sript and write to PDF device,
# or kint into a notebook document (Word, HTML or PDF)







