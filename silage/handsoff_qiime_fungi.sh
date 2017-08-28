#!/bin/bash

# the qiime miniconda environment should already have been activated, if relevant
source /home/mark/miniconda3/bin/activate qiime1

# this script now requires the following... make sure they're on the path
#    sem/parallel... avaialble for OSX?
#    in2csv/csvkit

# I'm not showing steps required to download 16S/ITSreference seqeunce files, etc, provided by Amy via Erica 

# top priorities:  
#    find appropriate flash overlap setting. esp. max, maybe min, too
#    find appropriate min coverage depth for core diversity (doesn't really impact anything "after"?) 
#    especially in context of iys:
#        what to do with unmapped?  
#        pick_open_reference_otus.py and map_reads_to_reference.py?
#    path to Rscript and handsoff_plots.R 
#    fastqc results maynnot all being moved to fastqc folder

# reset this POSIX env var in case getopts has been used previously in this shell
OPTIND=1
# initialize options-related variables if necessary

# usage statement (on -h or if problems with getopts?)
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

# show the contents of the config file
#    plus, show if any additional parapmeters were passed
[ "$1" = "--" ] && shift
echo "CONFIG_FILE=$CONFIG_FILE; unprocessed: $@"

###

# EXPECTED DATA ORGANIZATION & METADATA
# all R1 and R2 files for a given gene 
# from a given project/study should be in one folder
# with nothing else

# MAY NEED TO MOVE SOME MAPPING FIEL COLUMSN AND INSERT BARCODE PLACEHOLDERS

# UNNECESSARY ... look like Erica will be able to share data and metadat fiels as requested
#   if data files aren't organized as above and project/sample metadata isn't complete
#   try making a "fasta" file out of the filenames
#   list2fasta4mafft.py
#   and run through mafft in full text mode to better understand setup

# read project/gene specific information
source ${CONFIG_FILE}

# file\ path.xlsx can be converted to a csv file with in2csv, 
#   and then converted to tab delim in real time as below
# subseqeucnt steps may barf if there is a blank line at the end?
in2csv ${PROJ_GENE_HOME}/file\ path.xlsx | sed 's/,/\t/g' | sed 's/^\s*$//' | cut -f5,6 > ${PROJ_GENE_HOME}/min_file_path.txt

declare -A myArray
while IFS=$'\t' read number fn
do
  myArray[$number]=$fn
done < ${PROJ_GENE_HOME}/min_file_path.txt

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

###
### had been using the mapping file plus a sample metadata fiel that I created
### now use mapping plus ${PROJ_GENE_HOME}/min_file_path.txt ... see above for conversion from xlsx
###

# initial fastqc
# one loop each for r1 and r2...
# verbose, but better thread management?
if [ $INITIAL_FASTQC == "TRUE" ]
  then
    while read p; do
      # make this a function?
      # pass output as env vars
      CURRENT_PAIR=${PROJ_GENE_HOME}/${DATA_DIR}/${p} ;
      CURRENT_R1=${CURRENT_PAIR}${R1}${FILE_EXT} ;
      CURRENT_R2=${CURRENT_PAIR}${R2}${FILE_EXT} ;
      # ls -l ${CURRENT_R1} ;
      # ls -l ${CURRENT_R2} ;
      echo `date` "fastqc'ing ${CURRENT_PAIR}" | tee -a ${LOG_FILE}
      # omitting per-process logging
      # sending jobs to background (poor man's parallelization)
      # tune java heap memory?
      sem --no-notice -j-1 ${FASTQC_PATH} ${CURRENT_R1} ";"
     done < <(awk 'BEGIN {FS="\t"}; NR>1 {print $2}' ${PROJ_GENE_HOME}/min_file_path.txt | sed 's/...$//g' )

     sem --wait


    # harmless error if folder already exists?
    mkdir ${PROJ_GENE_HOME}/${DATA_DIR}/fastqc

    sleep 10
    cp ${PROJ_GENE_HOME}/${DATA_DIR}/*fastqc.* ${PROJ_GENE_HOME}/${DATA_DIR}/fastqc
    sleep 10
    rm -rf ${PROJ_GENE_HOME}/${DATA_DIR}/*fastqc.*
    sleep 10
    # fastqc results may not all being moved to fastqc folder

    while read p; do
      CURRENT_PAIR=${PROJ_GENE_HOME}/${DATA_DIR}/${p} ;
      CURRENT_R1=${CURRENT_PAIR}${R1}${FILE_EXT} ;
      CURRENT_R2=${CURRENT_PAIR}${R2}${FILE_EXT} ;
      echo `date` "fastqc'ing ${CURRENT_PAIR}" | tee -a ${LOG_FILE}
      sem --no-notice -j-1 ${FASTQC_PATH} ${CURRENT_R2} ";"
     done < <(awk 'BEGIN {FS="\t"}; NR>1 {print $2}' ${PROJ_GENE_HOME}/min_file_path.txt | sed 's/...$//g' )

     sem --wait


    # harmless error if folder already exists?
    mkdir ${PROJ_GENE_HOME}/${DATA_DIR}/fastqc

    sleep 10
    cp ${PROJ_GENE_HOME}/${DATA_DIR}/*fastqc.* ${PROJ_GENE_HOME}/${DATA_DIR}/fastqc
    sleep 10
    rm -rf ${PROJ_GENE_HOME}/${DATA_DIR}/*fastqc.*
    sleep 10
    # fastqc results may not all being moved to fastqc folder
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
    CURRENT_PAIR=${PROJ_GENE_HOME}/${DATA_DIR}/${p} ;
    CURRENT_R1=${CURRENT_PAIR}${R1}${FILE_EXT} ;
    CURRENT_R2=${CURRENT_PAIR}${R2}${FILE_EXT} ;
    # ls -l ${CURRENT_R1} ;
    # ls -l ${CURRENT_R2} ;

    echo `date` "flashing ${CURRENT_PAIR}" | tee -a ${LOG_FILE}
    FLASH_COMMAND="${FLASH_PATH} ${CURRENT_R1} ${CURRENT_R2} -M ${FLASH_OVERLAP} -d ${FLASHED_DIR} -o ${p}"
    echo ${FLASH_COMMAND}
    ${FLASH_COMMAND} | tee -a ${LOG_FILE}
    done < <(awk 'BEGIN {FS="\t"}; NR>1 {print $2}' ${PROJ_GENE_HOME}/min_file_path.txt | sed 's/...$//g' )

  mkdir ${EXTENDED_DIR}

  # segregate flash-extended files
  # what followup to perform on unextended reads?  what were they?  something like BLAST?
  mv ${FLASHED_DIR}/*.extendedFrags.fastq ${EXTENDED_DIR}
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
      sem --no-notice -j-1 ${FASTQC_PATH} $f  ";"
    done

    sem --wait

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
    while read p; 
    # for f in ${EXTENDED_DIR}/*.fastq
      do
      # NOPERIODS=`basename $f | sed 's/\./_/g'`

      # make this a function?
      # pass output as env vars

      # echo $p

      # echo ${myArray["$p"]}

      associated=`echo ${myArray["$p"]} | sed 's/...$//g'`

      echo `date` "'splitting' ${associated}" | tee -a ${LOG_FILE}

      # ls -l ${EXTENDED_DIR}/${associated}.extendedFrags.fastq

       # doesn't like being run as an environment variable ?!
     sem --no-notice -j-1  split_libraries_fastq.py -m ${QIIME_MAPPING_FILE} -i ${EXTENDED_DIR}/${associated}.extendedFrags.fastq  --sample_ids ${p} -o ${EXTENDED_DIR}/${associated}/ --barcode_type 'not-barcoded'  ";"
      # omit per-process logging?  doesn't generate any stdout log... saves to file?
      # send jobs to background (poor man's & parallelization)
      # or set number of threads?  by default appears to be single-threaded
      # java heap memory prob irrelevant?
  
     done < <(awk 'BEGIN {FS="\t"}; NR>1 {print $1}' ${PROJ_GENE_HOME}/min_file_path.txt )
  
    sem --wait

fi


if [ $DO_QIIME_PICK == "TRUE" ]
  then
    find ${EXTENDED_DIR} -name seqs.fna -exec cat {} \; > ${EXTENDED_DIR}/flash_trim_cat.fna

    echo `date` "start OTU picking" | tee -a ${LOG_FILE}

    # 16S
    # pick_open_reference_otus.py
    # -i ./14All.fna 
    # -m uclust 
    # -o./14All_otus 
    # -p ./otu_pick_para.txt
    # -r /Users/amybiddle/Downloads/gg_13_8_otus/rep_set/97_otus.fasta 

    # its... implict -m uclust?
    # pick_open_reference_otus.py \
    # -i its-soils-tutorial/seqs.fna \
    # -o otus/ \
    # -p its-soils-tutorial/params.txt \
    # -r its_12_11_otus/rep_set/97_otus.fasta \
    # --suppress_align_and_tree

    # -a = parllel
    # -f = force overwrite
    pick_open_reference_otus.py \
        -a \
        -f \
        -i ${EXTENDED_DIR}/flash_trim_cat.fna \
        -m uclust \
        -o ${PICKED_DIR} \
        -O ${THREAD_CT} \
        -p ${OTU_PICK_PARAM_FILE} \
        -r ${REF_SEQ_FILE} \
        --suppress_align_and_tree | tee -a ${LOG_FILE} 
  
    echo `date` "OTU picking complete" | tee -a ${LOG_FILE}

    # using force, make that optional
    # running 6 threads in parallel if possible, make that optyional
    # took ~ 9 minutes for the understudy 16s

    ## skip this?
    # map_reads_to_reference.py -i allfile.fna -o all_otu -r 99_otus.fasta -t 99_otu_taxonomy.txt -m usearch
fi


# fast
biom summarize-table -i ${PICKED_DIR}/otu_table_mc2_w_tax.biom | tee -a ${LOG_FILE}

if [ $DO_PYTHON_OTU_FILTERING == "TRUE" ]
then
# make these settings parameters
MIN_CT_FILTERED_FILE=${PICKED_DIR}/otu_table_no_sing_doub.biom
filter_otus_from_otu_table.py -i ${PICKED_DIR}/otu_table_mc2_w_tax.biom -o ${MIN_CT_FILTERED_FILE} -n 3 | tee -a ${LOG_FILE}
fi


if [ $DO_BIOM_CONVERSION == "TRUE" ]
then
biom convert \
-i ${MIN_CT_FILTERED_FILE} \
-o ${PICKED_DIR}/otu_table_json.biom \
--to-json \
--table-type="OTU table" | tee -a ${LOG_FILE} 
fi


if [ $DO_CORE_DIV == "TRUE" ]
then
echo `date` "start core diversity" | tee -a ${LOG_FILE}
# watch out for large negative eigenvalues
# 16S (Amy didn't muilti-thread, compare by treatment
#    DID surpress beta diversity analysis because she used a small input set
rm -rf ${PICKED_DIR}/corediv
core_diversity_analyses.py \
    -a \
    -c treatment \
    -e ${CORE_DIV_MIN_DEPTH} \
    -i ${MIN_CT_FILTERED_FILE} \
    -m ${QIIME_MAPPING_FILE} \
    -o ${PICKED_DIR}/corediv \
    -O ${THREAD_CT}  \
    --nonphylogenetic_diversity | tee -a ${LOG_FILE} 
# its
# core_diversity_analyses.py \
# -e 353 \
# -i otus/otu_table_no_sing_doub.biom \
# -m its-soils-tutorial/map.txt \
# -o cdout/ \
# --nonphylogenetic_diversity
echo `date` "core diversity complete" | tee -a ${LOG_FILE}
fi


if [ $DO_PHYLOSEQ == "TRUE" ]
then
# now start phyloseq analysis in R
#     just run as a background sript and write to PDF device,
#     or kint into a notebook document (Word, HTML or PDF, slides, more!)
#     could even puth the bash portions into the markdown... harder to control threading/forking?
echo `date` "start phyloseq analysis and plotting in R" | tee -a ${LOG_FILE}
Rscript /home/mark/gitstage/uderica/silage/handsoff_plots.R \
    --otufile=${PICKED_DIR}/otu_table_json.biom \
    --trefile=${PICKED_DIR}/rep_set.tre \
    --mapfile=${QIIME_MAPPING_FILE} \
    --prefix=${PHYLOSEQ_OUTPUT_PREFIX} >> ${LOG_FILE} 2>&1 
echo `date` "phyloseq complete" | tee -a ${LOG_FILE}
fi

