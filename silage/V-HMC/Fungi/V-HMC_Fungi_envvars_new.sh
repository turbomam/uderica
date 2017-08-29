export PROJ_GENE_TITLE="V-HMC_Fungi"
export PROJ_GENE_HOME="/home/mark/gitstage/uderica/silage/V-HMC/Fungi"

export R1="_R1"
export R2="_R2"
export FILE_EXT=".fastq"

export LOG_FILE="${PROJ_GENE_HOME}/${PROJ_GENE_TITLE}.log"

export REF_SEQ_FILE="/home/mark/gitstage/uderica/silage/its_12_11_otus/rep_set/97_otus.fasta"

export FASTQC_PATH="/home/mark/gitstage/uderica/silage/FastQC/fastqc"
export MULTIQC_PATH="multiqc"

export MULTIQC_CONFIG="/home/mark/gitstage/uderica/silage/multiqc_config.yaml"

# also sem (parallel)... avaialble for OSX? and in2csv from csvkit

export FLASH_PATH="/home/mark/gitstage/uderica/silage/FLASH-1.2.11/flash"
export MIN_FLASH_OVERLAP=50
export MAX_FLASH_OVERLAP=250

# work on this some more... how to pick best depth?
export CORE_DIV_MIN_DEPTH=10000
export THREAD_CT=6

export SLEEP_SECS=15

export DATA_DIR="raw_data_fungi"

export FLASHED_DIR="${PROJ_GENE_HOME}/${DATA_DIR}/flashed"
export EXTENDED_DIR="${FLASHED_DIR}/extended"
export PICKED_DIR="${EXTENDED_DIR}/flash_trim_cat_pick"

export OTU_PICK_PARAM_FILE="${PROJ_GENE_HOME}/otu_pick_param.txt"

# look here for sample ID, barcodes, primers, optionally experimental factors
# export QIIME_MAPPING_FILE="${PROJ_GENE_HOME}/map16S_MAM_some_factors.txt"
export QIIME_MAPPING_FILE="${PROJ_GENE_HOME}/map.txt"

# look here for sample infomration, except mapping info (esp. barcode and primers)
# will be looping though some column looking for IDs,
# which should then appear in files with names:
# ${PROJ_GENE_HOME}/${DATA_DIR}/${PAIR_PREFIX}<ID>${PAIR_SUFFIX}{${R1}|${R2}}${FILE_EXT}
# export SAMPLE_METADATA_FILE="${PROJ_GENE_HOME}/LalStressMicrobiomeInfo.txt"

export FILE_PATH_FILE="${PROJ_GENE_HOME}/file_path.txt"

# could probably just loop though mapping file instead
# add capability to look for ID:filename relationships in some lookup table

# not doing any error/dependency checking yet
# not doing and cleanup from previous run yet
# FALSE
export FRESH_LOG="TRUE"
export INITIAL_FASTQC="TRUE"
export INITIAL_MULTIQC="TRUE"
export DO_FLASH="TRUE"
export POST_FLASH_BOTHQC="TRUE"
# export POST_FLASH_FASTQC="TRUE"
# export POST_FLASH_MULTIQC="TRUE"
export DO_QIIME_SPLIT="TRUE"
export DO_QIIME_PICK="TRUE"
export DO_PYTHON_OTU_FILTERING="TRUE"
export DO_BIOM_CONVERSION="TRUE"
export DO_CORE_DIV="TRUE"
export DO_PHYLOSEQ="FALSE"




