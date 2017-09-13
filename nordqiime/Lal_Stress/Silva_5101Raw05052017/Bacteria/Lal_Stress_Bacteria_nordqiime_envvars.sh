export IN2CSV_PATH="/usr/local/bin/in2csv"
export FASTQC_PATH="/home/mark/gitstage/uderica/nordqiime/FastQC/fastqc"
export MULTIQC_PATH="/usr/local/bin/multiqc"

export MULTIQC_CONFIG="/home/mark/gitstage/uderica/nordqiime/multiqc_config.yaml"

export THREAD_COUNT="7"

export rootFold="/home/mark/gitstage/uderica/nordqiime"


# this is WAY easier and more robust if there are no space in the file and directory names

export projFold="/Lal_Stress/Silva_5101Raw05052017"

# export projName=`echo $projFold | sed -r 's/\s+/_/g' | sed 's/^\///'`

export projName="Lal_Stress"

# either bacteria or fungi
# would be nice to infer from the domain folder name
# will effect the otu picking and cordeiv settings
export lifeDomain="bacteria"
export domainFold="/Bacteria"
export fastqFold="/selected_raw"
export filePathFile="file path.xlsx"

export sampleIdCol="5"
export R1Col="6"
export R2Col="7"

export fastqFileSuffix=".fastq"

export FLASH_PATH="/home/mark/gitstage/uderica/nordqiime/FLASH-1.2.11/flash"
export MIN_FLASH_OVERLAP=50
export MAX_FLASH_OVERLAP=250

# export QIIME_MAPPING_FILE="${PROJ_GENE_HOME}/map.txt"
# haven't expliclitly checked the mapping file yet

export REF_SEQ_FILE="/home/mark/gitstage/uderica/nordqiime/gg_13_8_otus/rep_set/97_otus.fasta"

export BACT_ALIAS="${rootFold}${projFold}${domainFold}${fastqFold}"

# placeholder
# see main nordqiime script for notes on methodically calculating this
# export CORE_DIV_MIN_DEPTH=10000
# Min: 27347.0
# Max: 99927.0
# Median: 54502.000
# Mean: 55241.407
# Std. dev.: 16687.858
# mean - 2*sd
export CORE_DIV_MIN_DEPTH=18000

export CORE_DIV_EXP_FACTOR="treatment"

