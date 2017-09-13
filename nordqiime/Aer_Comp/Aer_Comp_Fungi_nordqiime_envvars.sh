export IN2CSV_PATH="/usr/local/bin/in2csv"
export FASTQC_PATH="/home/mark/gitstage/uderica/nordqiime/FastQC/fastqc"
export MULTIQC_PATH="/usr/local/bin/multiqc"

export MULTIQC_CONFIG="/home/mark/gitstage/uderica/nordqiime/multiqc_config.yaml"

export THREAD_COUNT="7"

export rootFold="/home/mark/gitstage/uderica/nordqiime"


# this is WAY easier and more robust if there are no space in the file and directory names
export projFold="/Aer_Comp"

export projName=`echo $projFold | sed -r 's/\s+/_/g' | sed 's/^\///'`

# either bacteria or fungi
# would be nice to infer from the domain folder name
# will effect the otu picking and cordeiv settings
export lifeDomain="fungi"
export domainFold="/Aer_Comp_Fungi"
export fastqFold="/raw_data"
export filePathFile="Aer Comp file path.xlsx"

export sampleIdCol="5"
export R1Col="6"
export R2Col="7"

export fastqFileSuffix=".fastq"

export FLASH_PATH="/home/mark/gitstage/uderica/nordqiime/FLASH-1.2.11/flash"
export MIN_FLASH_OVERLAP=50
export MAX_FLASH_OVERLAP=250

# export QIIME_MAPPING_FILE="${PROJ_GENE_HOME}/map.txt"
# haven't expliclitly checked the mapping file yet

export REF_SEQ_FILE="/home/mark/gitstage/uderica/nordqiime/its_12_11_otus/rep_set/97_otus.fasta"

# should globally renames this sto something like DOMAIN_ALIAS
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
export CORE_DIV_MIN_DEPTH=20000

export CORE_DIV_EXP_FACTOR="Treatment"

# export OTU_PICK_PARAM_FILE="${rootFold}/otu_pick_param.txt"

