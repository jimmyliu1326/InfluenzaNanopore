#!/usr/bin/env bash

usage() {
echo "
influenza_consensus.sh [options]

Required arguments:
-i|--input          .csv file with first column as sample name and second column as path to fastq, no headers required
-o|--output         Path to output directory
--db                Path to Centrifuge database for raw read taxonomic classification

Optional arguments:
-t|--threads        Number of threads [Default = 32]
-h|--help           Display help message
"
}

# get script directory
script_dir=$(dirname $(realpath $0))

# parse arguments
opts=`getopt -o hi:o:t: -l help,input:,output:,threads:,db: -- "$@"`
eval set -- "$opts"
if [ $? != 0 ] ; then echo "influenza_consensus: Invalid arguments used, exiting"; usage; exit 1 ; fi
if [[ $1 =~ ^--$ ]] ; then echo "influenza_consensus: Invalid arguments used, exiting"; usage; exit 1 ; fi

while true; do
    case "$1" in
        -i|--input) INPUT_PATH=$2; shift 2;;
        -o|--output) OUTPUT_PATH=$2; shift 2;;
        --db) DB_PATH=$2; shift 2;;
        -t|--threads) THREADS=$2; shift 2;;
        --) shift; break ;;
        -h|--help) usage; exit 0;;
    esac
done

if test -z $INPUT_PATH; then echo "influenza_consensus: Required argument -i is missing, exiting"; exit 1; fi
if test -z $OUTPUT_PATH; then echo "influenza_consensus: Required argument -o is missing, exiting"; exit 1; fi
if test -z $DB_PATH; then echo "influenza_consensus: Required argument --db is missing, exiting"; exit 1; fi

# check dependencies
medaka_consensus -h 2&>1 /dev/null
if [[ $? != 0 ]]; then echo "influenza_consensus: medaka cannot be called, check its installation"; exit 1; fi

seqtk 2&>1 /dev/null
if [[ $? != 1 ]]; then echo "influenza_consensus: samtools cannot be called, check its installation"; exit 1; fi

snakemake -h > /dev/null
if [[ $? != 0 ]]; then echo "influenza_consensus: snakemake cannot be called, check its installation"; exit 1; fi

racon -h > /dev/null
if [[ $? != 0 ]]; then echo "influenza_consensus: racon cannot be called, check its installation"; exit 1; fi

centrifuge -h > /dev/null
if [[ $? != 0 ]]; then echo "influenza_consensus: centrifuge cannot be called, check its installation"; exit 1; fi

# validate input files
if ! test -f $INPUT_PATH; then echo "influenza_consensus: Input sample file does not exist, exiting"; exit 1; fi

while read lines; do
  sample=$(echo $lines | cut -f1 -d',')
  path=$(echo $lines | cut -f2 -d',')
  if ! test -f $path; then
    echo "influenza_consensus: Fastq file for ${sample} cannot be found, check its path listed in the input file, exiting"
    exit 1
  fi
done < $INPUT_PATH

# create output directory if does not exist
if ! test -d $OUTPUT_PATH; then mkdir -p $OUTPUT_PATH; fi

# Set default number of threads if not specified
if test -z $THREADS; then THREADS=32; fi

# call snakemake
snakemake --snakefile $script_dir/SnakeFile --cores $THREADS \
  --config samples=$(realpath $INPUT_PATH) outdir=$OUTPUT_PATH pipeline_dir=$script_dir centrifuge_db=$DB_PATH