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
-s|--segment        Limit consensus sequence calling for specific Influenza A genomic segments with each segment number delimited by a comma (Example: -s 1,2,5,6)
--notrim            Disable adaptor trimming by Porechop
-h|--help           Display help message
"
}

# get script directory
script_dir=$(dirname $(realpath $0))

# parse arguments
opts=`getopt -o hi:o:t:s: -l help,input:,output:,threads:,db:,notrim,segment: -- "$@"`
eval set -- "$opts"
if [ $? != 0 ] ; then echo "influenza_consensus: Invalid arguments used, exiting"; usage; exit 1 ; fi
if [[ $1 =~ ^--$ ]] ; then echo "influenza_consensus: Invalid arguments used, exiting"; usage; exit 1 ; fi

while true; do
    case "$1" in
        -i|--input) INPUT_PATH=$2; shift 2;;
        -o|--output) OUTPUT_PATH=$2; shift 2;;
        --db) DB_PATH=$2; shift 2;;
        -t|--threads) THREADS=$2; shift 2;;
        --notrim) TRIM=0; shift 1;;
        -s|--segment) SEGMENTS=$2; shift 2;;
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
if [[ $? != 1 ]]; then echo "influenza_consensus: seqtk cannot be called, check its installation"; exit 1; fi

snakemake -h > /dev/null
if [[ $? != 0 ]]; then echo "influenza_consensus: snakemake cannot be called, check its installation"; exit 1; fi

racon -h > /dev/null
if [[ $? != 0 ]]; then echo "influenza_consensus: racon cannot be called, check its installation"; exit 1; fi

centrifuge -h > /dev/null
if [[ $? != 0 ]]; then echo "influenza_consensus: centrifuge cannot be called, check its installation"; exit 1; fi

porechop -h > /dev/null
if [[ $? != 0 ]]; then echo "influenza_consensus: porechop cannot be called, check its installation"; exit 1; fi

# validate segment input
for segment in $(echo $SEGMENTS | sed 's/,/\n/g'); do
  # check if individual elements are valid integers between [1:8]
  if [[ $(echo -n $segment | wc -m) -ne 1 ]]; then
    echo "Invalid characters passed to the -s argument, exiting"
    exit 1
  fi
  if ! [[ $segment =~ ^[1-8]$ ]]; then
    echo "Invalid characters passed to the -s argument, exiting"
    exit 1
  fi
done

# validate input samples.csv
if ! test -f $INPUT_PATH; then echo "influenza_consensus: Input sample file does not exist, exiting"; exit 1; fi

while read lines; do
  sample=$(echo $lines | cut -f1 -d',')
  path=$(echo $lines | cut -f2 -d',')
  if ! test -d $path; then
    echo "influenza_consensus: ${sample} directory cannot be found, check its path listed in the input file, exiting"
    exit 1
  fi
done < $INPUT_PATH

# validate database path
if ! test -f ${DB_PATH}.1.cf; then echo "influenza_consensus: Specified Centrifuge database does not exist, exiting"; exit 1; fi

# create output directory if does not exist
if ! test -d $OUTPUT_PATH; then mkdir -p $OUTPUT_PATH; fi

# Set default number of threads if not specified
if test -z $THREADS; then THREADS=32; fi

# Set default trim mode to true if not specified
if test -z $TRIM; then TRIM=1; fi

# Set default segments to produce to [1:8] if not specified
if test -z $SEGMENTS; then SEGMENTS="1,2,3,4,5,6,7,8,"; fi

# call snakemake
snakemake --dag --snakefile $script_dir/SnakeFile --cores $THREADS \
  --config samples=$(realpath $INPUT_PATH) \
  outdir=$OUTPUT_PATH \
  segments="$SEGMENTS" \
  pipeline_dir=$script_dir \
  centrifuge_db=$DB_PATH \
  trim=$TRIM

# clean up temporary directories
while read lines; do
  sample=$(echo $lines | cut -f1 -d',')
  for dir in fastq medaka centrifuge target racon porechop; do
    if test -d $OUTPUT_PATH/$sample/$dir; then
      rm -rf $OUTPUT_PATH/$sample/$dir
    fi
  done
done < $INPUT_PATH