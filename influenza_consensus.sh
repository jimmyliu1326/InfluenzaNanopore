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
-s|--segment        Target specific Influenza A genomic segments for consensus calling with each segment number delimited by a comma (Example: -s 1,2,5,6)
--subsample         Specify the target coverage for consensus calling [Default = 1000]
-m|--model          Specify the flowcell chemistry used for Nanopore sequencing {Options: r9, r10} [Default = r9]
--notrim            Disable adaptor trimming by Porechop
--keep-tmp          Keep all temporary files
-h|--help           Display help message
"
}

# define global variables
script_dir=$(dirname $(realpath $0))
script_name=$(basename $0 .sh)
DB_PATH="/db/viralRefSeq_InfA_custom"

# parse arguments
opts=`getopt -o hi:o:t:s:m: -l help,input:,output:,threads:,db:,notrim,segment:,model:,keep-tmp,subsample:,mode: -- "$@"`
eval set -- "$opts"
if [ $? != 0 ] ; then echo "influenza_consensus: Invalid arguments used, exiting"; usage; exit 1 ; fi
if [[ $1 =~ ^--$ ]] ; then echo "influenza_consensus: Invalid arguments used, exiting"; usage; exit 1 ; fi

while true; do
    case "$1" in
        -i|--input) INPUT_PATH=$2; shift 2;;
        -o|--output) OUTPUT_PATH=$2; shift 2;;
        --db) DB_PATH=$2; shift 2;;
        -t|--threads) THREADS=$2; shift 2;;
        --mode) MODE=$2; shift 2;;
        -m|--model) MODEL=$2; shift 2;;
        --subsample) SUBSAMPLE=$2; shift 2;;
        --notrim) TRIM=0; shift 1;;
        --keep-tmp) KEEP_TMP=1; shift 1;;
        -s|--segment) SEGMENTS=$2; shift 2;;
        --) shift; break ;;
        -h|--help) usage; exit 0;;
    esac
done

# check if required arguments are given
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

centrifuge -h > /dev/null
if [[ $? != 0 ]]; then echo "influenza_consensus: centrifuge cannot be called, check its installation"; exit 1; fi

porechop -h > /dev/null
if [[ $? != 0 ]]; then echo "influenza_consensus: porechop cannot be called, check its installation"; exit 1; fi

# validate model parameter input if specified
if ! test -z $MODEL; then
  # test if invalid characters used
  if ! [[ $MODEL =~ ^(r9|r10)$ ]]; then echo "Invalid model specification passed to the -m argument, exiting"; fi
  # set medaka model
  if [[ $MODEL == "r9" ]]; then MODEL="r941_min_high_g360"; else MODEL="r103_min_high_g360"; fi
else
  # Set default model if not specified
  if test -z $MODEL; then MODEL="r941_min_high_g360"; fi
fi

# validate segment parameter input if specified
if ! test -z $SEGMENTS; then
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
else
  # Set default segments to produce = [1:8] if not specified
  if test -z $SEGMENTS; then SEGMENTS="1,2,3,4,5,6,7,8"; fi
fi

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

# Set default keep temporary files (KEEP_TMP) to 0 if not specified
if test -z $KEEP_TMP; then KEEP_TMP=0; fi

# Set default subsample depth to 1000 if not specified
if test -z $SUBSAMPLE; then SUBSAMPLE=1000; fi

# Set default read selection mode to dynamic if not specified
if test -z $MODE; then MODE="dynamic"; fi

# validate read selection mode if specified
if ! [[ $MODE =~ ^(dynamic|static)$ ]]; then echo "${script_name}: Invalid mode option passed to the --mode argument, accepted values are dynamic/static, exiting"; exit 1; fi

# Remove existing analysis html report and summary statistics file in OUTPUT_PATH
if test -f $OUTPUT_PATH/InfA_analysis_viz.html; then rm $OUTPUT_PATH/InfA_analysis_viz.html; fi
if test -f $OUTPUT_PATH/summary_statistics.csv; then rm $OUTPUT_PATH/summary_statistics.csv; fi

while read lines; do
  sample=$(echo $lines | cut -f1 -d',')
  if test -d $OUTPUT_PATH/$sample; then
    rm -rf $OUTPUT_PATH/$sample
  fi
done < $INPUT_PATH

# call snakemake
snakemake --snakefile $script_dir/SnakeFile --cores $THREADS \
  --keep-going \
  --config samples=$(realpath $INPUT_PATH) \
  outdir=$(realpath $OUTPUT_PATH) \
  segments="$SEGMENTS" \
  pipeline_dir=$script_dir \
  centrifuge_db=$DB_PATH \
  trim=$TRIM \
  model=$MODEL \
  threads=$THREADS \
  subsample=$SUBSAMPLE \
  mode=$MODE

# clean up temporary directories
if [[ $KEEP_TMP -eq 0 ]]; then
  while read lines; do
    sample=$(echo $lines | cut -f1 -d',')
    for dir in nanoplot fastq medaka centrifuge binned_fastq porechop draft_consensus subsample_fastq; do
      if test -d $OUTPUT_PATH/$sample/$dir; then
        rm -rf $OUTPUT_PATH/$sample/$dir
      fi
    done
  done < $INPUT_PATH
fi

# check status of consensus building for each segment
IFS=","
read -r -a segment_array <<< $SEGMENTS
while read lines; do
  sample=$(echo "$lines" | cut -f1 -d',')
  echo -e "influenza_consensus: \033[0;33mConsensus Building Summary for $sample\033[0m"
  for i in ${segment_array[@]}; do
    if test -f $OUTPUT_PATH/$sample/consensus/segment_${i}.fa; then
    echo -e " \033[0;33mSegment ${i}\033[0m: \033[0;32mSUCCESS\033[0m"
    else 
    echo -e " \033[0;33mSegment ${i}\033[0m: \033[0;31mFAIL\033[0m"
    fi
  done
done < $INPUT_PATH

