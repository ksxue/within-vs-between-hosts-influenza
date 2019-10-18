#!/bin/bash

# This script is meant to be run from the top level of the Github repository.
# Script for batch submission of SplitReadPairs.sh jobs.

# Input parameters.
pipeline="analysis/DownloadDataCallVariants/split-read-pairs-Dinis.sh"
raw="$1" # Give the path to a file listing the locations of the FASTQ files.
samplesheet="$2" # Give the path to a samplesheet listing read1, read2, trimmed1, trimmed2, sample.
dir="$3" # Give the desired working directory for the trimmed FASTQ files.

# Extract the appropriate line of the list of run IDs to submit to the script.
fastq=`awk "NR==$SGE_TASK_ID" ${raw}`
outstem=`cut -f5 -d' ' ${samplesheet} | awk "NR==$SGE_TASK_ID"`
outstem=${dir}/${outstem}
${pipeline} ${fastq} ${outstem}