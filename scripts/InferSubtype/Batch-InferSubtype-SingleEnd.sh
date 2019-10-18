#!/bin/bash

# This script is meant to be run from the top level of the Github repository.
# Script for batch submission of InferSubtype jobs.

# Input parameters.
pipeline="scripts/InferSubtype/InferSubtype-SingleEnd.sh"
samplesheet="$1" # Give the path to a samplesheet listing rawread1, rawread2, trimmed1, trimmed2, samplename.
dir="$2" # Give the desired working directory for the trimmed FASTQ files.
references="$3" # Give the list of reference genomes to be mapped against.
numreads=$4 # Give the number of reads to use for test mapping.

# Extract the appropriate line of the list of run IDs to submit to the script.
fastq1=`cut -f1 -d' ' ${samplesheet} | awk "NR==$SGE_TASK_ID"`
fastq2=`cut -f2 -d' ' ${samplesheet} | awk "NR==$SGE_TASK_ID"`
sample=`cut -f5 -d' ' ${samplesheet} | awk "NR==$SGE_TASK_ID"`
${pipeline} ${fastq1} ${fastq2} ${sample} ${dir} ${references} ${numreads}