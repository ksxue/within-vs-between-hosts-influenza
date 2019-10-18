#!/bin/bash

# This script is meant to be run from the top level of the Github repository.
# Script for batch submission of SRADownloadFASTQ jobs.

# Input parameters.
pipeline="pipelines/TrimAdaptersNextera.sh"
samplesheet="$1" # Give the path to a samplesheet listing rawread1, rawread2, trimmed1, trimmed2, samplename.
dir="$2" # Give the desired working directory for the trimmed FASTQ files.
clean="$3" # If 0, then the script will run in its entirety, overwriting existing output.

# Extract the appropriate line of the list of run IDs to submit to the script.
fastq1=`cut -f1 -d' ' ${samplesheet} | awk "NR==$SGE_TASK_ID"`
fastq2=`cut -f2 -d' ' ${samplesheet} | awk "NR==$SGE_TASK_ID"`
sample=`cut -f5 -d' ' ${samplesheet} | awk "NR==$SGE_TASK_ID"`
${pipeline} ${fastq1} ${fastq2} ${sample} ${dir} ${clean}