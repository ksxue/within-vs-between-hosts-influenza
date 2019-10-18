#!/bin/bash

# Script for batch submission of SRADownloadFASTQ jobs.

# Input parameters.
pipeline="pipelines/SRADownloadFASTQ.sh"
dir="$1" # Give the desired working directory for the downloaded FASTQ files.
runs="$2" # File containing ONLY the list of sequencing run numbers to be submitted.
clean="$3" # If 0, then the script will run in its entirety, overwriting existing output.

# Extract the appropriate line of the list of run IDs to submit to the script.
run=`awk "NR==$SGE_TASK_ID" ${runs}`
${pipeline} ${dir} ${run} ${clean}