#!/bin/bash

# Script for batch submission of ExtractMutationFrequencies jobs.

# Input parameters.
pipeline="pipelines/ExtractMutationFrequencies.sh"
runs="$1" # File containing the list of alignments, trees, and output directories.
clean="$2" # If 0, then the script will run in its entirety, overwriting existing output.

# Extract the appropriate line of the list of run IDs to submit to the script.
tree=`cut -f1 -d' ' ${runs} | awk "NR==$SGE_TASK_ID"`
seq=`cut -f2 -d' ' ${runs} | awk "NR==$SGE_TASK_ID"`
outgroup=`cut -f3 -d' ' ${runs} | awk "NR==$SGE_TASK_ID"`
out=`cut -f4 -d' ' ${runs} | awk "NR==$SGE_TASK_ID"`
${pipeline} ${tree} ${seq} ${outgroup} ${out} ${clean}
