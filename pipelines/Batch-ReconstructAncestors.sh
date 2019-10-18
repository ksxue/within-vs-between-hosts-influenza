#!/bin/bash

# Script for batch submission of ReconstructAncestors jobs.

# Input parameters.
pipeline="pipelines/ReconstructAncestors.sh"
runs="$1" # File containing the list of alignments, trees, and output directories.
clean="$2" # If 0, then the script will run in its entirety, overwriting existing output.

# Extract the appropriate line of the list of run IDs to submit to the script.
alignment=`cut -f1 -d' ' ${runs} | awk "NR==$SGE_TASK_ID"`
tree=`cut -f2 -d' ' ${runs} | awk "NR==$SGE_TASK_ID"`
dir=`cut -f3 -d' ' ${runs} | awk "NR==$SGE_TASK_ID"`
${pipeline} ${alignment} ${tree} ${dir} ${clean}
