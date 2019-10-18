#!/bin/bash

# Script for batch submission of BuildTreeRAxML jobs.

# Input parameters.
pipeline="pipelines/BuildTreeRAxML.sh"
dir="$1" # Give the desired output directory for the trees that are constructed.
runs="$2" # File containing the list of sequences and their tree names.
clean="$3" # If 0, then the script will run in its entirety, overwriting existing output.

# Extract the appropriate line of the list of run IDs to submit to the script.
sequences=`cut -f1 -d' ' ${runs} | awk "NR==$SGE_TASK_ID"`
name=`cut -f2 -d' ' ${runs} | awk "NR==$SGE_TASK_ID"`
${pipeline} ${sequences} ${name} ${dir} ${clean}