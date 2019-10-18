#!/bin/bash

# This script is meant to be run from the top level of the Github repository.
# Script for batch submission of InferSubtype jobs.

# Input parameters.
pipeline="scripts/CallVariants/CallVariants-HardFilter-Longitudinal.R"
samplesheet="$1" # Give the path to a samplesheet listing rawread1, rawread2, trimmed1, trimmed2, samplename.
indir="$2" # Give the working directory for the input annotated variant files.
outdir="$3" # Give the desired working directory for the output files.
minfreq=$4 # Give the minimum frequency required to call a variant.
mincoverage=$5 # Give the minimum coverage at a site required to call a variant.

# Extract the appropriate line of the list of run IDs to submit to the script.
sample=`cut -f5 -d' ' ${samplesheet} | awk "NR==$SGE_TASK_ID"`
refname=`cut -f8 -d' ' ${samplesheet} | awk "NR==$SGE_TASK_ID"` # Extract reference name.
Rscript ${pipeline} ${indir}/${sample}-${refname}-realigned-annotated.summary \
  ${outdir}/${sample}-${refname}.variants ${minfreq} ${mincoverage}