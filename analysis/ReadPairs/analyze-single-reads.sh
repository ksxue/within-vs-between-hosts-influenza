#!/bin/bash

# Script is meant to be run from the top level of the Github repository.
# Script analyzes read ID information from remaining single-end reads.

# Input parameters and directory paths.
infile=$1 # Takes a FASTQ file containing single-end reads, each on a single, tab-delimited line.
outdir=$2 # Output directory.

# Parse out the sample IDs of each unpaired read.
zcat ${infile} \
  | awk '{ print $1 }' | awk -F'/' '{ print $2 }' \
  > ${outdir}/single_IDs_all.txt

# Sort by the identity of each unpaired read,
# and count the number of occurrences of each sample.
sort ${outdir}/single_IDs_all.txt | uniq -c > ${outdir}/single_counts.txt

