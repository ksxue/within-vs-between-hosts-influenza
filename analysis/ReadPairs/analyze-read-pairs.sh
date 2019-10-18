#!/bin/bash

# Script is meant to be run from the top level of the Github repository.
# Script analyzes read ID and run ID information from reconstructed read pairs.

# Input parameters and directory paths.
infile=$1 # Takes a FASTQ file containing both read pairs on a single, tab-delimited line.
outdir=$2 # Output directory.

# Parse out the sample IDs of each read 1 and read 2 sequence.
# Identify the read 1 IDs.
zcat ${infile} \
  | awk '{ print $1 }' | awk -F'/' '{ print $2 }' \
  > ${outdir}/read1IDs_all.txt
  
# Identify the read 2 IDs.
zcat ${infile} \
  | awk '{ print $5 }' | awk -F'/' '{ print $2 }' \
  > ${outdir}/read2IDs_all.txt


# Parse out the run and lane IDs from the beginning of the read identifier.
zcat ${infile} \
  | awk '{ print $1 }' | awk -F':' '{ OFS=":"; print $1,$2 }' \
  > ${outdir}/read1runIDs_all.txt
  
zcat ${infile} \
  | awk '{ print $5 }' | awk -F':' '{ OFS=":"; print $1,$2 }' \
  > ${outdir}/read2runIDs_all.txt

# Merge the files corresponding to read 1 and read 2 origins,
# sort by the identity of each read pair,
# and count the number of occurrences of each combination.
paste ${outdir}/read1runIDs_all.txt ${outdir}/read1IDs_all.txt \
  ${outdir}/read2runIDs_all.txt ${outdir}/read2IDs_all.txt |
  sort | uniq -c > ${outdir}/readIDrun_counts.txt

# Generate full read 1 and read 2, regularly formatted FASTQ files.
zcat ${infile} \
  | awk '{ OFS="\n"; print $1,$2,$3,$4 }' | gzip -f > ${outdir}/read1_all.fastq.gz
zcat ${infile} \
  | awk '{ OFS="\n"; print $5,$6,$7,$8 }' | gzip -f > ${outdir}/read2_all.fastq.gz
