#!/bin/bash

# Script is meant to be run from the top level of the Github repository.
# Script extracts consecutive read pairs from a file of sorted FASTQ reads.

# Input parameters and directory paths.
infile=$1 # Takes a header-sorted, gzipped FASTQ file with each read on a single, tab-delimited line.
outfile=$2 # Outputs a FASTQ file containing both read pairs on a single, tab-delimited line.

# Iterate through all pairs of consecutive reads in a sorted read file.
# Extract the read identifier, and identify pairs of reads with the same read ID
# but different pair numbers (i.e. 1 or 2).
prevRead=""
prevReadID=""
prevReadNum=""
prevReadRest=""
while read currRead currReadRest
do
  currReadID=${currRead%%/*}
  currReadNum=${currRead##*/}
  currReadNum=${currReadNum%%-*}
  if [ "${currReadID}" == "${prevReadID}" ] && \
     [ "${currReadNum}" != "${prevReadNum}" ]; then
	echo -e ${prevRead}'\t'${prevReadRest}'\t'${currRead}'\t'${currReadRest}
  fi
  prevRead=${currRead}
  prevReadID=${currReadID}
  prevReadNum=${currReadNum}
  prevReadRest=${currReadRest}
done < <( zcat ${infile} ) \
  > ${outfile}
gzip -f ${outfile}