#!/bin/bash

# Script is meant to be run from the top level of the Github repository.
# Script extracts unpaired reads from a file of sorted FASTQ reads.

# Input parameters and directory paths.
infile=$1 # Takes a header-sorted, gzipped FASTQ file with each read on a single, tab-delimited line.
outfile=$2 # Outputs a FASTQ file containing single-end reads with each read on a single, tab-delimited line.

# Iterate through consecutive groups of three reads in a sorted read file.
# Determine whether the read ID of the second read in the group
# matches the read ID of the first or third read. 
# Do not make use of pair information (i.e. 1 or 2).
# If the second read does not match the ID of the first or third read,
# then consider it a single-end read and output it to a separate file.
prevRead=""
prevReadID=""
prevReadRest=""
currRead=""
currReadID=""
currReadRest=""
while read nextRead nextReadRest
do
  nextReadID=${nextRead%%/*}
  if [ "${currReadID}" != "${prevReadID}" ] && \
     [ "${currReadID}" != "${nextReadID}" ]; then
    echo -e ${currRead}'\t'${currReadRest}
  fi
  prevRead=${currRead}
  prevReadID=${currReadID}
  prevReadRest=${currReadRest}
  currRead=${nextRead}
  currReadID=${nextReadID}
  currReadRest=${nextReadRest}
done < <( zcat ${infile} ) \
  > ${outfile}
# Determine whether the last line matches the ID of the previous line.
if [ "${currReadID}" != "${prevReadID}" ] ; then
  echo -e ${currRead}'\t'${currReadRest}
fi >> ${outfile}
gzip -f ${outfile}
