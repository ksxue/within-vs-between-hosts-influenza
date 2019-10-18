#!/bin/bash

# Script is meant to be run from the top level of the Github repository.
# Script extracts consecutive read pairs from a file of sorted FASTQ reads.

# Input parameters and directory paths.
project="flu-Poon"
analysis="ReadPairs"
rawdir="nobackup/flu-Poon/raw/rawfiles"
outdir="nobackup/${analysis}/"
mkdir -p ${outdir}


# For each FASTQ file in the dataset,
# merge the four lines of each FASTQ read into a single, tab-delimited line,
# sort the reads based on the header lines, 
# append the file name to the header for each read
# and concatenate the FASTQ samples for all files
while read fastq strain ID
do
  sample="${strain}-${ID}"
  echo ${rawdir}/${fastq}
  zless ${rawdir}/${fastq} | paste - - - - -d'\t' | sort -k1,1 | \
    sed "s/\t/-${sample}\t/" | gzip > ${outdir}/${sample}.fastq.gz
done < data/metadata/${project}/Synapse-samplenames.txt

# Merge all sorted files into a single sorted file
# for the entire study.
# Do this for only the first 100 reads of each file to start.
cmd="sort -m -k1,1"
while read fastq strain ID
do
  sample="${strain}-${ID}"
  cmd="$cmd <(gunzip -c '${outdir}/${sample}.fastq.gz' | head -n 100)"
done < data/metadata/${project}/Synapse-samplenames.txt
eval "$cmd" > ${outdir}/sorted_all_test.fastq
gzip ${outdir}/sorted_all_test.fastq

# Merge all sorted files into a single sorted file
# for the entire study.
cmd="sort -m -k1,1"
while read fastq strain ID
do
  sample="${strain}-${ID}"
  cmd="$cmd <( gunzip -c '${outdir}/${sample}.fastq.gz' )"
done < data/metadata/${project}/Synapse-samplenames.txt
eval "$cmd" > ${outdir}/sorted_all.fastq
gzip ${outdir}/sorted_all.fastq