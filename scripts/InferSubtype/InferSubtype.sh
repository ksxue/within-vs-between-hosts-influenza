#!/bin/bash

# Script is meant to be run from the top level of the Github repository.
# Given FASTQ files for a given sample and list of reference genomes,
# map the first n reads to each reference genome,
# summarize the mapping rates for each genome, and determine the best match.

# Input parameters.
fastq1="$1" # Give read 1 file.
fastq2="$2" # Give read 2 file.
samplename="$3" # Give the sample name.
dir="$4" # Give the desired location of intermediate output files.
references="$5" # Give list of references in path/name format.
numreads=$6 # Give the number of reads to be mapped.

# For each reference genome on the list,
# map the first n reads against each reference genome.
while read reference refname
do
  # Concatenate the sample and reference names
  # so that all files reflect both the strain and reference.
  sampleshort="$samplename"
  sample="$samplename-$refname"
  numreads=$[numreads*4] # Multiply number of reads by 4 to get number of lines.

  # Map the first n reads against the specified reference genome.
  bowtie2 --very-sensitive-local --un-conc-gz ${dir}/${sample}-unmapped \
	-X 2300 \
	-x ${reference} \
	-1 <( zcat ${fastq1} | head -n ${numreads} ) \
	-2 <( zcat ${fastq2} | head -n ${numreads} ) \
	-S ${dir}/${sample}.sam \
	2> ${dir}/${sample}.bt2.log
	
  # Parse .bt2.log file to obtain the mapping rate.
  # Output a file with a single line in format
  # FASTQ1 FASTQ2 sample reference refname maprate
  maprate="$(tail -n 1 ${dir}/${sample}.bt2.log | cut -f 1 -d' ')"
  maprate=${maprate%%\%}
  echo ${fastq1} ${fastq2} ${samplename} ${reference} ${refname} ${maprate} \
    > ${dir}/${sample}.temp
  
  # Remove other intermediate files.
  rm -f ${dir}/${sample}.bt2.log
  rm -f ${dir}/${sample}.sam
  rm -f ${dir}/${sample}-unmapped.1
  rm -f ${dir}/${sample}-unmapped.2
done < ${references}

# Summarize mapping information for all of the reference genomes.
# Sort in descending order of mapping rate.
cat ${dir}/${samplename}-*.temp | sort -r -n -k 6 > ${dir}/${samplename}.map

# Remove intermediate files.
rm -f ${dir}/${samplename}-*.temp
