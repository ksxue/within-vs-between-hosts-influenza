#!/bin/bash

# For data downloaded from the SRA, there are often multiple FASTQ
# files (generated from multiple sequencing runs) that correspond
# to a single biological sample.
#
# Given an SRA run table for a particular BioProject and a directory
# containing the FASTQ files downloaded from this BioProject,
# this script
# 1. renames FASTQ files from an SRR designation (column 10) to the listed 
#    Sample_Name_s from the SRA run table, in column 12 in cases in which
#    a sample has only a single SRR match (i.e. it is a singlet)
# 2. identifies samples with exactly two associated runs (SRR numbers),
#    merges the FASTQ files separately for read 1 and read 2,
#    and renames the merged file to match the listed sample name
# 3. outputs a list of samples with more than three associated runs
#
# This script does NOT handle samples with three or more associated runs
# in a comprehensive manner. It also does not download the FASTQ files
# associated with a particular BioProject. It assumes that all SRR files
# listed in the SRA run table have already been downloaded and are present
# in the same directory.
#
# Note that this script renames all singlets but copies doublets,
# meaning that there will remain SRR files to be deleted.
# Additionally, there may be exceptions (samples with >2 SRRs)
# that should be handled individually.
# As such, this script does not permanently rename or delete files
# with the exceptions of singlets, for which renaming is significantly
# faster than copying.
#
# This script is a variant of the typical ConsolidateRenameSRA script.
# It is modified to run on the SRA table that describes the data from
# McCrone et al. 2018.

# Input parameters.
samplesheet="$1" # Give the path to the SRA run table as downloaded from NCBI.
dir="$2" # Give the directory for the downloaded FASTQ files.

# Count the number of runs that correspond to each sample and
# store this information.
cut -f12 ${samplesheet} | tail -n +2 \
  | sort | uniq -c | sed -r 's/^( *[^ ]+) +/\1\t/' | \
  sort -k 1 > ${dir}/counts.txt

# Extract the list of samples that have only a single corresponding SRR run
# and rename the downloaded FASTQ files to match the sample name.
while read sample
do
  awk -v sample="$sample" -F'\t' '$12==sample' ${samplesheet} \
    | cut -f9,12
done < <( awk '$1==1' ${dir}/counts.txt | cut -f2 ) \
  > ${dir}/singlets.txt
  
while read run sample
do
  if [ -f ${dir}/${run}_1.fastq.gz ] && [ -f ${dir}/${run}_2.fastq.gz ];
  then
    echo ${sample}
    mv ${dir}/${run}_1.fastq.gz ${dir}/${sample}_1.fastq.gz
	mv ${dir}/${run}_2.fastq.gz ${dir}/${sample}_2.fastq.gz
  fi
done < ${dir}/singlets.txt

# Extract the names of samples with two SRRs,
# then determine the list of SRRs that correspond to those samples.
# The resulting file should consist of consecutive pairs of lines of the form
# SRR#### sample####
while read sample
do
  awk -v sample="$sample" -F'\t' '$12==sample' ${samplesheet} \
    | cut -f9,12
done < <( awk '$1==2' ${dir}/counts.txt | cut -f2 ) \
  > ${dir}/doublets.txt

# Split the doublets file into two files, each listing a single sample per line.
# Merge those files to take the form
# sample#### SRR#### SRR####
# Read the file of doublets and concatenate separate fastq files
# into a single file, named by the sample rather than the run.
while read sample run1 run2
do
  echo "$sample $run1 $run2"
  if [ ! -f ${dir}/${sample}_1.fastq.gz ] && [ ! -f ${dir}/${sample}_2.fastq.gz ];
  then
    cat ${dir}/${run1}_1.fastq.gz ${dir}/${run2}_1.fastq.gz > ${dir}/${sample}_1.fastq.gz
    cat ${dir}/${run1}_2.fastq.gz ${dir}/${run2}_2.fastq.gz > ${dir}/${sample}_2.fastq.gz
  fi
done < <(join <( awk 'NR%2==1{print $2,"\t",$1}' ${dir}/doublets.txt ) \
  <( awk 'NR%2==0{print $2,"\t",$1}' ${dir}/doublets.txt ))


# Extract the list of samples that have three or more corresponding SRR runs.
while read sample
do
  awk -v sample="$sample" -F'\t' '$12==sample' ${samplesheet} \
    | cut -f9,12
done < <( awk '$1>2' ${dir}/counts.txt | cut -f2 ) \
  > ${dir}/exceptions.txt
  
# Clean up file intermediates.
# rm -f ${dir}/counts.txt
# rm -f ${dir}/singlets.txt
# rm -f ${dir}/doublets.txt