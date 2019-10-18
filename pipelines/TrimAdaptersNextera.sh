#!/bin/bash

# This script takes raw sequencing reads,
# trims adapters, and outputs trimmed reads in the specified directory.

# Load dependencies.
module load modules modules-init modules-gs
module load python/3.6.4

# Input parameters.
fastq1="$1"
fastq2="$2"
sample="$3"
dir="$4"
clean="$5"

# Trim adapter sequences and bases below a quality threshold of 25.
# Also remove all reads that are shorter than 20 bases after trimming.
echo "Trim adapter sequences."
if [ ! -f ${dir}/${sample}_trimmed-R1.fastq.gz ] || [ ! -f ${dir}/${sample}_trimmed-R2.fastq.gz ] || \
  [ ${clean} -eq "0" ];
then
  cutadapt -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG  -A GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
    -q 25 -m 20 \
	-o ${dir}/${sample}_trimmed-R1.fastq.gz -p ${dir}/${sample}_trimmed-R2.fastq.gz \
	${fastq1} \
	${fastq2} \
	> ${dir}/${sample}.cutadapt.log
fi