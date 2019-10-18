#!/bin/bash

# This script takes raw sequencing reads,
# trims adapters, and outputs trimmed reads in the specified directory.

# Load dependencies.
module load modules modules-init modules-gs
module load ncbi-sra/2.8.2-1

# Input parameters.
dir="$1" # Give the desired working directory for the downloaded FASTQ files.
run="$2" # Give the SRR###### number for the run to download.
clean="$3" # If 0, then the script will run in its entirety, overwriting existing output.

echo "directory: ${dir}"
echo "run: ${run}"
echo "clean: ${clean}"

# Check that directory exists.
mkdir -p ${dir}

# Download reads.
# Split reads into two files representing paired-end reads.
# Also zip the resulting file.
if [ ! -f ${dir}/${run}_1.fastq.gz ] || [ ! -f ${dir}/${run}_2.fastq.gz ] || \
  [ ${clean} -eq "0" ];
then
  fastq-dump -I --split-files --skip-technical --gzip -O ${dir} ${run}
fi

