#!/bin/bash

# Script is meant to be run from the top level of the Github repository.
# Script organizes raw sequencing files for Dinis et al. 2016.


# Dinis et al. 2016, J Virol.
# Raw sequencing data was obtained through personal communication
# with the authors.
# Sequencing reads were provided in the form of a single FASTQ file
# per sample, presumably the result of exporting reads from an alignment file.

# We sorted the raw reads into read pairs based on read ID information.

# First, organize the appropriate directories.
project="flu-Dinis"
analysis="DownloadDataCallVariants"
dir="nobackup/${project}/raw"
samplesheet="data/metadata/${project}/samplesheet.txt"
mkdir -p nobackup/${project}/raw
mkdir -p nobackup/${project}/trimmed
mkdir -p data/metadata/${project}

# Create a file containing the locations of all of the raw reads.
rawdir="../160527-public-datasets/raw/2016-Dinis-JVirol"
ls -d ${rawdir}/H1N1/* > data/metadata/${project}/raw.txt
ls -d ${rawdir}/H3N2/* >> data/metadata/${project}/raw.txt

# Create a samplesheet containing the paths to
# raw and trimmed reads.
raw="nobackup/${project}/raw"
trimmed="nobackup/${project}/trimmed"
while read fastq
do
  sample=${fastq%%_201*}
  sample=${sample##*/}
  sample=${sample##*_}
  echo ${raw}/${sample}_1.fastq.gz ${raw}/${sample}_2.fastq.gz \
    ${trimmed}/${sample}_trimmed-R1.fastq.gz \
	${trimmed}/${sample}_trimmed-R2.fastq.gz \
	${sample} ${project}
done < data/metadata/${project}/raw.txt \
  > data/metadata/${project}/samplesheet.txt
  
# Split read pairs to form separate FASTQ files
# for read 1 and read 2.
numsamples="$(wc -l ${samplesheet} | cut -f1 -d' ')"

if [ ! "$(ls ${dir}/*.fastq.gz | wc -l)" -eq $(( numsamples * 2 )) ] && \
  [ -f data/metadata/${project}/samplesheet.txt ]; then
  echo "Reconstruct read pairs."
  # Submit batch jobs to reconstruct read pairs for all flu samples.
  qsub -cwd -N SplitReads \
    -t 1-${numsamples} -tc 250 \
    -o nobackup/sge/ -e nobackup/sge/ \
    analysis/${analysis}/batch-split-read-pairs-Dinis.sh \
    data/metadata/${project}/raw.txt data/metadata/${project}/samplesheet.txt \
    ${dir}  
fi


  