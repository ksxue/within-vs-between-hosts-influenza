#!/bin/bash

# Script is meant to be run from the top level of the Github repository.
# Script organizes raw sequencing files for Debbink et al. 2017.

# Input paths.
BatchSRADownloadFASTQ="pipelines/Batch-SRADownloadFASTQ.sh"
ConsolidateRenameSRA="scripts/ConsolidateRenameSRA.sh"
clean=1 # If 0, analysis will run from beginning, overwriting existing intermediates.

# Xue et al. 2017, eLife
# Download raw FASTQ files from the SRA Bioproject PRJNA364676.
# Use the RunInfo Table for this BioProject to find SRA accessions.
# Rename files to reflect biological sample names.
# Also generate a sample metadata sheet for downstream analyses.

# First, organize the appropriate metadata files.
project="flu-Xue-chronic"
dir="nobackup/${project}/raw"
SRARunTable="data/metadata/${project}/SraRunTable.txt"
numSRA="$(tail -n +2 ${SRARunTable} | wc -l | cut -f1 -d' ')"
runs="data/metadata/${project}/SRR.txt"
mkdir -p nobackup/${project}/raw
mkdir -p nobackup/${project}/trimmed
mkdir -p data/metadata/${project}

# Submit a batch job to download FASTQ files from the SRA.
mkdir -p nobackup/${project}/sge
if [ -z "$(ls ${dir})" ]; then
  echo "Download files from SRA."
  cut -f9 ${SRARunTable} | tail -n +2 > ${runs}
  qsub -cwd -N DownloadSRA \
    -t 1-${numSRA} -tc 250 \
    -o nobackup/${project}/sge/ -e nobackup/${project}/sge/ \
    ${BatchSRADownloadFASTQ} ${dir} \
    ${runs} ${clean}
fi

# Prepare dataset sample sheet.
echo "Generate sample data sheet."
while read sample SRR
do
  # Pad sample name with zeros if it is a number.
  samplefull=${sample}
  if [[ $sample =~ ^[0-9]+$ ]] ; then
   samplefull=`printf %03d $sample`
  fi
  echo "${dir}/${SRR}_1.fastq.gz ${dir}/${SRR}_2.fastq.gz nobackup/${project}/trimmed/${samplefull}_trimmed-R1.fastq.gz nobackup/${project}/trimmed/${samplefull}_trimmed-R2.fastq.gz ${samplefull} ${project}"
done < <(cut -f4,9 ${SRARunTable} | tail -n +2 | sort | uniq ) \
  > data/metadata/${project}/samplesheet.txt





