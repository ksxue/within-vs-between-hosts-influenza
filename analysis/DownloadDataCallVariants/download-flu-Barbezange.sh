#!/bin/bash

# Script is meant to be run from the top level of the Github repository.
# Script organizes raw sequencing files for Barbezange et al. 2018.
# https://www.frontiersin.org/articles/10.3389/fmicb.2018.02596/full

# Input paths.
BatchSRADownloadFASTQ="pipelines/Batch-SRADownloadFASTQ.sh"
ConsolidateRenameSRA="analysis/DownloadDataCallVariants/ConsolidateRenameSRA-flu-Barbezange.sh"
clean=1 # If 0, analysis will run from beginning, overwriting existing intermediates.

# Barbezange et al. 2018, Frontiers of Microbiology
# Download raw FASTQ files from the SRA Bioproject PRJEB11406.
# Use the RunInfo Table for this BioProject to find SRA accessions.
# Rename files to reflect biological sample names.
# Also generate a sample metadata sheet for downstream analyses.

# First, organize the appropriate metadata files.
project="flu-Barbezange"
dir="nobackup/${project}/raw"
SRARunTable="data/metadata/${project}/SraRunTable.txt"
numSRA="$(tail -n +2 ${SRARunTable} | wc -l | cut -f1 -d' ')"
runs="data/metadata/${project}/SRR.txt"
mkdir -p nobackup/${project}/raw
mkdir -p nobackup/${project}/sge
mkdir -p nobackup/${project}/trimmed
mkdir -p data/metadata/${project}

# Submit a batch job to download FASTQ files from the SRA.
mkdir -p nobackup/${project}/sge
if [ -z "$(ls ${dir})" ]; then
  echo "Download files from SRA."
  cut -f8 ${SRARunTable} | tail -n +2 > ${runs}
  qsub -cwd -N DownloadBarbezange \
    -t 1-${numSRA} -tc 250 \
    -o nobackup/${project}/sge/ -e nobackup/${project}/sge/ \
    ${BatchSRADownloadFASTQ} ${dir} \
    ${runs} ${clean}
fi

# Rename files to reflect sample names.
# In this case, each sample had exactly one corresponding run.
if [ "$(ls ${dir}/*.fastq.gz | wc -l)" -eq $(( numSRA )) ]; then
  
  echo "Rename files and consolidate FASTQs from multiple runs."
  ${ConsolidateRenameSRA} ${SRARunTable} ${dir}
	
  # Remove all remaining SRR* files.
  rm -f ${dir}/SRR*.fastq.gz
fi

# Prepare dataset sample sheet.
echo "Generate sample data sheet."
if [ ! -f data/metadata/${project}/samplesheet.txt ]; then
  while read sample
  do
    # Pad sample name with zeros if it is a number.
    samplefull=${sample}
    if [[ $sample =~ ^[0-9]+$ ]] ; then
      samplefull=`printf %03d $sample`
    fi
    echo "${dir}/${sample}_1.fastq.gz ${dir}/${sample}_2.fastq.gz ${dir}/${sample}_1.fastq.gz ${dir}/${sample}_2.fastq.gz Barbezange-${samplefull} ${project}"
  done < <(cut -f16 ${SRARunTable} | tail -n +2 | sort | uniq ) \
    > data/metadata/${project}/samplesheet.txt
fi



