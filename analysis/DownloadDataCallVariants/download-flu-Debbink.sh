#!/bin/bash

# Script is meant to be run from the top level of the Github repository.
# Script organizes raw sequencing files for Debbink et al. 2017.

# Input paths.
BatchSRADownloadFASTQ="pipelines/Batch-SRADownloadFASTQ.sh"
ConsolidateRenameSRA="scripts/ConsolidateRenameSRA.sh"
clean=1 # If 0, analysis will run from beginning, overwriting existing intermediates.

# Debbink et al. 2017, PLoS Pathogens
# Download raw FASTQ files from the SRA Bioproject PRJNA344659.
# Use the RunInfo Table for this BioProject to find SRA accessions.
# Rename files to reflect biological sample names.
# Also generate a sample metadata sheet for downstream analyses.

# First, organize the appropriate metadata files.
project="flu-Debbink"
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
  cut -f10 ${SRARunTable} | tail -n +2 > ${runs}
  qsub -cwd -N DownloadDebbink \
    -t 1-${numSRA} -tc 250 \
    -o nobackup/${project}/sge/ -e nobackup/${project}/sge/ \
    ${BatchSRADownloadFASTQ} ${dir} \
    ${runs} ${clean}
fi

# Rename files to reflect sample names.
# If multiple runs were performed for a single sample,
# then identify those samples and merge the resultant FASTQ files.
if [ "$(ls ${dir}/*.fastq.gz | wc -l)" -eq $(( numSRA * 2 )) ]; then
  
  echo "Rename files and consolidate FASTQs from multiple runs."
  ${ConsolidateRenameSRA} ${SRARunTable} ${dir}

  # Manually combine the remaining exception samples
  # that have three or more corresponding sequencing runs.
  # Then, clear remaining SRR files.
  cat ${dir}/SRR4309567_1.fastq.gz ${dir}/SRR4309568_1.fastq.gz \
    ${dir}/SRR4309569_1.fastq.gz > ${dir}/Brisbane_1.fastq.gz
  cat ${dir}/SRR4309567_2.fastq.gz ${dir}/SRR4309568_2.fastq.gz \
    ${dir}/SRR4309569_2.fastq.gz > ${dir}/Brisbane_2.fastq.gz
  cat ${dir}/SRR4309570_1.fastq.gz ${dir}/SRR4309571_1.fastq.gz \
    ${dir}/SRR4309572_1.fastq.gz > ${dir}/Cal-H3N2_1.fastq.gz
  cat ${dir}/SRR4309570_2.fastq.gz ${dir}/SRR4309571_2.fastq.gz \
    ${dir}/SRR4309572_2.fastq.gz > ${dir}/Cal-H3N2_2.fastq.gz
  
  rm -f ${dir}/SRR*
fi

# Prepare dataset sample sheet.
echo "Generate sample data sheet."
while read sample
do
  # Pad sample name with zeros if it is a number.
  samplefull=${sample}
  if [[ $sample =~ ^[0-9]+$ ]] ; then
   samplefull=`printf %03d $sample`
  fi
  echo "${dir}/${sample}_1.fastq.gz ${dir}/${sample}_2.fastq.gz nobackup/${project}/trimmed/${samplefull}_trimmed-R1.fastq.gz nobackup/${project}/trimmed/${samplefull}_trimmed-R2.fastq.gz ${samplefull} ${project}"
done < <(cut -f12 ${SRARunTable} | tail -n +2 | sort | uniq ) \
  > data/metadata/${project}/samplesheet.txt
  
# Download viral load values.
wget -O data/metadata/${project}/TitersStatus_2004-2005.csv \
  https://raw.githubusercontent.com/lauringlab/Fluvacs_paper/master/Titers_status_2004-2005.csv
wget -O data/metadata/${project}/TitersStatus_2005-2006.csv \
  https://raw.githubusercontent.com/lauringlab/Fluvacs_paper/master/Titers_status_2005-2006.csv
wget -O data/metadata/${project}/TitersStatus_2007-2008.csv \
  https://raw.githubusercontent.com/lauringlab/Fluvacs_paper/master/Titers_status_2007-2008.csv
# Append season to other metadata.
for f in data/metadata/${project}/TitersStatus*.csv
do
  season=${f##*_}
  season=${season%%.*}
  sed -i "s/$/,${season}/" ${f}
done
# Merge files for each season.
cat <( tail -n +2 data/metadata/${project}/TitersStatus*.csv ) > \
  data/metadata/${project}/TitersStatus.txt
# Format files by removing blank lines and spacer text.
sed -i '/^$/d' data/metadata/${project}/TitersStatus.txt
sed -i '/^=/d' data/metadata/${project}/TitersStatus.txt
# Add file header.
sed -i '1s/^/Sample CopyNumber VaccinationStatus Season\n/' \
  data/metadata/${project}/TitersStatus.txt
sed -i 's/,/ /g' data/metadata/${project}/TitersStatus.txt
# Remove intermediate files.
rm -f data/metadata/${project}/TitersStatus*.csv

# Download full sample metadata.
wget -O data/metadata/${project}/metadata-${project}.data \
  https://raw.githubusercontent.com/lauringlab/Fluvacs_paper/master/data/raw/all.meta.csv




