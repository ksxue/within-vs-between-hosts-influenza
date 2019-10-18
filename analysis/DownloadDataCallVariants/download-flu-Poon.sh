#!/bin/bash

# Script is meant to be run from the top level of the Github repository.
# Script organizes raw sequencing files from the SRA for Poon et al. 2016.


# Poon et al. 2016, Nature Genetics
# Raw sequencing data is available on the Synapse servers at
# https://www.synapse.org/#!Synapse:syn8033988
# Access requires registering for a free account.

# I manually downloaded the sequencing reads, README,
# sample annotations, and other associated files
# into the nobackup/flu-Poon/raw folder
# on April 23, 2018.
# I also manually converted the Excel file
# hongkong_transmission_pairs.xlsx to a text file,
# hongkong_transmission_pairs.txt and removed the first
# header line to facilitate subsequent analyses.

# First, organize the appropriate directories.
project="flu-Poon"
dir="nobackup/${project}/raw/rawfiles"
mkdir -p nobackup/${project}/raw
mkdir -p data/metadata/${project}


# Unzip the raw sequencing reads.
if [ -f nobackup/flu-Poon/raw/hongkong_rawdata.tar.gz ]; then
  echo "Unzipping sequencing reads."
  tar -xzvf nobackup/flu-Poon/raw/hongkong_rawdata.tar.gz \
    -C nobackup/flu-Poon/raw/
  rm -f nobackup/flu-Poon/raw/hongkong_rawdata.tar.gz
fi

# Copy metadata to tracked directory.
cp nobackup/flu-Poon/raw/bacid_blindednumber_virusname_filename.txt \
  data/metadata/flu-Poon/Synapse-raw.txt
cp nobackup/flu-Poon/raw/hongkong_transmission_pairs.txt \
  data/metadata/flu-Poon/Synapse-transmissions.txt

# Parse the metadata to create a list of strain and sample names.
paste <( cut -f4 -d, \
  data/metadata/flu-Poon/Synapse-raw.txt ) \
  <( cut -f3 -d, \
  data/metadata/flu-Poon/Synapse-raw.txt | \
  cut -f3 -d'/' | cut -f2,3 -d'-' ) \
  <( cut -f4 -d, \
  data/metadata/flu-Poon/Synapse-raw.txt | \
  cut -f5 -d_ ) \
  > data/metadata/flu-Poon/Synapse-samplenames.txt

# Create a samplesheet for the reads in single-end format.
# Note that the reads as downloaded were in single-end format.
while read fastq strain ID
do
  sample="${strain}-${ID}"
  echo ${dir}/${fastq} ${dir}/${fastq} \
    ${dir}/${fastq} ${dir}/${fastq} \
	${sample} ${project}
done < data/metadata/flu-Poon/Synapse-samplenames.txt \
  > data/metadata/flu-Poon/samplesheet-SingleEnd.txt
  
# Some samples appear to have sequencing replicates, and others do not.
# Using the available metadata, create a list of samples with replicates
# and a list of samples without replicates.
while read strain
do
  grep ${strain} data/metadata/flu-Poon/Synapse-samplenames.txt \
    | cut -f2,3 | tr '\t' '-'
done < <( cut -f2 data/metadata/flu-Poon/Synapse-samplenames.txt | \
  sort | uniq -c | awk '$1==2{print $2}' ) \
  > data/metadata/flu-Poon/Synapse-samplenames-replicates.txt
  
while read strain
do
  grep ${strain} data/metadata/flu-Poon/Synapse-samplenames.txt \
    | cut -f2,3 | tr '\t' '-'
done < <( cut -f2 data/metadata/flu-Poon/Synapse-samplenames.txt | \
  sort | uniq -c | awk '$1==1{print $2}' ) \
  > data/metadata/flu-Poon/Synapse-samplenames-noreplicates.txt