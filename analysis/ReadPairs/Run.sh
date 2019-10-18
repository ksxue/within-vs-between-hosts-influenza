#!/bin/bash

# Script is meant to be run from the top level of the Github repository.
# Script analyzes read pairing within samples in the Poon et al. raw data.

# Input parameters and directory paths.
project="flu-Poon"
analysis="ReadPairs"
rawdir="nobackup/flu-Poon/raw"
outdir="nobackup/${analysis}/"
mkdir -p ${outdir}
mkdir -p analysis/${analysis}/out

# For each FASTQ file, sort the reads based on read header.
# Then merge-sort all reads from the study into a single FASTQ file.
if [ ! -f ${outdir}/sorted_all.fastq.gz ]; then
  qsub -cwd -N MergeSortReads -o nobackup/sge -e nobackup/sge \
    -l m_mem_free=4G \
    analysis/ReadPairs/merge-sort-reads.sh
fi

# Iterate through the sorted reads.
# Identify paired reads and output them onto the same line.
if [ ! -f ${outdir}/paired_all.fastq.gz ] && \
  [ -f ${outdir}/sorted_all.fastq.gz ];
then
  qsub -cwd -N ExtractReadPairs -o nobackup/sge -e nobackup/sge \
  -l m_mem_free=4G \
  analysis/ReadPairs/extract-paired-reads.sh ${outdir}/sorted_all.fastq.gz \
    ${outdir}/paired_all.fastq
fi

# Iterate through the sorted reads.
# Identify unpaired reads and output them onto the same line.
if [ ! -f ${outdir}/single_all.fastq.gz ] && \
  [ -f ${outdir}/sorted_all.fastq.gz ];
then
  qsub -cwd -N ExtractUnpairedReads -o nobackup/sge -e nobackup/sge \
  analysis/ReadPairs/extract-single-reads.sh ${outdir}/sorted_all.fastq.gz \
    ${outdir}/single_all.fastq
fi

# Extract and count read IDs and run IDs for reconstructed read pairs.
# Also extract read pairs and format them in standard FASTQ format.
if [ ! -f ${outdir}/read2_all.fastq.gz ] && \
  [ -f ${outdir}/paired_all.fastq.gz ];
then
  qsub -cwd -N AnalyzeReadPairs -o nobackup/sge -e nobackup/sge \
  analysis/ReadPairs/analyze-read-pairs.sh ${outdir}/paired_all.fastq.gz \
    ${outdir}
fi

# Extract and count read IDs and run IDs for single-end reads.
if [ ! -f ${outdir}/single_counts.txt ] && \
  [ -f ${outdir}/single_all.fastq.gz ];
then
  qsub -cwd -N AnalyzeSingleReads -o nobackup/sge -e nobackup/sge \
  analysis/ReadPairs/analyze-single-reads.sh ${outdir}/single_all.fastq.gz \
    ${outdir}
fi

# Map the reconstructed read pairs.
if [ ! -f ${outdir}/alignment/paired_all.bt2.log ] && \
  [ -f ${outdir}/read2_all.fastq.gz ];
then
  qsub -cwd -N MapReadPairs -o nobackup/sge -e nobackup/sge \
    -pe serial 16 \
    analysis/ReadPairs/map-reconstructed-read-pairs.sh \
    ${outdir}/read1_all.fastq.gz ${outdir}/read2_all.fastq.gz \
    ${outdir}
fi

# Store the tallies of paired-end and single-end reads
# in a tracked directory.
mkdir -p analysis/${analysis}/out/
if [ -f ${outdir}/single_counts.txt ] && \
  [ -f ${outdir}/readIDrun_counts.txt ];
then
  cp -f ${outdir}/single_counts.txt analysis/${analysis}/out/
  cp -f ${outdir}/readIDrun_counts.txt analysis/${analysis}/out/
fi

# Summarize the number of paired and unpaired reads.
if [ -f ${outdir}/single_counts.txt ] && \
  [ -f ${outdir}/readIDrun_counts.txt ] && \
  [ ! -f analysis/${analysis}/ReadSummary.txt ];
then
  Rscript analysis/ReadPairs/AnalyzeReadPairs.R
fi