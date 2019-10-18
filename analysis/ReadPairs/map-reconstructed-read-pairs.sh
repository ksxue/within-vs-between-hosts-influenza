#!/bin/bash

# Script is meant to be run from the top level of the Github repository.
# Script maps reconstructed read pairs from pairing analysis.

# Input parameters and directory paths.
read1=$1 # Takes a FASTQ file containing read 1 of each pair.
read2=$2 # Takes a FASTQ file containing read 2 of each pair.
outdir=$3 # Output directory.

# Create a reference genome that combines the H3N2 and pdmH1N1 references.
# Map all paired reads to this pan-reference genome
# to determine what proportion align concordantly, i.e. are likely to be
# true read pairs.
mkdir -p ${outdir}/reference
# Concatenate reference genomes.
cat <( sed 's/>/>H3N2-/g' reference/flu-H3N2/H3N2-Brisbane-2007.fasta ) \
  <( sed 's/>/>pdmH1N1-/g' reference/flu-pdmH1N1/pdmH1N1-California-2009.fasta ) \
  > ${outdir}/reference/H3N2-pdmH1N1.fasta
# Build the bowtie2 index.
bowtie2-build ${outdir}/reference/H3N2-pdmH1N1.fasta \
  ${outdir}/reference/H3N2-pdmH1N1
  
# Map all reads.
mkdir -p ${outdir}/alignment
bowtie2 --very-fast-local --un-conc-gz ${outdir}/alignment/paired-unmapped \
  -X 2300 \
  -x ${outdir}/reference/H3N2-pdmH1N1 \
  -p 16 \
  -1 ${read1} \
  -2 ${read2} \
  -S ${outdir}/alignment/paired_all.sam \
  2> ${outdir}/alignment/paired_all.bt2.log
samtools view -bS -@ 16 ${outdir}/alignment/paired_all.sam \
  ${outdir}/alignment/paired_all.bam
