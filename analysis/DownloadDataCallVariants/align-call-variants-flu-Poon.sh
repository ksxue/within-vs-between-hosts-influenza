#!/bin/bash

# Script is meant to be run from the top level of the Github repository.
# Script aligns and calls variants in each of the sequenced samples
# for the Poon et al. 2016 Nature Genetics dataset.


# Poon et al. 2016, Nature Genetics
# Treat all reads as unpaired for this analysis.
# First, map the first 1000 reads for each sample against reference flu genomes
# to infer the best-match sample subtype.
# Next, use bowtie2 to align all reads against the reference genome.
# Infer the consensus sequence for each viral sample,
# and remap reads for each sample against the consensus genome for that sample.
# Then, summarize the number of reads supporting each nucleotide identity
# at each genome position, and annotate the variants with their effects on
# the protein sequence.
# In contrast to the standard analysis pipeline for paired-end reads,
# this single-end pipeline does not trim sequencing adapters
# (which are unknown for this dataset), and reads are not deduplicated
# because duplication rate is extremely high in the absence of read-pair data.
project="flu-Poon"
samplesheet="data/metadata/${project}/samplesheet-SingleEnd.txt"
numsamples="$(wc -l ${samplesheet} | cut -f1 -d' ')"
BatchInferSubtypeSingleEnd="scripts/InferSubtype/Batch-InferSubtype-SingleEnd.sh"
BatchAlignSummarizeRealignSingleEnd="pipelines/Batch-AlignSummarizeRealign-SingleEnd.sh"
BatchCallVariantsHardFilter="scripts/CallVariants/Batch-CallVariants-HardFilter.sh"

# Create output folder for intermediate files.
analysis="DownloadDataCallVariants"
dir="nobackup/${analysis}/${project}"
mkdir -p ${dir}

# Infer the subtype of each flu sample.
if [ ! "$(ls ${dir}/*.map | wc -l)" -eq ${numsamples} ]; then
  echo "Inferring sample subtypes."
  # Submit batch jobs to analyze the subtypes of each flu sample.
  qsub -cwd -N InferSubtype -l h_rt=1:00:00 \
    -hold_jid InferSubtype \
    -t 1-${numsamples} -tc 250 \
    -o nobackup/sge -e nobackup/sge \
    ${BatchInferSubtypeSingleEnd} ${samplesheet} ${dir} \
    scripts/InferSubtype/flu-subtypes.data 1000
fi

# Generate a samplesheet with reference genomes.
if [ "$(ls ${dir}/*.map | wc -l)" -eq ${numsamples} ] && \
  [ ! -f data/metadata/${project}/samplesheet-refs.txt ]; then
  echo "Generate sample sheet with inferred subtypes."
  # Generate a summary of the best-matching reference genomes
  # for each sample and the mapping rates to these genomes.
  while read fastq1 fastq2 trimmed1 trimmed2 sample study
  do
    head -n 1 ${dir}/${sample}.map
  done < ${samplesheet} > ${dir}/${project}-subtypes.summary
  cut -f3,5,6 -d' ' ${dir}/${project}-subtypes.summary \
    > analysis/${analysis}/${project}-subtypes.data

  # Merge information from the samplesheet and about sample subtypes.
  while read fastq1 fastq2 trimmed1 trimmed2 sample study \
    fastq12 fastq22 sample2 reference refname maprate
  do
    if [ "$sample" == "$sample2" ]
    then
      echo ${fastq1} ${fastq2} ${trimmed1} ${trimmed2} \
	    ${sample} ${study} ${reference} ${refname}
    else
      echo "Error in matching subtype to sample."
    fi
  done < <( paste -d' ' ${samplesheet} ${dir}/${project}-subtypes.summary ) \
    > ${dir}/${project}.samplesheet
	
  # Copy samplesheet with reference genomes to metadata folder.
  cp ${dir}/${project}.samplesheet \
    data/metadata/${project}/samplesheet-refs.txt
fi

# Align reads to reference genome, infer sample consensus,
# re-align reads to consensus genome, and tally and annotate variants.
if [ ! "$(ls ${dir}/*-realigned-annotated.summary | wc -l)" -eq ${numsamples} ] && \
  [ -f data/metadata/${project}/samplesheet-refs.txt ]; then
  echo "Aligning reads and calling variants."
  # Submit batch jobs to map all reads for each flu samples.
  # Map reads and call variants based on initial mapping to the
  # best-match reference genome determined above.
  qsub -cwd -N AlignReads -l h_rt=48:00:00 \
    -t 1-${numsamples} -tc 250 \
    -o nobackup/sge -e nobackup/sge \
    ${BatchAlignSummarizeRealignSingleEnd} \
    ${dir}/${project}.samplesheet ${dir} 1
fi

# Call variants using a hard filter.
# Classify a site as variable if there are at least 200 reads at that site
# and at least 3% of them support a non-consensus base.
if [ "$(ls ${dir}/variants/*.variants | wc -l)" -eq 0 ] && \
  [ "$(ls ${dir}/*-realigned-annotated.summary | wc -l)" -eq ${numsamples} ] && \
  [ -f data/metadata/${project}/samplesheet-refs.txt ]; then
  echo "Calling variants."
  mkdir -p ${dir}/variants
  # Submit batch jobs to call variants in each sample individually.
  qsub -cwd -N CallVariantsHardFilter -l h_rt=48:00:00 \
    -t 1-${numsamples} -tc 250 \
    -o nobackup/sge -e nobackup/sge \
    ${BatchCallVariantsHardFilter} \
    data/metadata/${project}/samplesheet-refs.txt ${dir} \
    ${dir}/variants 0.03 200 
fi

# Determine which samples do not produce variant files,
# i.e. presumably have no variable sites under the given criteria.
if [ ! "$(ls ${dir}/variants/*.variants | wc -l)" -eq ${numsamples} ] && \
  [ "$(ls ${dir}/variants/*.variants | wc -l)" > 0 ] && \
  [ ! -f ${dir}/variants/novariants.txt ] && \
  [ -f data/metadata/${project}/samplesheet-refs.txt ]; then
  echo "Identifying samples with no variable sites."
  while read fastq1 fastq1 trimmed1 trimmed2 sample project refpath ref
  do
    if [ ! -f ${dir}/variants/${sample}-${ref}.variants ]; then
	  echo ${sample}
	fi
  done < data/metadata/${project}/samplesheet-refs.txt \
    > ${dir}/variants/novariants.txt
fi

