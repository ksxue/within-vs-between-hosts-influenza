#!/bin/bash

# Script is meant to be run from the top level of the Github repository.
# Script aligns and calls variants in each of the sequenced samples
# for the Xue et al. 2017 eLife dataset.


# Xue et al. 2017 eLife
# Here, to analyze these sequencing reads,
# map the first 1000 reads for each sample against reference flu genomes.
# Trim the reads to remove any remaining Nextera adapters.
# Next, use bowtie2 to align all reads against the reference genome.
# Infer the consensus sequence for each viral sample,
# and remap reads for each sample against the consensus genome for that sample.
# Remove duplicate sequencing reads.
# Then, summarize the number of reads supporting each nucleotide identity
# at each genome position, and annotate the variants with their effects on
# the protein sequence.
project="flu-Xue-chronic"
samplesheet="data/metadata/${project}/samplesheet.txt"
numsamples="$(wc -l ${samplesheet} | cut -f1 -d' ')"
BatchInferSubtype="scripts/InferSubtype/Batch-InferSubtype.sh"
BatchAlignSummarizeRealignLongitudinal="pipelines/Batch-AlignSummarizeRealign-Longitudinal.sh"
BatchCallVariantsHardFilterLongitudinal="scripts/CallVariants/Batch-CallVariants-HardFilter-Longitudinal.sh"
BatchTrimAdaptersNextera="pipelines/Batch-TrimAdaptersNextera.sh"

# Create output folder for intermediate files.
analysis="DownloadDataCallVariants"
dir="nobackup/${analysis}/${project}"
mkdir -p ${dir}
mkdir -p nobackup/${project}/trimmed

# Infer the subtype of each flu sample.
if [ ! "$(ls ${dir}/*.map | wc -l)" -eq ${numsamples} ]; then
  echo "Inferring sample subtypes."
  # Submit batch jobs to analyze the subtypes of each flu sample.
  qsub -cwd -N InferSubtype -l h_rt=1:00:00 \
    -hold_jid InferSubtype \
    -t 1-${numsamples} -tc 250 \
    -o nobackup/sge -e nobackup/sge \
    ${BatchInferSubtype} ${samplesheet} ${dir} \
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

# Because this is a longitudinal dataset, we want to call variants
# relative to the patient consensus sequence at the first timepoint.
# Create a version of the samplesheet with this information.
if [ ! -f data/metadata/${project}/samplesheet-refs-timepoint.txt ]; then
  echo "Create samplesheet with initial timepoint reference genomes."
  while read fastq1 fastq2 trimmed1 trimmed2 sample study reference refname
  do
    # Parse the sample names and determine which timepoint is the initial one.
	patient="$(echo $sample | head -c 1)"
	replicate="${sample: -1}"
	timepoint="00"
	while read refpatient reftimepoint
	do
	  if [ "$patient" == "$refpatient" ]; then
	    timepoint=$reftimepoint
	  fi
	done < data/metadata/${project}/InitialTimepoints.txt
	echo ${fastq1} ${fastq2} ${trimmed1} ${trimmed2} \
	    ${sample} ${study} ${reference} ${refname} "${patient}${timepoint}A-NW-${replicate}-${refname}"
  done < data/metadata/${project}/samplesheet-refs.txt \
    > data/metadata/${project}/samplesheet-refs-timepoint.txt
fi


# Trim Nextera adapters from sequencing reads.
if [ -z "$(ls nobackup/${project}/trimmed/)" ] && \
  [ "$(ls nobackup/${project}/raw/*.fastq.gz | wc -l)" -eq $(( 2*numsamples )) ] && \
  [ -f data/metadata/${project}/samplesheet-refs.txt ]; then
  echo "Trimming reads."
  # Submit batch jobs to trim Nextera adapters for each flu samples.
  qsub -cwd -N TrimReads -l h_rt=48:00:00 \
    -t 1-${numsamples} -tc 250 \
    -o nobackup/sge -e nobackup/sge \
    ${BatchTrimAdaptersNextera} \
    ${samplesheet} nobackup/${project}/trimmed 0
fi

# Align reads to reference genome, infer sample consensus,
# re-align reads to consensus genome, and tally and annotate variants.
if [ ! "$(ls ${dir}/*-realigned-annotated.summary | wc -l)" -eq ${numsamples} ] && \
  [ "$(ls nobackup/${project}/trimmed/*.fastq.gz | wc -l)" -eq $(( 2*numsamples )) ] && \
  [ -f data/metadata/${project}/samplesheet-refs.txt ]; then
  echo "Aligning reads and calling variants."
  # Submit batch jobs to map all reads for each flu samples.
  # Map reads and call variants based on initial mapping to the
  # best-match reference genome determined above.
  qsub -cwd -N AlignReads -l h_rt=96:00:00 \
    -t 1-${numsamples} -tc 250 -l m_mem_free=12G \
    -o nobackup/sge -e nobackup/sge \
    ${BatchAlignSummarizeRealignLongitudinal} \
    data/metadata/${project}/samplesheet-refs-timepoint.txt ${dir} 1  
fi

# Call variants using a hard filter.
# Classify a site as variable if there are at least 200 reads at that site
# and at least 3% of them support a non-consensus base.
mkdir -p ${dir}/variants
if [ ! "$(ls ${dir}/variants/*.variants | wc -l)" -eq ${numsamples} ] && \
  [ -f data/metadata/${project}/samplesheet-refs.txt ] && \
  [ "$(ls ${dir}/*-realigned-annotated.summary | wc -l)" -eq ${numsamples} ]; then
  echo "Calling variants."
  # Submit batch jobs to call variants in each sample individually.
  qsub -cwd -N CallVariantsHardFilter -l h_rt=48:00:00 \
    -t 1-${numsamples} -tc 250 \
    -o nobackup/sge -e nobackup/sge \
    ${BatchCallVariantsHardFilterLongitudinal} \
    data/metadata/${project}/samplesheet-refs.txt \
	${dir} ${dir}/variants 0.005 200 
fi

# Determine which samples do not produce variant files,
# i.e. presumably have no variable sites under the given criteria.
if [ ! "$(ls ${dir}/variants/*.variants | wc -l)" -eq ${numsamples} ] && \
  [ ! -f ${dir}/variants/novariants.txt ] && \
  [ "$(ls ${dir}/variants/*.variants | wc -l)" > 0 ] && \
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

