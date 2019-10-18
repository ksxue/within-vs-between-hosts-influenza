#!/bin/bash

# Script is meant to be run from the top level of the Github repository.
# Script calls variants in all of the studies, concatenates them into a single file,
# and annotates them with their subtype and study of origin.

# Input parameters and directory paths.
analysis="CallWithinHostVariants"
analysisdir="analysis/${analysis}"
dir="nobackup/${analysis}"
mkdir -p ${dir}
mkdir -p analysis/${analysis}/out
BatchCallVariantsHardFilter="scripts/CallVariants/Batch-CallVariants-HardFilter.sh"
BatchCallVariantsHardFilterLongitudinal="scripts/CallVariants/Batch-CallVariants-HardFilter-Longitudinal.sh"
BatchSummarizeCoverage="scripts/SummarizeCoverage/Batch-SummarizeCoverage.sh"

if [ ! -f ${analysisdir}/variants-annotated-0.005.data.gz ]; then

  # Call variants using a hard filter.
  # Classify a site as variable if there are at least 400 reads at that site
  # and at least one of them supports a non-consensus base.
  projects=( "flu-Debbink" "flu-Dinis" "flu-McCrone" "flu-Xue-chronic" "flu-Barbezange" "flu-Xue-acute" )
  infectiontypes=( "acute" "acute" "acute" "chronic" "acute" "acute" )
  for i in $(seq 0 $((${#projects[@]}-1)) )
  do
    project=${projects[$i]}
	infectiontype=${infectiontypes[$i]}
    mkdir -p ${dir}/${project}
    samplesheet="data/metadata/${project}/samplesheet.txt"
    numsamples="$(wc -l ${samplesheet} | cut -f1 -d' ')"
    if [ ! "$(ls ${dir}/${project}/*.variants | wc -l)" -eq ${numsamples} ] && \
      [ -f data/metadata/${project}/samplesheet-refs.txt ]; then
	  echo "Calling variants, project ${project}, ${infectiontype} type"
	  if [ "${infectiontype}" == "acute" ]; then
        # Submit batch jobs to call variants in each sample individually.
        qsub -cwd -N CallVariantsHardFilter -l h_rt=48:00:00 \
          -t 1-${numsamples} -tc 250 \
          -o nobackup/sge -e nobackup/sge \
          ${BatchCallVariantsHardFilter} \
          data/metadata/${project}/samplesheet-refs.txt \
	      nobackup/DownloadDataCallVariants/${project} ${dir}/${project} 0 400
	  elif [ "${infectiontype}" == "chronic" ]; then
	    # Submit batch jobs to call variants in each sample individually.
		# Use the longitudinal variant of this script that calls variants
		# relative to the reference (original timepoint) sequence rather
		# than to the sample consensus.
        qsub -cwd -N CallVariantsHardFilter -l h_rt=48:00:00 \
          -t 1-${numsamples} -tc 250 \
          -o nobackup/sge -e nobackup/sge \
          ${BatchCallVariantsHardFilterLongitudinal} \
          data/metadata/${project}/samplesheet-refs.txt \
	      nobackup/DownloadDataCallVariants/${project} ${dir}/${project} 0 400
	  fi
    fi
  done

  # Concatenate all variants from each project.
  # Add a project identifier as the last column of each
  # variant call file.
  projects=( "flu-Debbink" "flu-Dinis" "flu-McCrone" "flu-Xue-chronic" )
  for i in $(seq 0 $((${#projects[@]}-1)) )
  do
    project=${projects[$i]}
	if [ ! -f ${dir}/${project}-annotated.variants.gz ]; then
      echo "Processing variants: ${project}"
	  cat ${dir}/${project}/*.variants > ${dir}/${project}.variants
      sed -i "s/$/\t${project}/g" ${dir}/${project}.variants
	  # Process the sequenced variants to annotate them with their study and subtype.
	  Rscript analysis/${analysis}/ProcessVariants.R \
	    ${dir}/${project}.variants ${dir}/${project}-annotated.variants ${project}
	  gzip ${dir}/${project}-annotated.variants
	  rm -f ${dir}/${project}.variants
	fi
  done

  # Concatenate all variants across all projects.
  if [ ! -f nobackup/${analysis}/variants-annotated.data.gz ]; then
    zcat ${dir}/flu-Debbink-annotated.variants.gz | head -n 1 \
	  > nobackup/${analysis}/variants-annotated.data
    for i in $(seq 0 $((${#projects[@]}-1)) )
    do
      project=${projects[$i]}
      echo "Concatenating variants: ${project}"
	  zcat ${dir}/${project}-annotated.variants | tail -n +2 \
	    >> nobackup/${analysis}/variants-annotated.data
	done
	gzip nobackup/${analysis}/variants-annotated.data
  fi
  
  # Generate a reduced file containing only variants >0.5% frequency.
  if [ ! -f nobackup/${analysis}/variants-annotated-0.005.data.gz ]; then
    zcat nobackup/${analysis}/variants-annotated.data.gz | \
	  awk '$18 > 0.005' > nobackup/${analysis}/variants-annotated-0.005.data
	gzip nobackup/${analysis}/variants-annotated-0.005.data
  fi
  
  # Copy this file to a location that is tracked by Github.
  cp nobackup/${analysis}/variants-annotated-0.005.data.gz \
    ${analysisdir}/variants-annotated-0.005.data.gz
fi

# Summarize sequencing coverage for all samples
# and concatenate coverage summaries into a single file.
if [ ! -f ${analysisdir}/coverage-annotated.data.gz ]; then
  # Summarize sequencing coverage.
  # Assess what proportion of sites have at least 400x coverage.
  projects=( "flu-Debbink" "flu-Dinis" "flu-McCrone" "flu-Xue-chronic" "flu-Barbezange" "flu-Xue-acute" )
  for i in $(seq 0 $((${#projects[@]}-1)) )
  do
    project=${projects[$i]}
    mkdir -p ${dir}/${project}
    samplesheet="data/metadata/${project}/samplesheet.txt"
    numsamples="$(wc -l ${samplesheet} | cut -f1 -d' ')"
    if [ ! "$(ls ${dir}/${project}/*.coverage | wc -l)" -eq ${numsamples} ] && \
      [ -f data/metadata/${project}/samplesheet-refs.txt ]; then
	  echo "Summarizing coverage, project ${project}"
	  # Submit batch jobs to call variants in each sample individually.
        qsub -cwd -N SummarizeCoverage -l h_rt=48:00:00 \
          -t 1-${numsamples} -tc 250 \
          -o nobackup/sge -e nobackup/sge \
          ${BatchSummarizeCoverage} \
          data/metadata/${project}/samplesheet-refs.txt \
	      nobackup/DownloadDataCallVariants/${project} ${dir}/${project} 400
    fi
  done
  
  # Concatenate all variants from each project.
  # Add a project identifier as the last column of each
  # variant call file.
  projects=( "flu-Debbink" "flu-Dinis" "flu-McCrone" "flu-Xue-chronic" )
  for i in $(seq 0 $((${#projects[@]}-1)) )
  do
    project=${projects[$i]}
	if [ ! -f ${dir}/${project}-annotated.coverage.gz ]; then
      echo "Processing coverage summaries: ${project}"
	  cat ${dir}/${project}/*.coverage > ${dir}/${project}.coverage
      sed -i "s/$/\t${project}/g" ${dir}/${project}.coverage
	  # Process the coverage summaries to annotate them with their study and subtype.
	  Rscript analysis/${analysis}/ProcessCoverageSummaries.R \
	    ${dir}/${project}.coverage ${dir}/${project}-annotated.coverage ${project}
	  gzip ${dir}/${project}-annotated.coverage
	  rm -f ${dir}/${project}.coverage
	fi
  done

  # Concatenate all coverage summaries across all projects.
  if [ ! -f nobackup/${analysis}/coverage-annotated.data.gz ]; then
    zcat ${dir}/flu-Debbink-annotated.coverage.gz | head -n 1 \
	  > nobackup/${analysis}/coverage-annotated.data
    for i in $(seq 0 $((${#projects[@]}-1)) )
    do
      project=${projects[$i]}
      echo "Concatenating coverage summaries: ${project}"
	  zcat ${dir}/${project}-annotated.coverage | tail -n +2 \
	    >> nobackup/${analysis}/coverage-annotated.data
	done
	gzip nobackup/${analysis}/coverage-annotated.data
  fi
  
  
  # Copy this file to a location that is tracked by Github.
  cp nobackup/${analysis}/coverage-annotated.data.gz \
    ${analysisdir}/coverage-annotated.data.gz
fi