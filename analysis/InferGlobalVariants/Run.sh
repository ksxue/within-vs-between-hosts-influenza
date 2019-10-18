#!/bin/bash

# Script is meant to be run from the top level of the Github repository.

module load mafft/7.407

analysis="InferGlobalVariants"
RetainUnpassagedSequences="analysis/${analysis}/RetainUnpassagedSequences.R"
ExtractCodingSequences="analysis/${analysis}/ExtractCodingSequences.R"
CalculateDifferences="analysis/${analysis}/CalculateSequenceDifferences.R"
IdentifyOutlierSequences="analysis/${analysis}/IdentifyOutlierSequences.R"
IdentifyGISAIDGlobalVariants="analysis/${analysis}/IdentifyGISAIDGlobalVariants.R"

SplitSequencesBySeason="analysis/${analysis}/SplitSequencesBySeason.R"
BatchBuildTreeRAxML="pipelines/Batch-BuildTreeRAxML.sh"
RootTree="scripts/RootTree.py"
BatchReconstructAncestors="pipelines/Batch-ReconstructAncestors.sh"
BatchExtractMutationFrequencies="pipelines/Batch-ExtractMutationFrequencies.sh"
CalculateNSProportions="analysis/${analysis}/CalculateNSProportions"

mkdir -p nobackup/${analysis}
dir="nobackup/${analysis}"
mkdir -p analysis/${analysis}/out

##################################
# Analysis of H3N2 sequences.
##################################

# For each influenza gene, concatenate all downloaded sequences
# into a single file, analyze the passage history, and
# retain only sequences that are from unpassaged isolates.
seqdir="data/GISAID/H3N2"
# Genes=( "1-PB2" )
# Coding=( "1-PB2" )
Genes=( "1-PB2" "2-PB1" "3-PA" "4-HA" "5-NP" "6-NA" "7-M" "8-NS" )
Coding=( "1-PB2" "2-PB1" "3-PA" "4-HA" "5-NP" "6-NA" "7-M1" "8-NS1" )
NumGenes=${#Genes[@]}
for i in $(seq 0 $((${#Genes[@]}-1)) )
do
  # Concatenate all downloaded sequences into a single file.
  if [ ! -f ${dir}/H3N2-${Genes[$i]}.fasta ]; then
    echo "Concatenate sequences."
    echo ${Genes[$i]}
    cat ${seqdir}/H3N2-${Genes[$i]}-*.fasta > \
      ${dir}/H3N2-${Genes[$i]}.fasta
  fi
  
  # Retain only sequences that are marked as unpassaged.
  if [ ! -f ${dir}/H3N2-${Genes[$i]}-unpassaged.fasta ]; then
    echo "Remove passaged sequences."
	echo ${Genes[$i]}
    Rscript ${RetainUnpassagedSequences} \
	  ${dir}/H3N2-${Genes[$i]}.fasta ${dir}/H3N2-${Genes[$i]}-unpassaged.fasta
  fi
  
  # Concatenate the reference sequence for a gene with the unpassaged sequences.
  # Then, use MAFFT to build a multiple sequence alignment.
  if [ ! -f ${dir}/H3N2-${Coding[$i]}-unpassaged-ref.aligned ]; then
    echo "Align sequences."
	echo ${Genes[$i]}
	cat reference/flu-H3N2/H3N2-Victoria-2011-${Coding[$i]}-coding.fasta \
	  ${dir}/H3N2-${Genes[$i]}-unpassaged.fasta \
	  > ${dir}/H3N2-${Genes[$i]}-unpassaged-ref.fasta
    mafft ${dir}/H3N2-${Genes[$i]}-unpassaged-ref.fasta \
	  > ${dir}/H3N2-${Coding[$i]}-unpassaged-ref.aligned
	rm -f ${dir}/H3N2-${Genes[$i]}-unpassaged-ref.fasta
  fi
  
  # Extract the coding sequences based on the reference coding sequence.
  if [ ! -f ${dir}/H3N2-${Coding[$i]}-unpassaged-ref.aligned.coding ]; then
    echo "Extract coding sequences."
	echo ${Genes[$i]}
	Rscript ${ExtractCodingSequences} ${dir}/H3N2-${Coding[$i]}-unpassaged-ref.aligned \
	  ${dir}/H3N2-${Coding[$i]}-unpassaged-ref.aligned.coding
  fi
  
  # Calculate S and NS differences between each sequence and the reference coding sequence.
  if [ ! -f ${dir}/H3N2-${Coding[$i]}-unpassaged-ref.distances ]; then
    echo "Calculate sequence distances."
	echo ${Genes[$i]}
	Rscript ${CalculateDifferences} \
	  ${dir}/H3N2-${Coding[$i]}-unpassaged-ref.aligned.coding \
	  ${dir}/H3N2-${Coding[$i]}-unpassaged-ref.distances
  fi
  
  # Identify outlier sequences.
  if [ ! -f ${dir}/H3N2-${Coding[$i]}-unpassaged-ref.distances.annotated ]; then
    echo "Identify outlier sequences."
	echo ${Genes[$i]}
	Rscript ${IdentifyOutlierSequences} \
	  ${dir}/H3N2-${Coding[$i]}-unpassaged-ref.distances \
	  ${dir}/H3N2-${Coding[$i]}-unpassaged-ref.exclusions \
	  ${dir}/H3N2-${Coding[$i]}-unpassaged-ref.distances.annotated
  fi
  
  # Summarize the variant counts from the aligned FASTA file.
  # Remove outlier sequences, and only analyze coding variation.
  if [ ! -f ${dir}/H3N2-${Coding[$i]}-unpassaged.summary ]; then
    echo "Summarize global variants."
    echo ${Genes[$i]}
	Rscript ${IdentifyGISAIDGlobalVariants} \
	  ${dir}/H3N2-${Coding[$i]}-unpassaged-ref.aligned.coding \
	  ${dir}/H3N2-${Coding[$i]}-unpassaged.summary ${Coding[$i]} \
	  ${dir}/H3N2-${Coding[$i]}-unpassaged-ref.exclusions
  fi
  
  # Concatenate global variant data.
  if [ ! -f ${dir}/H3N2-unpassaged.summary.gz ]; then
    echo "Concatenate global variants."
    echo ${Genes[$i]}
	cat <( tail -n +2 ${dir}/H3N2-${Coding[$i]}-unpassaged.summary ) \
	  >> ${dir}/H3N2-unpassaged.summary.temp
  fi
  
  # Concatenate the sequence distance data for quality checks.
  if [ ! -f ${dir}/H3N2-unpassaged.distances.gz ]; then
    echo "Concatenate global distances."
    echo ${Genes[$i]}
	cat <( tail -n +2 ${dir}/H3N2-${Coding[$i]}-unpassaged-ref.distances.annotated ) \
	  >> ${dir}/H3N2-unpassaged.distances.temp
  fi
  
done

if [ ! -f ${dir}/H3N2-unpassaged.summary.gz ]; then
  # Compress the concatenated global variant data.
  cat <( head -n 1 ${dir}/H3N2-${Coding[$i]}-unpassaged.summary ) \
    ${dir}/H3N2-unpassaged.summary.temp > ${dir}/H3N2-unpassaged.summary
  rm -f ${dir}/H3N2-unpassaged.summary.temp
  gzip -f ${dir}/H3N2-unpassaged.summary
  # Copy the summary to a tracked directory.
  cp ${dir}/H3N2-unpassaged.summary.gz \
    analysis/${analysis}/out/H3N2-unpassaged.summary.gz
fi

if [ ! -f ${dir}/H3N2-unpassaged.distances.gz ]; then
  # Compress the concatenated sequence distance data.
  cat <( head -n 1 ${dir}/H3N2-${Coding[$i]}-unpassaged-ref.distances.annotated ) \
    ${dir}/H3N2-unpassaged.distances.temp > ${dir}/H3N2-unpassaged.distances
  rm -f ${dir}/H3N2-unpassaged.distances.temp
  gzip -f ${dir}/H3N2-unpassaged.distances
  # Copy the summary to a tracked directory.
  cp ${dir}/H3N2-unpassaged.distances.gz \
    analysis/${analysis}/out/H3N2-unpassaged.distances.gz  
fi


##################################
# Analysis of H3N2 phylogenies.
##################################

for i in $(seq 0 $((${#Genes[@]}-1)) )
do

  # Split the sequences by season from 2016 onwards.
  if [ ! -f ${dir}/H3N2-${Coding[$i]}-unpassaged-ref.aligned.coding.2016 ]; then
    echo "Split sequences by season."
	echo ${Genes[$i]}
    Rscript ${SplitSequencesBySeason} \
	  ${dir}/H3N2-${Coding[$i]}-unpassaged-ref.aligned.coding 2016 \
	  ${dir}/H3N2-${Coding[$i]}-unpassaged-ref.exclusions
  fi

done

# Submit batch job to build RAxML trees.
if [ ! -f ${dir}/tree/trees.txt ]; then
  # Create a samplesheet for the batch RAxML submission script.
  for file in ${dir}/H3N2-*.coding.201*
  do
    gene=${file##*H3N2-}
    gene=${gene%%-unpassaged*}
    year=${file##*.}
    echo ${file} H3N2-${gene}-${year}
  done > ${dir}/tree/trees.txt
  
  echo "Submit jobs, build trees."
  numsamples="$(wc -l ${dir}/tree/trees.txt | cut -f1 -d' ')"
  qsub -cwd -N BuildTree -l h_rt=96:00:00 \
    -t 1-${numsamples} -tc 250 -l m_mem_free=48G \
    -o nobackup/sge -e nobackup/sge \
    ${BatchBuildTreeRAxML} \
    ${dir}/tree ${dir}/tree/trees.txt 1
fi

numsamples="$(wc -l ${dir}/tree/trees.txt | cut -f1 -d' ')"
# Root trees using Victoria/2011 sequence as outgroup.
if [ ! -f ${dir}/tree/RAxML_bestTree.H3N2-1-PB2-2016.rooted ] && \
  [ "$(ls ${dir}/tree/*bestTree* | wc -l)" -eq ${numsamples} ]; then
  echo "Root trees."
  for file in ${dir}/tree/RAxML_bestTree.H3N2-*-201*
  do
    # Extract the name of the outgroup, which is the first sequence.
    stem=${file##*.}
	outgroup="$( head -n 1 ${dir}/tree/${stem}.sequences.RAxML )"
	outgroup=${outgroup##*>}
	echo "Rooting ${stem}"
    python ${RootTree} ${file} ${outgroup} ${file}.rooted
  done
fi


# Perform ancestral reconstruction for RAxML trees.
if [ ! -f ${dir}/tree/ancestors.txt ] && \
  [ "$(ls ${dir}/tree/*bestTree*.rooted | wc -l)" -eq ${numsamples} ]; then
  
  # Set up run parameters for batch job submission.
  while read sequences name
  do
    mkdir -p ${dir}/tree/${name}
    echo ${dir}/tree/${name}.sequences.RAxML ${dir}/tree/RAxML_bestTree.${name}.rooted ${dir}/tree/${name}
  done < ${dir}/tree/trees.txt > ${dir}/tree/ancestors.txt
  
  echo "Submit jobs, reconstruct ancestors."
  qsub -cwd -N ReconstructAncestors -l h_rt=96:00:00 \
    -t 1-${numsamples} -tc 250 -l m_mem_free=32G \
    -o nobackup/sge -e nobackup/sge \
    ${BatchReconstructAncestors} \
    ${dir}/tree/ancestors.txt 1
fi

# Extract frequencies of mutations from the reconstructed trees.
if [ ! -f ${dir}/tree/mutations.txt ] && \
  [ "$(ls ${dir}/tree/*/annotated_tree.nexus | wc -l)" -eq ${numsamples} ]; then
  echo "Extract mutation frequencies from reconstructed trees."
  for directory in ${dir}/tree/*/
  do
	# Extract the name of the ougroup.
	stem=${directory##*tree/}
	stem=${stem%%/*}
	echo ${stem}
	outgroup="$( head -n 1 ${dir}/tree/${stem}.sequences.RAxML )"
	outgroup=${outgroup##*>}
	echo ${directory}annotated_tree.nexus \
	  ${directory}ancestral_sequences.fasta ${outgroup} ${directory}frequencies.mutations \
	  >> ${dir}/tree/mutations.txt
  done
  qsub -cwd -N ExtractMutationFrequencies -l h_rt=96:00:00 \
    -t 1-${numsamples} -tc 250 -l m_mem_free=16G \
    -o nobackup/sge -e nobackup/sge \
    ${BatchExtractMutationFrequencies} \
    ${dir}/tree/mutations.txt 1
fi

# Concatenate the mutations from the reconstructed trees.
if [ ! -f analysis/${analysis}/out/H3N2-unpassaged.mutations.gz ] && \
  [ "$(ls ${dir}/tree/*/frequencies.mutations | wc -l)" -eq ${numsamples} ]; then
  echo "Concatenate mutations"
  for directory in ${dir}/tree/*/
  do
    # Extract the name of the gene and year.
	stem=${directory##*tree/}
	stem=${stem%%/}
	echo ${stem}
	cat <( zless -S ${directory}/frequencies.mutations | sed "s/$/ ${stem}/g" | tail -n +2 ) \
	  >> ${dir}/tree/H3N2-unpassaged.mutations.temp
  done
  cat <( head -n 1 ${dir}/tree/H3N2-1-PB2-2016/frequencies.mutations | sed "s/$/ Sequences/g" ) \
    ${dir}/tree/H3N2-unpassaged.mutations.temp > ${dir}/tree/H3N2-unpassaged.mutations
  rm -f ${dir}/tree/H3N2-unpassaged.mutations.temp
  gzip -f ${dir}/tree/H3N2-unpassaged.mutations
  mv ${dir}/tree/H3N2-unpassaged.mutations.gz \
    analysis/${analysis}/out/H3N2-unpassaged.mutations.gz
fi

# Summarize the proportion of NS variants in the within- and between-host data.
Rscript ${CalculateNSProportions}
