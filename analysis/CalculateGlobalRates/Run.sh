#!/bin/bash

# Script is meant to be run from the top level of the Github repository.

module load mafft/7.407

analysis="CalculateGlobalRates"
CalculateDifferences="analysis/${analysis}/CalculateSequenceDifferences.R"
SubsampleByYear="analysis/${analysis}/SubsampleByYear.R"
ExtractHAAntigenicSites="analysis/${analysis}/ExtractHAAntigenicSites.R"
ExtractNASurfaceSites="analysis/${analysis}/ExtractNASurfaceSites.R"
CalculateGlobalRates="analysis/${analysis}/CalculateGlobalRates.R"
CombineRates="analysis/${analysis}/CombineRates.R"

mkdir -p nobackup/${analysis}
dir="nobackup/${analysis}"
mkdir -p analysis/${analysis}/out


##################################
# Analysis of H3N2 sequences.
##################################

# For each influenza gene, concatenate all downloaded sequences from 1968 to 2017
# into a single file and subsample to a maximum of 50 sequences per year.
# Calculate the distance between each sequence and the 2007 reference sequence.
seqdir="data/GISAID/H3N2"
Genes=( "1-PB2" "2-PB1" "3-PA" "4-HA" "5-NP" "6-NA" "7-M" "8-NS" )
Coding=( "1-PB2" "2-PB1" "3-PA" "4-HA" "5-NP" "6-NA" "7-M1" "8-NS1" )
NumGenes=${#Genes[@]}
for i in $(seq 0 $((${#Genes[@]}-1)) )
do
  # Concatenate all downloaded sequences into a single file.
  if [ ! -f ${dir}/H3N2-${Genes[$i]}.fasta ]; then
    echo "Concatenate sequences, ${Genes[$i]}"
    cat ${seqdir}/H3N2-${Genes[$i]}-*.fasta > \
      ${dir}/H3N2-${Genes[$i]}.fasta
  fi
  # Subsample the sequences to a maximum of 50 per year.
  if [ ! -f ${dir}/H3N2-${Genes[$i]}-subsampled.fasta ]; then
    echo "Subsample sequences, ${Genes[$i]}"
    Rscript ${SubsampleByYear} \
	  ${dir}/H3N2-${Genes[$i]}.fasta ${dir}/H3N2-${Genes[$i]}-subsampled.fasta 50
  fi
  # Concatenate the reference sequence for a gene with the subsampled sequences.
  # Then, use MAFFT to build a multiple sequence alignment.
  if [ ! -f ${dir}/H3N2-${Coding[$i]}-subsampled-ref.aligned ]; then
    echo "Align sequences, ${Genes[$i]}"
	cat reference/flu-H3N2/H3N2-Victoria-2011-${Coding[$i]}-coding.fasta \
	  ${dir}/H3N2-${Genes[$i]}-subsampled.fasta \
	  > ${dir}/H3N2-${Genes[$i]}-subsampled-ref.fasta
    mafft ${dir}/H3N2-${Genes[$i]}-subsampled-ref.fasta \
	  > ${dir}/H3N2-${Coding[$i]}-subsampled-ref.aligned
	rm -f ${dir}/H3N2-${Genes[$i]}-subsampled-ref.fasta
  fi
  # Extract sets of sequences corresponding to HA antigenic and non-antigenic sites.
  # Use the Wolf definition of antigenic sites.
  if [ ! -f ${dir}/H3N2-4-HA-nonantigenic-Wolf-subsampled-ref.aligned ]; then
    echo "Extract antigenic and non-antigenic sites."
	Rscript ${ExtractHAAntigenicSites} ${dir}/H3N2-4-HA-subsampled-ref.aligned \
	  ${dir}/H3N2-4-HA-antigenic-Wolf-subsampled-ref.aligned \
	  ${dir}/H3N2-4-HA-nonantigenic-Wolf-subsampled-ref.aligned
  fi
  
  # Extract sets of sequences corresponding to NA surface and non-surface sites.
  # Use the Bhatt definition of surface sites.
  if [ ! -f ${dir}/H3N2-6-NA-nonsurface-Bhatt-subsampled-ref.aligned ]; then
    echo "Extract surface and non-surface sites."
	Rscript ${ExtractNASurfaceSites} ${dir}/H3N2-6-NA-subsampled-ref.aligned \
	  ${dir}/H3N2-6-NA-surface-Bhatt-subsampled-ref.aligned \
	  ${dir}/H3N2-6-NA-nonsurface-Bhatt-subsampled-ref.aligned
  fi
done

# Expand the list of genes to include antigenic and non-antigenic sites.
# Also include surface and non-surface sites.
# Then calculate the distance of each sequence from several reference sequences.
Genes=( "1-PB2" "2-PB1" "3-PA" "4-HA" "4-HA-antigenic-Wolf" "4-HA-nonantigenic-Wolf" "5-NP" "6-NA" "6-NA-surface-Bhatt" "6-NA-nonsurface-Bhatt" "7-M" "8-NS" )
Coding=( "1-PB2" "2-PB1" "3-PA" "4-HA" "4-HA-antigenic-Wolf" "4-HA-nonantigenic-Wolf" "5-NP" "6-NA" "6-NA-surface-Bhatt" "6-NA-nonsurface-Bhatt" "7-M1" "8-NS1" )
NumGenes=${#Genes[@]}
for i in $(seq 0 $((${#Genes[@]}-1)) )
do
  # Calculate the distance of each sequence
  # from a 1968 reference sequence.
  if [ ! -f ${dir}/H3N2-${Coding[i]}-subsampled-distances-1968.data ]; then
    echo "Calculate sequence differences, ${Genes[$i]}"
	Rscript ${CalculateDifferences} \
	  ${dir}/H3N2-${Coding[$i]}-subsampled-ref.aligned \
	  reference/flu-H3N2/H3N2-Aichi-1968-${Coding[$i]}-coding.fasta \
	  1968 ${dir}/H3N2-${Coding[i]}-subsampled-distances-1968.data
  fi
  # Calculate the distance of each sequence
  # from a 1999 reference sequence.
  if [ ! -f ${dir}/H3N2-${Coding[i]}-subsampled-distances-1999.data ]; then
    echo "Calculate sequence differences, ${Genes[$i]}"
	Rscript ${CalculateDifferences} \
	  ${dir}/H3N2-${Coding[$i]}-subsampled-ref.aligned \
	  reference/flu-H3N2/H3N2-Moscow-1999-${Coding[$i]}-coding.fasta \
	  1999 ${dir}/H3N2-${Coding[i]}-subsampled-distances-1999.data
  fi
  # Calculate the distance of each sequence
  # from a 2007 reference sequence.
  if [ ! -f ${dir}/H3N2-${Coding[i]}-subsampled-distances-2007.data ]; then
    echo "Calculate sequence differences, ${Genes[$i]}"
	Rscript ${CalculateDifferences} \
	  ${dir}/H3N2-${Coding[$i]}-subsampled-ref.aligned \
	  reference/flu-H3N2/H3N2-Brisbane-2007-${Coding[$i]}-coding.fasta \
	  2007 ${dir}/H3N2-${Coding[i]}-subsampled-distances-2007.data
  fi
done

# Perform linear regressions to estimate the rates of evolution
# of each gene in 10-year intervals from the 1968 reference sequence.
# Trim the first 3 years of sequences following the pandemic
# to reduce the influence of reassortment.
if [ ! -f analysis/${analysis}/out/GlobalRates-H3N2-1968.data ]; then
  echo "Estimate global evolutionary rates, 1968"
  Rscript ${CalculateGlobalRates} ${dir} 1968 10 3 \
    analysis/${analysis}/out/
fi
# Perform linear regressions to estimate the rates of evolution
# of each gene in 20-year intervals from the 1999 reference sequence.
# This interval is necessary to estimate evolutionary rates of
# more slowly-evolving genes.
if [ ! -f analysis/${analysis}/out/GlobalRates-H3N2-1999.data ]; then
  echo "Estimate global evolutionary rates, 1999"
  Rscript ${CalculateGlobalRates} ${dir} 1999 20 0 \
    analysis/${analysis}/out/
fi
# Perform linear regressions to estimate the rates of evolution
# of each gene in 10-year intervals from the 2007 reference sequence.
# This interval is necessary to estimate the evolutionary rates of
# more rapidly evolving genes without saturating available sites.
# Trim the 2 years following the reference sequence.
if [ ! -f analysis/${analysis}/out/GlobalRates-H3N2-2007.data ]; then
  echo "Estimate global evolutionary rates, 2007"
  Rscript ${CalculateGlobalRates} ${dir} 2007 10 0 \
    analysis/${analysis}/out/
fi

# Aggregate rates of evolution calculated through different methods.
echo "Combine rates of evolution calculated through different methods."
Rscript ${CombineRates}
