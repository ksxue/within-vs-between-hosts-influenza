#!/bin/bash

# Script is meant to be run from the top level of the Github repository.
# Script runs sub-scripts to conduct comparisons of evolutionary rates within and between hosts.

# Download sequencing data from the SRA.
# Align reads and parse BAM files to produce variant summaries.
analysis/DownloadDataCallVariants/Run.sh

# Call within-host variants in each sample.
# Identify all sites with coverage >400x and non-consensus frequency >0%
# for an overall analysis of within-host variants.
# Also produce a subset of variants with frequency >0.5%.
analysis/CallWithinHostVariants/Run.sh

# Estimate rates of evolution in acute infections.
# Perform preliminary analyses like parsing the sample metadata,
# identifying excluded samples, parsing antigenic site classifications,
# and identifying the number of available S, NS, and Stop sites.
# Calculate per-site divergence for each sample.
# Use this information and the sample DPIs to estimate rates of divergence.
analysis/CalculateAcuteRates/Run.sh

# Estimate rates of evolution in chronic infections 
# by calculating sample divergence and using linear regression.
analysis/CalculateChronicRates/Run.sh

# Estimate global rates of evolution
# by analyzing GISAID sequences collected between 1999 and 2017,
# calculating their S and NS distances from a reference sequence,
# and performing a linear regression of sequence distances vs collection time.
analysis/CalculateGlobalRates/Run.sh

# Analyze the distribution of within-host diversity
# at different variant-frequency thresholds.
analysis/CompareVariantCallingThresholds/Run.sh

# Calculate the proportion of nonsynonymous mutations
# in the global flu population using GISAID sequences 
# collected between 2015 and 2018.
# Align the sequences, generate a phylogeny, reconstruct ancestral states,
# and traverse the tree to calculate the number of descendants of each lineage.
analysis/InferGlobalVariants/Run.sh

# Plot figures analyzing within- and between-host variation.
Rscript analysis/figures-WvB/PlotFigures-WvB.R

