#!/bin/bash

# Script is meant to be run from the top level of the Github repository.
# Script runs sub-scripts to calculate rates of divergence in acute infections.

analysis="CalculateAcuteRates"

# Parse sample metadata to organize information about DPI, viral load,
# and other sample characteristics.
Rscript analysis/${analysis}/ParseSampleMetadata.R
# Parse lists of samples to identify samples that should be excluded
# from subsequent analyses.
Rscript analysis/${analysis}/IdentifySampleExclusions.R
# Analyze a reference genome and calculate the number of available
# sites for S, NS, and stop mutations.
Rscript analysis/${analysis}/CalculateAvailableSites.R
# Analyze a reference genome and calculate the lengths of each coding sequence.
Rscript analysis/${analysis}/CalculateCodingSequenceLengths.R
# Parse lists of antigenic sites in H3 HA.
Rscript analysis/${analysis}/ImportAntigenicSites.R

# Import variant data and calculate the viral divergence for each sample.
Rscript analysis/${analysis}/CalculateAcuteSampleDivergence.R
# Use the metadata generated previously to convert viral divergence per sample
# into average rates of divergence for each viral gene.
# Also estimate rates through linear regression.
Rscript analysis/${analysis}/CalculateAcuteRates.R
