#!/bin/bash

# Script is meant to be run from the top level of the Github repository.
# Script runs sub-scripts to analyze within-host variation
# across different variant-frequency thresholds.

analysis="CompareVariantCallingThresholds"

# Analyze the number of variants at different variant-frequency thresholds.
# Compare the number of variants at the first, second, and third codon positions.
Rscript analysis/${analysis}/CompareVariantCallingThresholds.R
