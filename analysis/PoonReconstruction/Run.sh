#!/bin/bash

# Script is meant to be run from the top level of the Github repository.
# Script analyzes within-host variants in the Poon et al. raw data
# with the goal of reproducing Figure 2 from the published analyses.

# Input parameters and directory paths.
project="flu-Poon"
analysis="PoonReconstruction"
outdir="nobackup/${analysis}/"
mkdir -p ${outdir}
mkdir -p analysis/${analysis}/out

# Output a list of variable sites for each flu subtype,
# retaining only sites that are predicted as variable in both sequencing replicates.
# Note that this excludes some samples without sequencing replicates.
Rscript analysis/${analysis}/CallReplicateVariants-flu-Poon.R
Rscript analysis/${analysis}/ReconstructFigure2.R
