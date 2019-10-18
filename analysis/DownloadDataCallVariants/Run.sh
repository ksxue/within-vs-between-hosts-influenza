#!/bin/bash

# Script is meant to be run from the top level of the Github repository.
# Script downloads raw sequencing files from the SRA for each project
# and renames them appropriately.
# Script then aligns and calls variants in each of the sequenced samples.

analysis="DownloadDataCallVariants"
mkdir -p nobackup/${analysis}

# Dinis et al. 2016, J. Virol.
# Raw sequencing data was obtained through personal communication
# with the manuscript's authors.
# Files were obtained in the form of a single FASTQ file per sample.
analysis/${analysis}/download-flu-Dinis.sh
analysis/${analysis}/align-call-variants-flu-Dinis.sh

# Debbink et al. 2017, PLoS Pathogens
# Raw sequencing data was obtained from the SRA.
analysis/${analysis}/download-flu-Debbink.sh
analysis/${analysis}/align-call-variants-flu-Debbink.sh

# McCrone et al. 2018, eLife
# Raw sequencing data was obtained from the SRA.
analysis/${analysis}/download-flu-McCrone.sh
analysis/${analysis}/align-call-variants-flu-McCrone.sh

# Xue et al. 2017, eLife
# Raw sequencing data was obtained from the SRA.
analysis/${analysis}/download-flu-Xue-chronic.sh
analysis/${analysis}/align-call-variants-flu-Xue-chronic.sh

