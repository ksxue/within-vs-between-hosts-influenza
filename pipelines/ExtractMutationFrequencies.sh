#!/bin/bash

# Script is meant to be run from the top level of the Github repository.

# Input parameters.
tree="$1" # Give the file containing the annotated tree.
seqs="$2" # Give the ancestral sequences that corresponds to it.
outgroup="$3" # Give the desired outgroup.
out="$4" # Give the path to the output file.
clean="$5" # If 0, then the script will run in its entirety, overwriting existing output.


if [ ! -f ${out} ] || [ ${clean} -eq "0" ];
then

  module load python/3.6.4
  python3 scripts/ExtractMutationFrequencies.py ${tree} ${seqs} ${outgroup} ${out}
  
fi

