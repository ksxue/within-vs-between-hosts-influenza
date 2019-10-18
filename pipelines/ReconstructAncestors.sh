#!/bin/bash

# Script is meant to be run from the top level of the Github repository.

# Input parameters.
alignment="$1" # Give the file of aligned sequences.
tree="$2" # Give the tree that corresponds to it.
dir="$3" # Give the desired output directory.
clean="$4" # If 0, then the script will run in its entirety, overwriting existing output.

# Check that the output directory exists.
mkdir -p ${dir}

if [ ! -f ${dir}/annotated_tree.nexus ] || [ ${clean} -eq "0" ];
then

  treetime ancestral --aln ${alignment} --tree ${tree} --outdir ${dir}
  
fi

