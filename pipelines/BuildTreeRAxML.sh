#!/bin/bash

# Script is meant to be run from the top level of the Github repository.

# Input parameters.
sequences="$1" # Give the file of aligned sequences.
name="$2" # Give the desired sequence name.
dir="$3" # Give the desired output directory.
clean="$4" # If 0, then the script will run in its entirety, overwriting existing output.

# Check that the output directory exists.
mkdir -p ${dir}

if [ ! -f ${dir}/RAxML_bestTree.${name} ] || [ ${clean} -eq "0" ];
then

  # Format sequences for use with RAxML.
  # In particular, remove characters in the sequence names that may cause problems.
  less -S ${sequences} |  tr ' ' '_' | tr ":,();[]\'" "_" \
    > ${dir}/${name}.sequences.RAxML
  
  # Build a tree from the sequences using RAxML.
  raxmlHPC-SSE3 -s ${dir}/${name}.sequences.RAxML \
    -n ${name} -m GTRCAT -p 1
  
  # Once RAxML has finished running, move the outputs
  # to the desired output directory.
  mv -f RAxML*${name} ${dir}  
fi

