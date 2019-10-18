#!/bin/bash

# Script is meant to be run from the top level of the Github repository.
# Script runs sub-scripts to calculate rates of divergence in chronic infections.

analysis="CalculateChronicRates"

# Estimate rates of evolution in chronic infections 
# by calculating sample divergence and using linear regression.
Rscript analysis/${analysis}/CalculateChronicRates.R
