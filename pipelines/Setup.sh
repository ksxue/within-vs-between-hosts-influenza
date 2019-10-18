#!/bin/bash

# This script is meant to be run from the top level of the Github repository.
# This script sets up some of the main analyses 
# by preparing indexed reference genomes, compiled scripts, and so on.


##################################################################
# Index reference genomes and compile scripts.
##################################################################

module load bowtie2/2.2.3
bowtie2-build reference/flu-H3N2/H3N2-Brisbane-2007.fasta reference/flu-H3N2/H3N2-Brisbane-2007
bowtie2-build reference/flu-H3N2/H3N2-Victoria-2011.fasta reference/flu-H3N2/H3N2-Victoria-2011
bowtie2-build reference/flu-pdmH1N1/pdmH1N1-California-2009.fasta reference/flu-pdmH1N1/pdmH1N1-California-2009
bowtie2-build reference/flu-seasonalH1N1/H1N1-Boston-2007.fasta reference/flu-seasonalH1N1/H1N1-Boston-2007

module load mpc/0.8.2
module load mpfr/3.1.0
module load gmp/5.0.2
module load gcc/4.9.1
g++ -O scripts/SummarizeBAM.cpp -o bin/SummarizeBAM-1.21
g++ -O scripts/AnnotateVariants.cpp -o bin/AnnotateVariants-1.1
