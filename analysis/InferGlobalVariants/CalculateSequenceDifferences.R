require(tidyverse)
require(Biostrings)
require(seqinr)
require(purrr)

# This script takes in a set of coding sequences in FASTA format that have been
# downloaded from GISAID and aligned to a reference coding sequence,
# and it calculates the S and NS distances between that each sequence and the reference.

args = commandArgs(trailingOnly=TRUE)

# Verify that the correct number of arguments is given.
if(length(args)!=2){
  stop("These arguments must be supplied: input FASTA file, output summary file.", 
       call.=FALSE)
}

FASTA <- args[1]
OutputFile <- args[2]

# Import functions for calculating sequence differences.
source("analysis/CalculateGlobalRates/CountSequenceDifferenceFunctions.R")

# Read in FASTA file.
fastaFile <- readDNAStringSet(FASTA)
Data <- data.frame(names(fastaFile), paste(fastaFile)) %>%
  `colnames<-`(c("Header","Sequence"))
rm(fastaFile)

# Identify the reference sequence.
RefSeq <- as.character(Data[1,2])

# Calculate the number of S and NS differences between each sequence
# and the reference.
Data$Sequence <- as.character(Data$Sequence)
Distances <- sapply(Data$Sequence, function(x) {
  CountSeqDifferences(RefSeq,x)})
Data$S <- Distances[1,]
Data$NS <- Distances[2,]

# Export table of sequence distances from the reference.
write.table(Data %>% dplyr::select(-Sequence), OutputFile,
            row.names=FALSE, quote=TRUE)