require(tidyverse)
require(Biostrings)
require(seqinr)
require(purrr)

# This script takes in a set of sequences in FASTA format that have been
# downloaded from GISAID and aligned to a reference sequence,
# and it prunes the alignment to retain only the coding sequences.
# It then outputs these coding sequences in the form of a FASTA file.

args = commandArgs(trailingOnly=TRUE)

# Verify that the correct number of arguments is given.
if(length(args)!=2){
  stop("These arguments must be supplied: input FASTA file, output FASTA file.", 
       call.=FALSE)
}

FASTA <- args[1]
OutputFile <- args[2]

# Read in FASTA file.
fastaFile <- readDNAStringSet(FASTA)
Data <- data.frame(names(fastaFile), paste(fastaFile)) %>%
  `colnames<-`(c("Header","Sequence"))
rm(fastaFile)

# Truncate alignment to retain only coding region of the gene.
# First, extract the aligned reference sequence
# and determine which alignment positions do not contain gaps in the reference,
# i.e. are coding.
# Split the alignment into separate columns for each position.
# Remove any whitespace that may have been added while loading the file.
# Then, retain only positions that have been determined to be coding.
# Merge the remaining alignment positions back into a coding sequence.
AlignedRef <- strsplit(as.character(((head(Data,1) %>% dplyr::select(Sequence))$Sequence)),"")[[1]]
MaxSeqLength <- max(sapply(as.character(Data$Sequence), nchar))
Data <- Data %>% 
  separate(Sequence, into=as.character(seq(0,MaxSeqLength)), sep="|")
CodingPositions <- which(AlignedRef!="-")
NoncodingPositions <- which(AlignedRef=="-")
Data <- Data %>% dplyr::select(Header, as.character(CodingPositions))
Data <- Data %>% unite("Sequence", as.character(CodingPositions), sep="")

# Remove any exact duplicates of name and sequence.
Data <- Data %>% distinct()
# Remove any sequences that have the same name but not the same sequence.
Data <- Data %>% 
  group_by(Header) %>% filter(n()==1)

# Export remaining sequences.
Sequences <- DNAStringSet(Data$Sequence)
names(Sequences) <- Data$Header
writeXStringSet(Sequences, OutputFile)
