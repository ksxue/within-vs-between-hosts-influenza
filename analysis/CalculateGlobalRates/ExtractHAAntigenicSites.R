require(tidyverse)
require(Biostrings)
require(seqinr)
require(purrr)
require(generics)

source("analysis/CalculateAcuteRates/ImportAntigenicSites.R")

# This script takes a set of H3 HA sequences that have been downloaded
# from GISAID and aligned to a reference sequence.
# Using the set of (computationally numbered) HA sites annotated as
# antigenic sites defined by Wolf, it then splits the sequence alignment into two sets
# of antigenic and non-antigenic sites.

args = commandArgs(trailingOnly=TRUE)

# Verify that the correct number of arguments is given.
if(length(args)!=3){
  stop("These arguments must be supplied: input FASTA file, 
       output antigenic FASTA file, output non-antigenic FASTA file.", 
       call.=FALSE)
}

FASTA <- args[1]
OUTANTIGENIC <- args[2]
OUTNONANTIGENIC <- args[3]

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


# Extract the antigenic and non-antigenic regions of the gene.
# First, split the alignment into separate columns for each position again.
AlignedRef <- strsplit(as.character(((head(Data,1) %>% dplyr::select(Sequence))$Sequence)),"")[[1]]
Data <- Data %>% 
  separate(Sequence, into=as.character(seq(0,MaxSeqLength)), sep="|")
# Then, extract the list of nucleotide positions in the antigenic sites.
# Note that the current antigenic site list has codons in computational numbering.
AntigenicNt <- sort(c(WolfAntigenicSites$Codon*3-2, 
                 WolfAntigenicSites$Codon*3-1,
                 WolfAntigenicSites$Codon*3))
NonantigenicNt <- setdiff(seq(1:length(AlignedRef)*3), AntigenicNt)
# Select the positions corresponding to antigenic and non-antigenic sites.
DataAntigenic <- Data %>% dplyr::select(Header, as.character(AntigenicNt)) %>%
  unite("Sequence", as.character(AntigenicNt), sep="")
DataNonantigenic <- Data %>% dplyr::select(Header, as.character(NonantigenicNt)) %>%
  unite("Sequence", as.character(NonantigenicNt), sep="")

# The information contained in the FASTA headers downloaded from GISAID
# are listed in the fields below.
# The fields are delimited by the | character.
FASTAHeader <- c("IsolateName","DNAAccession","Type","Lineage","CollectionDate",
                 "Segment","PassageHistory","OriginatingLab","SubmittingLab")

# Modify the "segment" field of the GISAID header to reflect the
# antigenic and nonantigenic segments.
DataAntigenic <- DataAntigenic %>%
  separate(Header, into=FASTAHeader, sep=" \\| ") %>%
  mutate(Segment="HA-antigenic-Wolf") %>%
  unite(Header, FASTAHeader, sep=" | ")
DataNonantigenic <- DataNonantigenic %>%
  separate(Header, into=FASTAHeader, sep=" \\| ") %>%
  mutate(Segment="HA-nonantigenic-Wolf") %>%
  unite(Header, FASTAHeader, sep=" | ")

# Export remaining sequences.
Sequences <- DNAStringSet(DataAntigenic$Sequence)
names(Sequences) <- DataAntigenic$Header
writeXStringSet(Sequences, OUTANTIGENIC)

Sequences <- DNAStringSet(DataNonantigenic$Sequence)
names(Sequences) <- DataNonantigenic$Header
writeXStringSet(Sequences, OUTNONANTIGENIC)
