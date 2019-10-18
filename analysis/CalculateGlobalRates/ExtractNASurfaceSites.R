require(tidyverse)
require(Biostrings)
require(seqinr)
require(purrr)
require(generics)

source("analysis/CalculateAcuteRates/ImportAntigenicSites.R")

# This script takes a set of N2 NA sequences that have been downloaded
# from GISAID and aligned to a reference sequence.
# Using the set of (computationally numbered)NA sites annotated as
# surface-exposed sites defined by Bhatt, it then splits the sequence alignment into two sets
# of surface-exposed and non-surface-exposed sites.

args = commandArgs(trailingOnly=TRUE)

# Verify that the correct number of arguments is given.
if(length(args)!=3){
  stop("These arguments must be supplied: input FASTA file, 
       output antigenic FASTA file, output non-antigenic FASTA file.", 
       call.=FALSE)
}

FASTA <- args[1]
OUTSURFACE <- args[2]
OUTNONSURFACE <- args[3]

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


# Extract the surface and non-surface regions of the gene.
# First, split the alignment into separate columns for each position again.
AlignedRef <- strsplit(as.character(((head(Data,1) %>% dplyr::select(Sequence))$Sequence)),"")[[1]]
Data <- Data %>% 
  separate(Sequence, into=as.character(seq(0,MaxSeqLength)), sep="|")
# Then, extract the list of nucleotide positions in the surface-exposed sites.
# Note that the current antigenic site list has codons in computational numbering.
SurfaceNt <- sort(c(BhattSurfaceExposedSites$Codon*3-2, 
                    BhattSurfaceExposedSites$Codon*3-1,
                    BhattSurfaceExposedSites$Codon*3))
NonsurfaceNt <- setdiff(seq(1:length(AlignedRef)*3), SurfaceNt)
# Select the positions corresponding to surface and non-surface sites.
DataSurface <- Data %>% dplyr::select(Header, as.character(SurfaceNt)) %>%
  unite("Sequence", as.character(SurfaceNt), sep="")
DataNonsurface <- Data %>% dplyr::select(Header, as.character(NonsurfaceNt)) %>%
  unite("Sequence", as.character(NonsurfaceNt), sep="")

# The information contained in the FASTA headers downloaded from GISAID
# are listed in the fields below.
# The fields are delimited by the | character.
FASTAHeader <- c("IsolateName","DNAAccession","Type","Lineage","CollectionDate",
                 "Segment","PassageHistory","OriginatingLab","SubmittingLab")

# Modify the "segment" field of the GISAID header to reflect the
# surface and nonsurface segments.
DataSurface <- DataSurface %>%
  separate(Header, into=FASTAHeader, sep=" \\| ") %>%
  mutate(Segment="NA-surface-Bhatt") %>%
  unite(Header, FASTAHeader, sep=" | ")
DataNonsurface <- DataNonsurface %>%
  separate(Header, into=FASTAHeader, sep=" \\| ") %>%
  mutate(Segment="NA-nonsurface-Bhatt") %>%
  unite(Header, FASTAHeader, sep=" | ")

# Export remaining sequences.
Sequences <- DNAStringSet(DataSurface$Sequence)
names(Sequences) <- DataSurface$Header
writeXStringSet(Sequences, OUTSURFACE)

Sequences <- DNAStringSet(DataNonsurface$Sequence)
names(Sequences) <- DataNonsurface$Header
writeXStringSet(Sequences, OUTNONSURFACE)
