require(tidyverse)
require(Biostrings)
require(seqinr)
require(purrr)

# This script takes in a set of sequences in FASTA format that have been
# downloaded from GISAID.
# It parses the GISAID-generated FASTA header for information on collection year.
# It then randomly subsamples the sequences per year to up to
# the specified number of sequences.
# It outputs the sequences in FASTA format for further analysis.

args = commandArgs(trailingOnly=TRUE)

# Verify that the correct number of arguments is given.
if(length(args)!=3){
  stop("These arguments must be supplied: input FASTA file,
       output FASTA file, and maximum number of sequences per year.", 
       call.=FALSE)
}

FASTA <- args[1]
OutputFile <- args[2]
MaxNumSeqYear <- as.integer(args[3])

set.seed(0)

# The information contained in the FASTA headers downloaded from GISAID
# are listed in the fields below.
# The fields are delimited by the | character.
FASTAHeader <- c("IsolateName","DNAAccession","Type","Lineage","CollectionDate",
                 "Segment","PassageHistory","OriginatingLab","SubmittingLab")

# Read in FASTA file.
fastaFile <- readDNAStringSet(FASTA)
Data <- data.frame(names(fastaFile), paste(fastaFile)) %>%
  `colnames<-`(c("Header","Sequence"))
rm(fastaFile)

# Parse the FASTA header.
# Also parse the collection date to obtain the year value.
# Retain sequences whose month and day of collection is not known,
# grouping only by year instead.
Data <- Data %>%
  separate(Header, into=FASTAHeader, sep=" \\| ") %>%
  separate(CollectionDate, into=c("CollectionDate"), sep=" ", extra="drop") %>%
  mutate(Year=round(as.numeric(CollectionDate)))
# Convert the collection date to a numeric value.
Data$CollectionDate <- as.numeric(as.character(Data$CollectionDate))

# Subsample the data to retain up to MaxNumSeqYear unpassaged sequences per year.
# Retain the number of sequences that are available if less than the maximum.
# Remove sequences whose year is uninterpretable, i.e. NA.
Data <- Data %>% filter(!is.na(Year)) %>% group_by(Year) %>% 
  mutate(SampleSize=min(n(), MaxNumSeqYear)) %>%
  ungroup() %>% group_by(Year, SampleSize) %>%
  nest() %>%
  mutate(samp = map2(data, SampleSize, sample_n)) %>%
  select(Year, samp) %>%
  unnest()

# Export remaining sequences.
Data <- Data %>% unite(Header, FASTAHeader, sep=" | ")
Sequences <- DNAStringSet(Data$Sequence)
names(Sequences) <- Data$Header
writeXStringSet(Sequences, OutputFile)
