require(tidyverse)
require(Biostrings)
require(seqinr)

args = commandArgs(trailingOnly=TRUE)

# Verify that exactly two arguments are given.
if(length(args)!=2){
  stop("Two arguments must be supplied: input FASTA file, output summary file.", 
       call.=FALSE)
}

FASTA <- args[1]
OutputFile <- args[2]

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
Data <- Data %>%
  separate(Header, into=FASTAHeader, sep=" \\| ")

# Retain only sequences that have not been passaged.
# Discard the reference sequence at this point.
Data <- Data %>%
  filter(PassageHistory %in% c("P0 ","original ","Clinical Specimen ",
                               "passage details:Original ", "passage details: Original ",
                               "Original Specimen ", "clinical specimen ", "Clinical specimen ",
                               "Original specimen ", "Original ", "direct ", "Direct ",
                               "original specimen ", "passage details: direct clinical sample ",
                               "isolated directly from host; no passage ",
                               "passage details: original specimen ",
                               "passage details: primary clinical sample ",
                               "passage details: no passage ",
                               "Clinical sample ", "Original Samlpe ", "Original sample ",
                               "clinical sample ", "CLINICAL SPECIMEN ", "Clincal Specimen ",
                               "passage details: clinical sample ", "PCR on original specimen ",
                               "p0 ", "P0 (H3N2) ", "Original specimen uncultured in VTM "))

# Export remaining sequences.
Data <- Data %>% unite(Header, FASTAHeader, sep=" | ")
Sequences <- DNAStringSet(Data$Sequence)
names(Sequences) <- Data$Header
writeXStringSet(Sequences, OutputFile)
