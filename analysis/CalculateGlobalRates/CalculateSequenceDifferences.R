require(tidyverse)
require(Biostrings)
require(seqinr)
require(purrr)

# This script takes in a set of sequences in FASTA format that have been
# downloaded from GISAID and aligned to a reference sequence, along with
# a reference sequence from which all subsequent distances are calculated.
# It removes sequences that have been passaged from subsequent analyses,
# and for each sequence, it outputs the number of S and NS differences
# between each sequence and the reference sequence.
# This information is output as a tab-delimited file,
# with one row for each sequence, and the dates and numbers of differences.

args = commandArgs(trailingOnly=TRUE)

# Verify that the correct number of arguments is given.
if(length(args)!=4){
  stop("These arguments must be supplied: input FASTA file, reference FASTA file,
       reference sequence year of origin,
       output summary file.", 
       call.=FALSE)
}

FASTA <- args[1]
REF <- args[2]
REFYEAR <- as.integer(args[3])
OutputFile <- args[4]

# Import the functions for calculating sequence differences.
source("analysis/CalculateGlobalRates/CountSequenceDifferenceFunctions.R")


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
Data <- Data %>% dplyr::select(FASTAHeader, Year, as.character(CodingPositions))
Data <- Data %>% unite("Sequence", as.character(CodingPositions), sep="")

# Retain only sequences that were collected during or after
# the year of the reference sequence.
# This will remove the initial sequence to which others were aligned,
# so wait to perform this step until after the alignment has been truncated.
Data <- Data %>% filter(Year>=REFYEAR)

# Read in the reference sequence file.
fastaFile <- readDNAStringSet(REF)
Ref <- data.frame(names(fastaFile), paste(fastaFile)) %>%
  `colnames<-`(c("Header","Sequence"))
rm(fastaFile)
RefSeq <- as.character(Ref$Sequence[1])
# If the sequence being analyzed is NS, then truncate
# the reference to include only the first 693 bases
# to account for a length polymorphism between early and late
# NS sequences.
if("NS" %in% Data$Segment){
  RefSeq <- substr(RefSeq,1,693)
}

# Calculate the number of S and NS differences between each sequence
# and the reference.
Distances <- sapply(Data$Sequence, function(x) {
  CountSeqDifferences(RefSeq,x)})
Data$S <- Distances[1,]
Data$NS <- Distances[2,]
 
# Export data fields containing collection year,
# S substitutions from reference, and NS substitutions from reference
# for each sequence.
write.table(Data %>% dplyr::select(Segment, Year, DNAAccession, CollectionDate, S, NS) %>%
              filter(!is.na(Year)) %>% mutate(RefYear=REFYEAR), OutputFile,
            quote=FALSE, row.names=FALSE)
