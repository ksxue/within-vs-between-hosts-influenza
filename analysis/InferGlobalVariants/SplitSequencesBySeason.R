suppressMessages(require(tidyverse))
suppressMessages(require(Biostrings))
suppressMessages(require(seqinr))
suppressMessages(require(foreach))

args = commandArgs(trailingOnly=TRUE)

# Verify that exactly three arguments are given.
if(length(args)!=3){
  stop("These arguments must be supplied: input FASTA file, first year, and list of sequence exclusions.", 
       call.=FALSE)
}

FASTA <- args[1]
YEAR <- as.integer(args[2])
Exclusions <- args[3]


# Read in file format information.
source("analysis/FileFormats.R")

# Read in FASTA file.
fastaFile <- readDNAStringSet(FASTA)
Data <- data.frame(names(fastaFile), paste(fastaFile)) %>%
  `colnames<-`(c("Header","Sequence"))
rm(fastaFile)

# Identify the reference sequence.
Ref <- Data[1,]

# Read in list of exclusions and omit outlier sequences.
SequenceExclusions <- read.table(Exclusions, header=TRUE, stringsAsFactors = FALSE)
Data <- Data %>% filter(!(Header %in% SequenceExclusions$Header))

# Parse the FASTA header.
Data <- Data %>%
  separate(Header, into=GISAIDFASTAHeader, sep=" \\| ") %>%
  mutate(Year=round(as.numeric(CollectionDate)))
# Exclude the reference sequence.
Data <- Data %>% filter(!is.na(Year))

# Export sequences by year, starting from the indicated year.
foreach(year=seq(YEAR,max(Data$Year))) %do% {
  # Select only sequences within the specified year.
  DataYear <- Data %>% unite(Header, GISAIDFASTAHeader, sep=" | ") %>%
    filter(Year==year) %>% dplyr::select(-Year)
  # Add on the reference sequence.
  DataYear <- rbind(Ref, DataYear)
  # Export the sequences.
  Sequences <- DNAStringSet(DataYear$Sequence)
  names(Sequences) <- DataYear$Header
  writeXStringSet(Sequences, 
                  paste0(FASTA,".",as.character(year)))  
}

