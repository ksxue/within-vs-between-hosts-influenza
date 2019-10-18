require(tidyverse)
require(Biostrings)
require(stringr)
require(foreach)
require(seqinr)

args = commandArgs(trailingOnly=TRUE)

# Verify that exactly three arguments are given.
if(length(args)!=4){
  stop("These arguments must be supplied: input FASTA file, output summary file,
        gene name, and list of sequence exclusions.", 
       call.=FALSE)
}

FASTA <- args[1]
OutputFile <- args[2]
Gene <- args[3]
Exclusions <- args[4]

# Read in file format information.
source("analysis/FileFormats.R")

# Read in FASTA file.
fastaFile <- readDNAStringSet(FASTA)
Data <- data.frame(names(fastaFile), paste(fastaFile)) %>%
  `colnames<-`(c("Header","Sequence"))

# Read in list of exclusions and omit outlier sequences.
SequenceExclusions <- read.table(Exclusions, header=TRUE, stringsAsFactors = FALSE)
Data <- Data %>% filter(!(Header %in% SequenceExclusions$Header))

# Parse the FASTA header.
Data <- Data %>%
  separate(Header, into=GISAIDFASTAHeader, sep=" \\| ") %>%
  mutate(Year=round(as.numeric(CollectionDate)))
# Exclude the reference sequence.
Data <- Data %>% filter(!is.na(Year))

# Summarize base counts at each position in the remaining sequences.
# First, exclude metadata and separate sequences into codons.
MaxCodingSeqLength <- max(sapply(as.character(Data$Sequence), nchar))
Sequences <- Data %>%
  dplyr::select(Year, Sequence)
Sequences <- Sequences %>%
  separate(Sequence, into=as.character(seq(1,MaxCodingSeqLength/3)), sep=seq(3,MaxCodingSeqLength-1,3))
# Convert data to tidy format.
Sequences <- Sequences %>% gather(Position, Codon, -Year)
Sequences$Position <- as.numeric(Sequences$Position)
# Summarize base frequencies for each year.
Sequences <- Sequences %>% group_by(Year, Position, Codon) %>%
  summarize(Count=n())
# Translate each recorded codon into an amino acid using the seqinr package.
Sequences$AA <- sapply(Sequences$Codon, function(x) seqinr::translate(s2c(x)))
# Remove codons that contain N bases, which are translated as X.
Sequences <- Sequences %>% filter(AA!="X")
# Determine the consensus (most frequent) amino acid at each position
# in each year.
Sequences <- Sequences %>% group_by(Year, Position) %>%
  mutate(ConsensusCodon=Codon[base::which.max(Count)],
         ConsensusAA=AA[base::which.max(Count)],
         Freq=Count/sum(Count))

# Export codon frequency tallies.
write.table(Sequences %>% mutate(Gene=Gene), OutputFile,
            quote=FALSE, row.names=FALSE)
