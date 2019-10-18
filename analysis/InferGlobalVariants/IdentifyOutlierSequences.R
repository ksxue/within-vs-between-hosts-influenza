require(tidyverse)


# This script takes in a table of sequence names plus S and NS distances.
# It parses the sequence names (from GISAID) to identify the year in which they were collected,
# and it identifies and excludes sequences with abnormally high distances from the reference
# relative to the other sequences collected in that year.

args = commandArgs(trailingOnly=TRUE)

# Verify that the correct number of arguments is given.
if(length(args)!=3){
  stop("These arguments must be supplied: input distance summaries, 
       output list of outlier sequences,
       output list of sequence distances.", 
       call.=FALSE)
}

InputFile <- args[1]
OutputFileSequences <- args[2]
OutputFileDistances <- args[3]

# Read in file formats.
source("analysis/FileFormats.R")
source("analysis/PlotThemes.R")

# Read in sequence distances.
Data <- read.table(InputFile, header=TRUE, stringsAsFactors = FALSE)

# Parse the FASTA header.
Data <- Data %>%
  separate(Header, into=GISAIDFASTAHeader, sep=" \\| ") %>%
  mutate(Year=round(as.numeric(CollectionDate)))

# Tidy the data.
Data <- Data %>%
  gather(S, NS, key="Syn", value="Distance")

# Use the S and NS sequence distances to identify sequences that are
# more than three IQRs away from the median using either metric.
# If the IQR is 0, then use 1 for the IQR instead.
Data <- Data %>%
  group_by(Syn, Year) %>%
  mutate(Median=median(Distance),
         IQR=quantile(Distance,0.75)-quantile(Distance,0.25)) %>%
  mutate(Outlier=ifelse(Distance>Median+3*max(IQR,1),"yes",
                        ifelse(Distance<Median-3*max(IQR,1),"yes","no"))) %>%
  ungroup() %>% dplyr::select(-Median, -IQR)

# Export a list of sequences to exclude.
Data <- Data %>% unite(Header, GISAIDFASTAHeader, sep=" | ")
write.table(Data %>% filter(Outlier=="yes") %>% dplyr::select(Header), 
            OutputFileSequences,
            row.names=FALSE, quote=TRUE)

# Export a list of sequence distances.
write.table(Data, OutputFileDistances, row.names=FALSE, quote=TRUE)
