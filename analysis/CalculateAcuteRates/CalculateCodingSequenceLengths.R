library(ggplot2)
library(tidyverse)
library(cowplot)
library(foreach)
library(Cairo)
library(Biostrings)
library(seqinr)

# Output directory for generated data.
outdir<-"analysis/CalculateAcuteRates/out/"

# Import information on file formatting.
source("analysis/FileFormats.R")

# Import subtype reference genomes and calculate lengths of coding sequences.
CodingSequenceLengths <- rbind(
  read.table("reference/flu-H3N2/H3N2-Victoria-2011.bed", header=FALSE, stringsAsFactors = FALSE) %>%
    `colnames<-`(BEDCols) %>% mutate(Subtype="H3N2"),
  read.table("reference/flu-pdmH1N1/pdmH1N1-California-2009.bed", header=FALSE, stringsAsFactors = FALSE) %>%
    `colnames<-`(BEDCols) %>% mutate(Subtype="pdmH1N1"),
  read.table("reference/flu-seasonalH1N1/H1N1-Boston-2007.bed", header=FALSE, stringsAsFactors = FALSE) %>%
    `colnames<-`(BEDCols) %>% mutate(Subtype="seasonalH1N1")
)
CodingSequenceLengths$Length <- sapply(CodingSequenceLengths$BlockSizes,
                                       function(x) sum(sapply(strsplit(x,split=","), as.integer)))
CodingSequenceLengths <- CodingSequenceLengths %>%
  dplyr::select(Subtype,Chr,Gene,Length)

# Export list of coding sequence lengths.
write.table(CodingSequenceLengths, paste0(outdir,"CodingSequenceLengths.data"),
            quote=FALSE, row.names=FALSE)
