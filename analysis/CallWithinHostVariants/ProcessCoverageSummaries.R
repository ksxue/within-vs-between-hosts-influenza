#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

require(tidyverse)

source("analysis/FileFormats.R")
SummaryFile <- args[1]
OutputFile <- args[2]
FocusProject <- args[3]

# Read in within-host variants called from all studies.
# Add some necessary metadata like sample subtype.
Data <- read.table(SummaryFile, header=FALSE, stringsAsFactors = FALSE)
colnames(Data) <- c(CovSummaryCols,"Project")

# Determine the subtype of the originating samples.
fluDebbink <- read.table("data/metadata/flu-Debbink/samplesheet-refs.txt",
                         stringsAsFactors = FALSE)
fluDinis <- read.table("data/metadata/flu-Dinis/samplesheet-refs.txt",
                       stringsAsFactors = FALSE)
fluMcCrone <- read.table("data/metadata/flu-McCrone/samplesheet-refs.txt",
                         stringsAsFactors = FALSE)
fluXueChronic <- read.table("data/metadata/flu-Xue-chronic/samplesheet-refs.txt",
                         stringsAsFactors = FALSE)
fluBarbezange <- read.table("data/metadata/flu-Barbezange/samplesheet-refs.txt",
                            stringsAsFactors = FALSE)
fluXueAcute <- read.table("data/metadata/flu-Xue-acute/samplesheet-refs.txt",
                          stringsAsFactors = FALSE)
Metadata <- rbind(fluDebbink, fluDinis, fluMcCrone, fluXueChronic, fluBarbezange, fluXueAcute)
colnames(Metadata) <- MetadataCols
Metadata <- Metadata %>% select(Sample, Project, Subtype) %>%
  filter(Project==FocusProject)
Data <- full_join(Data %>% mutate(Sample=str_pad(Sample, 3, pad="0")),
                  Metadata, by=c("Sample","Project"))


write.table(Data, OutputFile, quote=FALSE, row.names=FALSE)