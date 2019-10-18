library(tidyverse)

# Import information on file formats.
source("analysis/FileFormats.R")

# Import the metadata for the studies analyzed here.
Metadata <- rbind(
  read.table("data/metadata/flu-Dinis/samplesheet-refs.txt", header=FALSE, stringsAsFactors = FALSE),
  read.table("data/metadata/flu-Debbink/samplesheet-refs.txt", header=FALSE, stringsAsFactors = FALSE),
  read.table("data/metadata/flu-McCrone/samplesheet-refs.txt", header=FALSE, stringsAsFactors = FALSE),
  read.table("data/metadata/flu-Xue-chronic/samplesheet-refs.txt", header=FALSE, stringsAsFactors = FALSE)
)
colnames(Metadata) <- MetadataCols

# Import and clean metadata on the timing of sample collection
# for each sample in this analysis.
# Note that not all studies have sample-specific data.
# Read in the metadata files, interpret the column names,
# and remove information on samples that do not have associated sequencing data.
# Generate separate columns for viral load measured in genome copies/uL
# and Ct, which corresponds to viral load but in an unclear way.
MetadataDPI <- rbind(read.csv("data/metadata/flu-Debbink/metadata-flu-Debbink.data",
                              header=TRUE, stringsAsFactors = FALSE) %>%
                       mutate(Project="flu-Debbink") %>% 
                       dplyr::select(Id, season, Day.of.Infection.sample.collected, Project, Copy_num) %>%
                       mutate(Subtype="H3N2") %>% 
                       dplyr::rename(Season=season, DPI=Day.of.Infection.sample.collected, 
                                     Sample=Id, ViralLoad=Copy_num) %>%
                       # convert the sample names to padded, 3-digit designators like in other analyses
                       mutate(Sample=str_pad(as.character(Sample),3,side="left",pad="0"), Ct=NA) %>%
                       filter(Sample %in% (Metadata %>% filter(Project=="flu-Debbink"))$Sample),
                     read.csv("data/metadata/flu-McCrone/metadata.data",
                              header=TRUE, stringsAsFactors = FALSE) %>%
                       mutate(Project="flu-McCrone") %>%
                       dplyr::select(LAURING_ID, season, DPI, Project, pcr_result, gc_ul) %>%
                       dplyr::rename(Season=season, Subtype=pcr_result, Sample=LAURING_ID,
                                     ViralLoad=gc_ul) %>%
                       filter(Sample %in% (Metadata %>% filter(Project=="flu-McCrone"))$Sample) %>%
                       mutate(Subtype=substr(Subtype, 3, 7)) %>% # format the subtype designation
                       mutate(Subtype=ifelse(Subtype=="H1N1","pdmH1N1",Subtype), Ct=NA))

