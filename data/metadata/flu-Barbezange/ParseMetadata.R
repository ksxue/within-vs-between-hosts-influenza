library(tidyverse)
library(lubridate)

# Read in raw sample metadata.
Metadata <- read.table("data/metadata/flu-Barbezange/metadata-raw.data",
                       header=TRUE, stringsAsFactors = FALSE)
# Unify the spelling on the sample groups.
Metadata <- Metadata %>% 
  mutate(Group=ifelse(Group=="Severe","severe",Group)) %>%
  rename(VirusIdentifier=Sample)

# Import the SRA run table.
SRA <- read.table("data/metadata/flu-Barbezange/SraRunTable.txt",
                  header=TRUE, stringsAsFactors = FALSE, sep="\t")
SRA <- SRA %>% dplyr::select(virus_identifier, host_subject_id) %>%
  rename(VirusIdentifier=virus_identifier,
         Sample=host_subject_id)

# Join the virus identifier to the sample name used
# in my analysis of the sequencing data.
Metadata <- left_join(Metadata, SRA,
                      by=c("VirusIdentifier"))

# Calculate days post-infection for each sample.
Metadata <- Metadata %>% 
  mutate(DPI=ydm(paste(SamplingYear,SamplingDate,sep="-"))-
           ydm(paste(SamplingYear,SymptomOnsetDate,sep="-")))

# Format values of viral load.
Metadata <- Metadata %>%
  mutate(ViralLoad=as.double(ViralLoad))

# Export formatted metadata for the Barbezange et al. study.
write.table(Metadata, "data/metadata/flu-Barbezange/metadata.data",
            row.names=FALSE, quote=FALSE)