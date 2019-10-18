library(ggplot2)
library(tidyverse)
library(cowplot)
library(foreach)
library(Cairo)
library(Biostrings)
library(seqinr)

# Output directory for generated data.
outdir<-"analysis/CalculateAcuteRates/out/"

# Import the list of coding sequence lengths.
source("analysis/CalculateAcuteRates/CalculateCodingSequenceLengths.R")

# Import the list of antigenic sites from Wolf et al.
# These sites are listed in H3 numbering, so add 16
# to give the computational numbering that is used here.
WolfAntigenicSites <- read.table("reference/flu-H3N2/H3_wolf_epitope_sites.data",
                                 header=FALSE, stringsAsFactors = FALSE)
colnames(WolfAntigenicSites) <- c("Codon")
WolfAntigenicSites <- WolfAntigenicSites %>%
  mutate(Codon=Codon+16)

# Calculate the length of the H3 HA antigenic region.
# This is based on the Wolf et al. list of antigenic sites.
CodingSequenceLengths <- rbind(CodingSequenceLengths,
                               c("H3N2","4-HA","4-HA-antigenic-Wolf",
                                 3*nrow(WolfAntigenicSites)),
                               c("H3N2","4-HA","4-HA-nonantigenic-Wolf",
                                 (CodingSequenceLengths %>% 
                                    filter(Subtype=="H3N2", Gene=="4-HA"))$Length - 
                                           3*nrow(WolfAntigenicSites)))
CodingSequenceLengths$Length <- as.integer(CodingSequenceLengths$Length)

# Import the list of antigenic sites from Wiley et al.
# and Skehel et al., as annotated by Jesse Bloom.
# These sites are listed in H3 numbering.
# Convert them to computational site numbering.
WileySkehelAntigenicSites <- read.table("reference/flu-H3N2/H3_antigenic_sites_WileySkehel.txt",
                                        header=FALSE, stringsAsFactors = FALSE)
colnames(WileySkehelAntigenicSites) <- c("Mutation","Site")
WileySkehelAntigenicSites <- WileySkehelAntigenicSites %>%
  mutate(AA=substr(Mutation,1,3), Codon=as.integer(substr(Mutation,4,6))) %>%
  dplyr::select(-Mutation,AA) %>%
  mutate(Codon=Codon+16)
WileyAntigenicSites <- WileySkehelAntigenicSites %>%
  filter(Site %in% c("A","B","C","D","E"))

# Calculate the length of the H3 HA antigenic region
# based on the Wiley et al. list of antigenic sites.
CodingSequenceLengths <- rbind(CodingSequenceLengths,
                               c("H3N2","4-HA","4-HA-antigenic-Wiley",
                                 3*n_distinct(WileyAntigenicSites$Codon)),
                               c("H3N2","4-HA","4-HA-nonantigenic-Wiley",
                                 (CodingSequenceLengths %>% 
                                    filter(Subtype=="H3N2", Gene=="4-HA"))$Length - 
                                           3*n_distinct(WileyAntigenicSites$Codon)))
CodingSequenceLengths$Length <- as.integer(CodingSequenceLengths$Length)


# Import the list of NA surface-exposed sites from Bhatt et al. 2012.
# These sites are listed in computational numbering, so no offsets are necessary.
BhattSurfaceExposedSites <- read.table("reference/flu-H3N2/N2_Bhatt_SurfaceExposedSites.data",
                                 header=FALSE, stringsAsFactors = FALSE)
colnames(BhattSurfaceExposedSites) <- c("Codon")

# Calculate the length of the NA surface-exposed sites.
# This is based on the Bhatt et al. list of surface-exposed sites.
CodingSequenceLengths <- rbind(CodingSequenceLengths,
                               c("H3N2","6-NA","6-NA-surface-Bhatt",
                                 3*nrow(BhattSurfaceExposedSites)),
                               c("H3N2","6-NA","6-NA-nonsurface-Bhatt",
                                 (CodingSequenceLengths %>% 
                                    filter(Subtype=="H3N2", Gene=="6-NA"))$Length - 
                                   3*nrow(BhattSurfaceExposedSites)))
CodingSequenceLengths$Length <- as.integer(CodingSequenceLengths$Length)
