library(ggplot2)
library(tidyverse)
library(cowplot)
library(foreach)
library(Cairo)

# Output directory for generated data.
outdir<-"analysis/CalculateAcuteRates/out/"

# Import the list of coding sequence lengths.
source("analysis/CalculateAcuteRates/CalculateCodingSequenceLengths.R")
# Import the sample metadata.
source("analysis/CalculateAcuteRates/ParseSampleMetadata.R")
# Import sets of antigenic site classifications.
source("analysis/CalculateAcuteRates/ImportAntigenicSites.R")
# Import the set of available sites of each mutation class.
AvailableSites <- read.table("analysis/CalculateAcuteRates/out/AvailableSitesByType.data",
                             header=TRUE, stringsAsFactors = FALSE)

# Import raw variant calls above a frequency of 0.5%.
DataWithinRaw <- read.table("analysis/CallWithinHostVariants/variants-annotated-0.005.data.gz",
                         header=TRUE, stringsAsFactors = FALSE)
# Exclude samples from chronic infections.
DataWithinRaw <- DataWithinRaw %>% 
  filter(Project!="flu-Xue-chronic")
# Classify variants into mutation classes.
DataWithinRaw <- DataWithinRaw %>%
  mutate(Syn=ifelse(RefAA==AltAA,"S",
             ifelse(AltAA=="*","Stop","NS")))

# Annotate HA sites within antigenic regions.
DataWithinRaw <- 
  rbind(DataWithinRaw,
        DataWithinRaw %>% 
          filter(Subtype=="H3N2", Gene=="4-HA",
                 Codon %in% WolfAntigenicSites$Codon) %>%
          mutate(Gene="4-HA-antigenic-Wolf"),
        DataWithinRaw %>% 
          filter(Subtype=="H3N2", Gene=="4-HA",
                 !(Codon %in% WolfAntigenicSites$Codon)) %>%
          mutate(Gene="4-HA-nonantigenic-Wolf"))
DataWithinRaw <- 
  rbind(DataWithinRaw,
        DataWithinRaw %>% 
          filter(Subtype=="H3N2", Gene=="4-HA",
                 Codon %in% WileyAntigenicSites$Codon) %>%
          mutate(Gene="4-HA-antigenic-Wiley"),
        DataWithinRaw %>% 
          filter(Subtype=="H3N2", Gene=="4-HA",
                 !(Codon %in% WileyAntigenicSites$Codon)) %>%
          mutate(Gene="4-HA-nonantigenic-Wiley"))

# Annotate NA sites within surface-exposed regions.
DataWithinRaw <- 
  rbind(DataWithinRaw,
        DataWithinRaw %>% 
          filter(Subtype=="H3N2", Gene=="6-NA",
                 Codon %in% BhattSurfaceExposedSites$Codon) %>%
          mutate(Gene="6-NA-surface-Bhatt"),
        DataWithinRaw %>% 
          filter(Subtype=="H3N2", Gene=="6-NA",
                 !(Codon %in% BhattSurfaceExposedSites$Codon)) %>%
          mutate(Gene="6-NA-nonsurface-Bhatt"))


# Choose variant-frequency thresholds to analyze.
MinFreqs <- c(0.005,0.01,0.02)

# Write a function to calculate the raw sample divergence
# in each mutation class by summing the variant frequencies
# in each sample in each class of sites.
# This function does not yet perform any normalizations for
# available sites, DPI, etc.
CalculateSampleDivergence <- function(minfreq){
  
  # Filter out all variants above the variant-call threshold.
  DataWithin <- DataWithinRaw %>% 
    filter(Freq>minfreq)
  # Calculate the sample divergence for each gene and mutation class.
  DataWithinDiv <- DataWithin %>%
    group_by(Project, Subtype, Sample, Gene, Syn) %>%
    summarize(Div=sum(Freq))
  # Identify samples that are not represented in this dataframe
  # because they have zero variants.
  # Add them to the dataframe in the form of a dummy value.
  DataWithinDiv <- rbind(DataWithinDiv %>% ungroup(),
                         Metadata %>% filter(Project!="flu-Xue-chronic") %>%
                           filter(!(Sample %in% DataWithinDiv$Sample)) %>%
                           dplyr::select(Project, Subtype, Sample) %>%
                           mutate(Gene="4-HA",Syn="NS",Div=0))
  # Fill in zero values for samples, genes, and mutation classes
  # that have zero reported variants using the complete function.
  DataWithinDiv <- DataWithinDiv %>%
    ungroup() %>% group_by(Project, Subtype) %>%
    complete(Sample, Gene, Syn, fill=list(Div=0))
  # Add on a tag signifying the frequency threshold at which
  # these divergence values were calculated.
  DataWithinDiv <- DataWithinDiv %>%
    mutate(MinFreq=minfreq)
  
  return(DataWithinDiv)
}

# Calculate the list of sample divergence values for each
# minimum frequency threshold.
DataWithinDivRaw <- lapply(MinFreqs, CalculateSampleDivergence) %>%
  bind_rows()

# Add metadata about the length of each coding sequence.
DataWithinDivRaw <- left_join(DataWithinDivRaw,
                              CodingSequenceLengths,
                              by=c("Subtype","Gene"))
# Add metadata about the proportion of sites of each mutation type.
DataWithinDivRaw <- left_join(DataWithinDivRaw,
                              AvailableSites, by=c("Syn"))
# Add metadata about sample DPI, viral load, and Ct where available.
DataWithinDivRaw <- left_join(DataWithinDivRaw,
                              MetadataDPI %>% dplyr::select(-Season), 
                              by=c("Project","Subtype","Sample"))

# Normalize sample divergence based on the number of available sites.
DataWithinDiv <- DataWithinDivRaw %>%
  mutate(DivPerSite=Div/(Length*PercentSites))

# Export divergence calculations.
write.table(DataWithinDiv, paste0(outdir, "SampleDivPerSite.data"),
            quote=FALSE, row.names=FALSE)
