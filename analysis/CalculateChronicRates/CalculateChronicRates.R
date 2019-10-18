library(tidyverse)
library(broom)

# Import information on coding sequence lengths.
CodingSequenceLengths <- read.table("analysis/CalculateAcuteRates/out/CodingSequenceLengths.data",
                                    header=TRUE, stringsAsFactors = FALSE)
# Import information on excluded samples.
source("analysis/CalculateAcuteRates/IdentifySampleExclusions.R")
# Import information on the number of available sites of each mutation type.
AvailableSites <- read.table("analysis/CalculateAcuteRates/out/AvailableSitesByType.data",
                             header=TRUE, stringsAsFactors = FALSE)

outdir <- "analysis/CalculateChronicRates/out/"

# Import within-host variants above a frequency of 0.5%.
Data <- read.table("analysis/CallWithinHostVariants/variants-annotated-0.005.data.gz",
                   header=TRUE, stringsAsFactors = FALSE)
# Exclude samples from acute infections.
Data <- Data %>% 
  filter(Project=="flu-Xue-chronic")
# Classify variants into mutation classes.
Data <- Data %>%
  mutate(Syn=ifelse(RefAA==AltAA,"S",
                    ifelse(AltAA=="*","Stop","NS")))

# Choose variant-frequency thresholds to analyze.
MinFreqs <- c(0.005,0.01,0.02)

# Write a function to calculate the raw sample divergence
# in each mutation class by summing the variant frequencies
# in each sample in each class of sites.
# This function does not yet perform any normalizations for
# available sites, DPI, etc.
CalculateSampleDivergence <- function(minfreq){
  
  # Filter out all variants above the variant-call threshold.
  DataWithin <- Data %>% 
    filter(Freq>minfreq)
  # Calculate the sample divergence for each gene and mutation class.
  DataWithinDiv <- DataWithin %>%
    group_by(Project, Subtype, Sample, Gene, Syn) %>%
    summarize(Div=sum(Freq))
  # Identify samples that are not represented in this dataframe
  # because they have zero variants.
  # Add them to the dataframe in the form of a dummy value.
  DataWithinDiv <- rbind(DataWithinDiv %>% ungroup(),
                         Metadata %>% filter(Project=="flu-Xue-chronic") %>%
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

# Normalize sample divergence based on the number of available sites.
DataWithinDiv <- DataWithinDivRaw %>%
  mutate(DivPerSite=Div/(Length*PercentSites))

# Export the sample divergence values.
write.table(DataWithinDiv, paste0(outdir, "SampleDivPerSite-chronic.data"),
            quote=FALSE, row.names=FALSE)

# Exclude samples that were determined to have poor replicate correlations
# in the Xue et al. 2017 paper.
DataWithinDiv <- DataWithinDiv %>%
  filter(!(Sample %in% SampleExclusions$Sample))

# Parse the sample names for each patient to extract patient ID and timepoint.
DataWithinDiv <- DataWithinDiv %>%
  mutate(Patient=substr(Sample,1,1),
         Timepoint=as.integer(substr(Sample,2,3)),
         Replicate=as.integer(substr(Sample,nchar(Sample),nchar(Sample))))

# Perform linear regression for each patient, replicate, gene, and mutation type
# to determine the average evolutionary rates.
DataWithinRates <- DataWithinDiv %>%
  group_by(MinFreq, Patient, Replicate, Gene, Syn) %>%
  do(tidy(lm(DivPerSite ~ Timepoint, .))) %>%
  filter(term=="Timepoint") %>%
  dplyr::select(-term, -statistic, -p.value) %>%
  dplyr::rename(MeanDivPerSitePerDay=estimate,
         SEDivPerSitePerDay=std.error)
# Export chronic within-host evolutionary rates by patient.
write.table(DataWithinRates,
            paste0(outdir,"ChronicRatesByPatient.data"),
            quote=FALSE, row.names=FALSE)

# Perform linear regression for each replicate, gene, and mutation type
# while averaging across patients to determine the average evolutionary rates.
DataWithinRatesAllPatients <- 
  DataWithinDiv %>%
  group_by(MinFreq, Replicate, Gene, Syn) %>%
  do(tidy(lm(DivPerSite ~ Timepoint, .))) %>%
  filter(term=="Timepoint") %>%
  dplyr::select(-term, -statistic, -p.value) %>%
  dplyr::rename(MeanDivPerSitePerDay=estimate,
         SEDivPerSitePerDay=std.error)
# Export chronic within-host evolutionary rates across all patients.
write.table(DataWithinRatesAllPatients,
            paste0(outdir,"ChronicRatesAllPatients.data"),
            quote=FALSE, row.names=FALSE)
