library(tidyverse)
library(foreach)

source("analysis/FileFormats.R")
# Import sample metadata.
source("analysis/CalculateAcuteRates/ParseSampleMetadata.R")
# Import list of excluded samples.
source("analysis/CalculateAcuteRates/IdentifySampleExclusions.R")

# List output directory path.
outdir <- "analysis/CompareVariantCallingThresholds/out/"

# Import all variant calls above 0% frequency.
Data <- read.table("nobackup/CallWithinHostVariants/variants-annotated.data.gz",
                   header=TRUE, stringsAsFactors = FALSE)

# Analyze acute variants only.
Data <- Data %>% filter(Project != "flu-Xue-chronic")

# Summarize the numbers of variants at each codon site
# at each variant frequency.
# Exclude samples that were omitted for having large numbers of variants
# or for being plasmid controls or parts of longitudinal pairs.
# Also exclude sites in the M1/M2 and NS1/NEP genes,
# since the interpretation of codon position is unclear.
DataCodonPos <- Data %>% 
  filter(!Sample %in% SampleExclusions$Sample,
         Gene %in% FluGenesNonOverlapping) %>%
  mutate(FreqBin=floor(log10(Freq)*5)/5) %>%
  mutate(CodonPos=CodonPos+1) %>%
  group_by(Subtype, FreqBin, CodonPos) %>% 
  summarize(NumVariants=n()) %>% ungroup() %>%
  group_by(Subtype, FreqBin) %>% 
  mutate(PercentVariants=NumVariants/sum(NumVariants))
# Export the table of variant proportions.
write.table(DataCodonPos, paste0(outdir,"VariantCodonPositions.data"),
            quote=FALSE, row.names=FALSE)

# Summarize the percent NS variants at each variant frequency.
# Exclude samples that were omitted for having large numbers of variants
# or for being plasmid controls or parts of longitudinal pairs.
# Also exclude sites in the M1/M2 and NS1/NEP genes,
# since the interpretation of mutation effect is unclear.
DataSyn <- Data %>% 
  filter(!Sample %in% SampleExclusions$Sample,
         Gene %in% FluGenesNonOverlapping) %>%
  mutate(FreqBin=10^floor(log10(Freq)*2)/2) %>%
  mutate(Syn=ifelse(RefAA==AltAA,"S","NS")) %>%
  group_by(Subtype, FreqBin, Syn) %>% 
  summarize(NumVariants=n()) %>% ungroup() %>%
  group_by(Subtype, FreqBin) %>% 
  mutate(PercentVariants=NumVariants/sum(NumVariants))
# Export the table of variant proportions.
write.table(DataCodonPos, paste0(outdir,"VariantSvsNS.data"),
            quote=FALSE, row.names=FALSE)


# Summarize the numbers of variants in each sample
# at each variant freqency.
# Count each genome position as a single variant,
# i.e. don't double-count variants in overlapping genes.
MinFreqs <- seq(0,0.02,0.001)
DataNumVariants <- foreach(minfreq=MinFreqs, .combine="rbind") %do% {
  DataVariant <- Data %>% 
    filter(Freq>minfreq) %>% mutate(MinFreq=minfreq)
  DataVariant <- DataVariant %>%
    group_by(MinFreq, Project, Subtype, Sample) %>%
    summarize(NumVariants=n_distinct(GenomePos)) %>% ungroup()
  # If samples have no variants, then make the implicit missing values explicit.
  DataVariant <- rbind(DataVariant,
        Metadata %>% filter(Project!="flu-Xue-chronic") %>%
          filter(!(Sample %in% DataVariant$Sample)) %>%
          dplyr::select(Project, Subtype, Sample) %>%
          mutate(NumVariants=0, MinFreq=0))
  return(DataVariant)
}

# Export the number of variants per sample across various thresholds.
write.table(DataNumVariants, paste0(outdir,"NumVariantsPerSample.data"),
            row.names=FALSE, quote=FALSE)
