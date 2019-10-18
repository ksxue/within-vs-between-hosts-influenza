library(ggplot2)
library(tidyverse)
library(cowplot)

# Import sample metadata.
source("analysis/CalculateAcuteRates/ParseSampleMetadata.R")
source("analysis/FileFormats.R")

# Output directory for generated data.
outdir<-"analysis/CalculateAcuteRates/out/"

set.seed(0)


# Import a list of sample exclusions for the Xue et al. 2017 datas --------


# Import a list of excluded samples from the Xue et al. 2017 dataset.
# These samples have low sequencing coverage, high replicate variability,
# or are technical controls.
XueChronicExclusions <- read.table("data/metadata/flu-Xue-chronic/exclusions.txt",
                                   header=FALSE, stringsAsFactors = FALSE)
colnames(XueChronicExclusions) <- c("Sample")
XueChronicExclusions <- XueChronicExclusions %>%
  mutate(Project="flu-Xue-chronic", Reason="LowCoverageHighVariability")


# Import a list of longitudinal pairs from the McCrone et al. data --------

# Import lists of the samples in the McCrone study that are from the same individual.
# Exclude the first sample from each of these longitudinal pairs,
# since they were typically self-collected swabs.
# (The second sample of the pair was clinically collected.)
McCroneLongitudinal <- read.table("data/metadata/flu-McCrone/LongitudinalSamples.data",
                                  header=TRUE, stringsAsFactors = FALSE)
McCroneLongitudinalEarly <- McCroneLongitudinal %>%
  group_by(ENROLLID) %>% arrange(ENROLLID, collect) %>% 
  filter(row_number()==1) %>% dplyr::rename(Sample=LAURING_ID) %>% ungroup() %>% 
  mutate(Project="flu-McCrone") %>% dplyr::select(Project, Sample) %>%
  mutate(Reason="LongitudinalFirstTimepoint")


# Identify samples that represent plasmid controls ------------------------


# Identify other samples to be excluded from downstream analyses.
# These are mostly plasmid controls that were included in the SRA datasets.
ControlSamples <- Metadata %>%
  filter(Sample %in% c("Brisbane", "Cal-H3N2", "Cali_pool","PC1A","Perth_mp","Vic_pool")) %>%
  dplyr::select(Project, Sample) %>%
  mutate(Reason="PlasmidControl")


# Identify samples with abnormally high within-host variation -------------

# Identify samples with abnormally high within-host variation
# to be excluded from subsequent analyses.
# Do this only for the acute infection samples, since the 
# chronic infection samples from Xue et al. 2017 have exclusions identified separately.
# Identify the top 10% of samples in terms of number of variants above 0.5% frequency
# for each subtype and each study.
DataWithin <- read.table("analysis/CallWithinHostVariants/variants-annotated-0.005.data.gz",
                         header=TRUE, stringsAsFactors = FALSE)
DataWithin <- DataWithin %>% filter(Project!="flu-Xue-chronic")

# Identify the the top X% of samples from each study
# by number of variants at a 0.5% calling threshold.
truncate <- 0.1
# Determine what proportion of samples to remove from each analysis.
NumSamplesToRemove <- DataWithin %>% group_by(Project, Subtype) %>%
  summarize(NumSamples=n_distinct(Sample)) %>%
  mutate(NumSamples=ceiling(truncate*NumSamples)) %>%
  arrange(Project, Subtype)
# Identify the top X% of samples by number of variants from each study.
TopSamples <- DataWithin %>%
  filter(Freq>0.005) %>% group_by(Project, Subtype, Sample) %>%
  summarize(NumVariants=n()) %>% 
  ungroup() %>% group_by(Project, Subtype) %>%
  arrange(Project, Subtype, desc(NumVariants)) %>%
  nest() %>%
  mutate(n=NumSamplesToRemove$NumSamples) %>%
  mutate(Sample=map2(data, n, top_n)) %>%
  dplyr::select(Project, Subtype, Sample) %>% unnest()
# Retain only the project and sample names.
TopSamples <- TopSamples %>% dplyr::select(Project, Sample) %>%
  mutate(Reason="TopNumbersOfVariants")

rm(DataWithin)


# Identify samples with low sequencing coverage ---------------------------

# Import summary of sequencing coverage for each sample.
DataCoverage <- read.table("analysis/CallWithinHostVariants/coverage-annotated.data.gz",
                   header=TRUE, stringsAsFactors = FALSE)
# Focus on samples in acute infections (chronic-infection samples have been assessed separately).
# Omit genes that were not sequenced in the Dinis et al. dataset.
DataCoverage <- DataCoverage %>%
  filter(Project!="flu-Xue-chronic", 
         !(Project=="flu-Dinis" & Chr!="4-HA"))
# Identify samples that have less than 80% coverage on at least one gene.
LowCoverageSamples <- DataCoverage %>%
  filter(PercentSitesAboveMin<0.8) %>%
  dplyr::select(Project,Subtype,Sample) %>% distinct() %>%
  mutate(Reason="LowCoverage")

# Combine the lists of sample exclusions identified through each method --------

# Combine the lists of excluded samples.
SampleExclusions <- rbind(XueChronicExclusions, 
                          McCroneLongitudinalEarly, 
                          ControlSamples,
                          TopSamples,
                          LowCoverageSamples %>% dplyr::select(-Subtype))
SampleExclusions <- SampleExclusions %>% distinct()

# Export the list of sample exclusions.
write.table(SampleExclusions %>% dplyr::select(-Reason), 
            paste0(outdir,"SampleExclusions.data"),
            row.names=FALSE, quote=FALSE)
# Export annotated list of sample exclusions.
write.table(SampleExclusions, 
            paste0(outdir,"SampleExclusions-annotated.data"),
            row.names=FALSE, quote=FALSE)
