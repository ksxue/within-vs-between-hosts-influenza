require(tidyverse)
require(foreach)

set.seed(0)

# Import file formats.
source("analysis/FileFormats.R")


# Summarize NS proportions, global variants -------------------------------

# Read in global variant spectrum data.
DataGlobalMutations <- read.table("analysis/InferGlobalVariants/out/H3N2-unpassaged.mutations.gz",
                         header=TRUE, stringsAsFactors = FALSE)
# Parse global mutation data to infer gene and year.
DataGlobalMutations <- DataGlobalMutations %>%
  separate(Sequences,sep="-",into=c("Subtype","Num","Gene","Year")) %>%
  mutate(Gene=paste0(Num,"-",Gene)) %>% dplyr::select(-Subtype, -Num)
# Calculate the proportion of mutations that are NS
# in each gene, year, and frequency bin.
DataGlobal <- DataGlobalMutations %>%
  mutate(Freq=NumDescendants/TotalSeqs,
         BinFreq=10^(floor(log10(Freq)*2)/2),
         Syn=ifelse(AAAnc==AADer,"S","NS")) %>%
  group_by(Year,Gene,BinFreq, Syn) %>%
  dplyr::count() %>% ungroup() %>%
  spread(Syn,n) %>% replace_na(list(NS=0,S=0)) %>%
  mutate(ProportionNS=NS/(NS+S),
         NumVariants=NS+S)


# Summarize NS proportions, within-host variants --------------------------
  
# Read in within-host variant spectrum data.
DataWithinVariants <- read.table("analysis/CallWithinHostVariants/variants-annotated-0.005.data.gz",
                         header=TRUE, stringsAsFactors = FALSE)
# Read in sample exclusions.
SampleExclusions <- read.table("analysis/CalculateAcuteRates/out/SampleExclusions.data",
                               header=TRUE, stringsAsFactors = FALSE)
# Exclude chronic infection samples and excluded samples.
DataWithinVariants <- DataWithinVariants %>%
  filter(Project!="flu-Xue-chronic", Subtype=="H3N2",
         !(Sample %in% SampleExclusions$Sample))
# Annotate variants based on mutation effect.
# Do not include stop variants as a separate class in this analysis.
DataWithinVariants <- DataWithinVariants %>%
  mutate(Syn=ifelse(RefAA==AltAA,"S","NS"))

# For each sample, summarize the number of S and NS variants.
DataWithinSamples <- DataWithinVariants %>%
  mutate(BinFreq=10^(floor(log10(Freq)*2)/2)) %>%
  group_by(Sample, Gene, Syn, BinFreq) %>%
  dplyr::count() %>%
  spread(Syn,n) %>% replace_na(list(NS=0,S=0)) %>%
  ungroup() %>%
  complete(BinFreq, Gene, Sample, fill=list(NS=0,S=0))

# Bootstrap the samples and calculate the resulting proportion
# of NS variants across the bootstrapped dataset.
DataWithinBootstraps <- foreach(i=seq(1:100), .combine="rbind") %do% {
  DataWithinSamples %>% group_by(BinFreq, Gene) %>%
    sample_frac(1,replace=TRUE) %>%
    summarize(NS=sum(NS),S=sum(S),ProportionNS=NS/(NS+S),
              NumVariants=NS+S)
}


# Combine within- and between-host NS proportions -------------------------

# Combine and export within- and between-host NS proportions.
write.table(
  rbind(DataGlobal %>% ungroup(), 
      DataWithinBootstraps %>% ungroup() %>% mutate(Year="within-host")),
  "analysis/InferGlobalVariants/out/H3N2-unpassaged-NSproportions.data",
  col.names = TRUE, row.names = FALSE, quote = FALSE
  )



# Calculate global fixation NS proportions --------------------------------


# Import proportion of available NS sites.
AvailableSites <- read.table("analysis/CalculateAcuteRates/out/AvailableSitesByType.data",
                             header=TRUE, stringsAsFactors = FALSE)
# Classify stop sites as NS sites.
AvailableNS <- (AvailableSites %>% 
                  spread(Syn,PercentSites) %>%
                  mutate(NS=NS+Stop))$NS[1]

# Import global rates of evolution, i.e. rates of fixation.
GlobalFixationRates <- read.table("analysis/CalculateGlobalRates/out/CombinedRates.data",
                                  header=TRUE, stringsAsFactors = FALSE)
GlobalFixationRates <- GlobalFixationRates %>% filter(Scale=="global")
# Calculate the proportion of NS substitutions relative to total substitutions.
# Take into account the differences in the number of available sites.
GlobalFixationRates <- left_join(GlobalFixationRates,
                                 AvailableSites, by="Syn") %>%
  dplyr::select(-SEDivPerSitePerDay) %>%
  mutate(MeanDivPerDayTotal=MeanDivPerSitePerDay*PercentSites) %>%
  dplyr::select(Gene, Syn, MeanDivPerDayTotal) %>%
  spread(Syn, MeanDivPerDayTotal) %>%
  mutate(ProportionNSFixed=NS/(NS+S))

# Export global fixation rates.
write.table(GlobalFixationRates,
            "analysis/InferGlobalVariants/out/H3N2-global-fixation-rates.data",
            row.names=FALSE, quote=FALSE)
