library(tidyverse)

# Output directory.
outdir <- "analysis/CalculateGlobalRates/out/"

# Import evolutionary rates calculated for different scales
# through different methods.
# Standardize the columns and annotations for each set of rates.
# Columns are Subtype, Gene, Syn, MeanDivPerSitePerDay, SEDivPerSitePerDay, Scale, MinFreq, Type

# Import evolutionary rates in acute infections, 
# calculated separately for each study.
AcuteRatesByStudy <-
  read.table("analysis/CalculateAcuteRates/out/AcuteRatesByStudy.data",
             header=TRUE, stringsAsFactors = FALSE) %>%
  mutate(Type=paste0(Project,"-",MinFreq), Scale="acute") %>% 
  dplyr::select(-Project)
# Import evolutionary rates in acute infections, 
# calculated across all studies using the point method.
AcuteRatesAll <-
  read.table("analysis/CalculateAcuteRates/out/AcuteRatesAllStudies.data",
             header=TRUE, stringsAsFactors = FALSE) %>%
    mutate(Type=paste0("point-all-",MinFreq), Scale="acute")
# Import evolutionary rates in acute infections,
# calculated across all studies using the regression method.
AcuteRatesRegression <-
  read.table("analysis/CalculateAcuteRates/out/AcuteRatesAllStudies-Regression.data",
             header=TRUE, stringsAsFactors = FALSE) %>% 
    mutate(Type=paste0("regress-all-",MinFreq), Scale="acute")

# Import evolutionary rates in chronic infections,
# calculated separately for each patient through linear regression.
# Analyze only rates from the first sequencing replicate.
ChronicRatesByPatient <-
  read.table("analysis/CalculateChronicRates/out/ChronicRatesByPatient.data",
             header=TRUE, stringsAsFactors = FALSE) %>%
    mutate(Subtype="H3N2", Type=paste0(Patient,"-",MinFreq), Scale="chronic") %>%
    filter(Replicate==1) %>% dplyr::select(-Patient, -Replicate)
# Import evolutionary rates in chronic infections,
# calculated across all patients through linear regression.
# Analyze only rates from the first sequencing replicate.
ChronicRatesAll <-
  read.table("analysis/CalculateChronicRates/out/ChronicRatesAllPatients.data",
             header=TRUE, stringsAsFactors = FALSE) %>%
  mutate(Subtype="H3N2", Type=paste0("chronic-",MinFreq), Scale="chronic") %>%
    filter(Replicate==1) %>% dplyr::select(-Replicate)
  
# Import evolutionary rates in global infections,
# calculated from 2007 to the 2017 using linear regression.
GlobalRates2007 <-
  read.table("analysis/CalculateGlobalRates/out/GlobalRates-H3N2-2007.data",
             header=TRUE, stringsAsFactors = FALSE) %>%
    mutate(MinFreq="global", Type=paste0("global"), Subtype="H3N2", Scale="global") %>%
    dplyr::select(-Interval, -MeanDivPerSitePerYear, -SEDivPerSitePerYear) %>%
  filter(Gene %in% c("4-HA","4-HA-antigenic-Wolf","4-HA-nonantigenic-Wolf",
                     "6-NA","6-NA-surface-Bhatt","6-NA-nonsurface-Bhatt","8-NS1"))

# Import evolutionary rates in global infections,
# calculated from 1999 to the 2017 using linear regression.
# This longer interval of time is necessary to accurately estimate
# the evolution of the slower genes.
GlobalRates1999 <-
  read.table("analysis/CalculateGlobalRates/out/GlobalRates-H3N2-1999.data",
             header=TRUE, stringsAsFactors = FALSE) %>%
  mutate(MinFreq="global", Type=paste0("global"), Subtype="H3N2", Scale="global") %>%
  dplyr::select(-Interval, -MeanDivPerSitePerYear, -SEDivPerSitePerYear) %>%
  filter(!(Gene %in% c("4-HA","4-HA-antigenic-Wolf","4-HA-nonantigenic-Wolf",
                       "6-NA","6-NA-surface-Bhatt","6-NA-nonsurface-Bhatt","8-NS1")))
  
# Aggregate all rates calculated through all methods.
RatesAll <- rbind(AcuteRatesByStudy, AcuteRatesAll, AcuteRatesRegression,
                  ChronicRatesByPatient, ChronicRatesAll, GlobalRates2007, GlobalRates1999)

# Export the list of all rates calculated through all methods.
write.table(RatesAll, paste0(outdir,"CombinedRates.data"),
            quote=FALSE, row.names=FALSE)
