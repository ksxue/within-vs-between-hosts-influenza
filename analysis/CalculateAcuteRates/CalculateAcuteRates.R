library(ggplot2)
library(tidyverse)
library(cowplot)
library(foreach)
library(Cairo)
library(broom)

# Output directory for generated data.
outdir<-"analysis/CalculateAcuteRates/out/"

# Import list of sample exclusions.
SampleExclusions <- read.table("analysis/CalculateAcuteRates/out/SampleExclusions.data",
                               header=TRUE, stringsAsFactors = FALSE)
# Import plot themes.
source("analysis/PlotThemes.R")

# Import the data on acute sample divergence in each site category.
# These rates have been normalized for the number of available sites,
# but not the sample collection timing.
Data <- read.table("analysis/CalculateAcuteRates/out/SampleDivPerSite.data",
                   header=TRUE, stringsAsFactors = FALSE)

# Exclude samples according to previously defined criteria.
Data <- Data %>%
  filter(!(Sample %in% SampleExclusions$Sample))
# Exclude seasonal and pandemic H1N1 samples.
Data <- Data %>%
  filter(Subtype=="H3N2")

# Calculate the average DPI for each study among the remaining samples.
DataDPI <- Data %>%
  group_by(Project, Subtype, Sample) %>%
  summarize(DPI=mean(DPI)) %>%
  ungroup() %>% filter(!is.na(DPI)) %>% 
  group_by(Project, Subtype) %>%
  summarize(StudyMeanDPI=mean(DPI))

# Import information on study mean DPI for the Dinis et al. study.
DinisDPI <- read.table("data/metadata/flu-Dinis/metadata-H3N2-2012-2013.txt",
                       header=TRUE, stringsAsFactors = FALSE) %>%
  mutate(Project="flu-Dinis", Subtype="H3N2") %>%
  group_by(Project,Subtype) %>% 
  summarize(StudyMeanDPI=sum((DPI)*NumPatients)/sum(NumPatients))
DataDPI <- bind_rows(DataDPI,DinisDPI)

# Add information on study mean DPI to the sample divergence data.
Data <- left_join(Data, DataDPI, by=c("Project","Subtype"))
# If no DPI information is available for a sample,
# then use the study mean DPI instead.
Data <- Data %>% 
  mutate(DPIInterpolated=ifelse(is.na(DPI),StudyMeanDPI,DPI))

# Calculate average evolutionary rates in each gene
# by normalizing the per-site divergence in each mutation class
# by the DPI for that sample.
# Note that we substitute the study mean DPI when the sample DPI is unavailable.
# Note that this excludes the Dinis et al. study,
# which does not contain sample-specific DPI information.
# We add 2 to each sample DPI to account for the time
# between infection and symptom onset.
# Calculate rates separately for each study.
write.table(Data %>% filter(DPIInterpolated!=-2) %>%
              mutate(DivPerSitePerDay=DivPerSite/(DPIInterpolated+2)) %>%
              group_by(MinFreq, Project, Subtype, Gene, Syn) %>%
              summarize(MeanDivPerSitePerDay=mean(DivPerSitePerDay),
                        SEDivPerSitePerDay=sd(DivPerSitePerDay)/sqrt(n())),
            paste0(outdir,"AcuteRatesByStudy.data"),
            row.names=FALSE, quote=FALSE)
# Calculate rates averaged across samples from all studies.
write.table(Data %>% filter(DPIInterpolated!=-2) %>%
              mutate(DivPerSitePerDay=DivPerSite/(DPIInterpolated+2)) %>%
              group_by(MinFreq, Subtype, Gene, Syn) %>%
              summarize(MeanDivPerSitePerDay=mean(DivPerSitePerDay),
                        SEDivPerSitePerDay=sd(DivPerSitePerDay)/sqrt(n())),
            paste0(outdir,"AcuteRatesAllStudies.data"),
            row.names=FALSE, quote=FALSE)

# Perform linear regression to estimate evolutionary rates
# for each gene and mutation class.
# Remove all samples that do not have a DPI value recorded
# rather than using the interpolated rates.
write.table(Data %>%
              filter(!is.na(DPI)) %>%
              group_by(MinFreq, Subtype, Gene, Syn) %>%
              do(tidy(lm(DivPerSite ~ DPI, .))) %>%
              filter(term=="DPI") %>%
              dplyr::select(-term, -statistic, -p.value) %>%
              dplyr::rename(MeanDivPerSitePerDay=estimate, SEDivPerSitePerDay=std.error),
            paste0(outdir, "AcuteRatesAllStudies-Regression.data"),
            row.names=FALSE, quote=FALSE)
