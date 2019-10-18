require(tidyverse)
require(cowplot)
require(ggplot2)
require(purrr)
require(readr)
require(broom)
require(Cairo)

source("analysis/PlotThemes.R")
source("analysis/FileFormats.R")
source("analysis/CalculateAcuteRates/ImportAntigenicSites.R")

options(bitmapType='cairo')

args = commandArgs(trailingOnly=TRUE)

# Verify that the correct number of arguments is given.
if(length(args)!=5){
  stop("These arguments must be supplied: input directory,
       reference sequence year, year interval for calculation,
       years following the reference to exclude from the analysis,
       output directory.", 
       call.=FALSE)
}

INDIR <- args[1]
REFYEAR <- as.numeric(args[2]) # year of reference sequence from which distances were calculated
INTERVAL <- as.numeric(args[3]) # lengths in years of interval over which distances should be calculated
TRIM <- as.numeric(args[4]) # years to exclude from the analysis following the reference sequence
outdir <- args[5] 

# Read in the files tallying S and NS differences for each H3N2 gene over time.
DataPath <- INDIR 
files <- dir(DataPath, pattern=paste0("*distances-",REFYEAR,"*"))
Data <- files %>%
  map(~ read.table(file.path(DataPath, .),
                   header=TRUE, stringsAsFactors = FALSE)) %>%
  purrr::reduce(rbind)

# Remap the gene names to reflect genome ordering.
Data <- Data %>% ungroup() %>%
  mutate(Segment=ifelse(is.na(Segment),"NA",Segment))
SegmentToGene <- c("1-PB2","2-PB1","3-PA","4-HA",
                   "4-HA-antigenic-Wolf","4-HA-nonantigenic-Wolf",
                   "5-NP","6-NA","6-NA-surface-Bhatt","6-NA-nonsurface-Bhatt",
                   "7-M1","8-NS1")
names(SegmentToGene) <- c("PB2","PB1","PA","HA",
                          "HA-antigenic-Wolf","HA-nonantigenic-Wolf",
                          "NP","NA","NA-surface-Bhatt","NA-nonsurface-Bhatt","MP","NS")
Data <- Data %>% mutate(Gene=SegmentToGene[Segment])

# Remove sequences whose distances are listed as NA
Data <- Data %>% filter(!is.na(S), !is.na(NS))
# Remove sequences whose year is uninterpretable or before the reference year.
# Also remove sequences collected in the specified interval following the reference year.
Data <- Data %>% filter(!is.na(Year), Year >= REFYEAR+TRIM)
# Add information about sequence subtypes.
Data <- Data %>% mutate(Subtype="H3N2")

# Tidy the data.
Data <- Data %>% 
  gather(S, NS, key="Syn", value="Distance") %>%
  dplyr::select(-Segment, -DNAAccession, -RefYear)

# Group sequences into 3-year intervals.
# Mark outliers that are far from the median
# across the 3-year window.
Data <- Data %>% mutate(Interval3=floor((Year-REFYEAR)/3)) %>%
  group_by(Gene, Syn, Interval3) %>%
  mutate(Median=median(Distance),
         IQR=quantile(Distance,0.75)-quantile(Distance,0.25)) %>%
  mutate(Outlier=ifelse(Distance>Median+3*max(IQR,1),"yes",
                 ifelse(Distance<Median-3*max(IQR,1),"yes","no"))) %>%
  ungroup() %>% dplyr::select(-Interval3, -Median, -IQR)

# Output raw sequence distances,
# with certain sequences annotated as outliers.
write.table(Data, paste0(outdir,REFYEAR,"-SequenceDistances-raw.data"),
            quote=FALSE, row.names=FALSE)

# Remove outliers from subsequent analyses.
DataCleaned <- Data %>% filter(Outlier=="no") %>% ungroup() %>%
  dplyr::select(-Outlier)

# Normalize distance so that it is calculated per available site instead of
# per gene.
# Import information about the coding sequence lengths to do this.
DataCleaned <- left_join(DataCleaned,
                         CodingSequenceLengths %>%
                           dplyr::select(Subtype, Gene, Length),
                         by=c("Subtype","Gene"))
# Also import information on the ratio of available sites for each mutation type.
AvailableSites <- read.table("analysis/CalculateAcuteRates/out/AvailableSitesByType.data",
                             header=TRUE, stringsAsFactors = FALSE)
DataCleaned <- left_join(DataCleaned, AvailableSites, by=c("Syn"))
# Normalize the sequence distances.
DataCleaned <- DataCleaned %>% 
  mutate(DivPerSite=Distance/(Length*PercentSites)) %>%
  dplyr::select(-Length, -PercentSites)

# Export the sequence distances once outliers have been removed
# and the distances have been normalized to the number of available sites.
write.table(DataCleaned, paste0(outdir,REFYEAR,"-SequenceDistances.data"),
            row.names=FALSE, quote=FALSE)

# Perform a linear regression of evolutionary distance by
# sequence collection date.
# Do this for intervals of the specified length.
# Exclude intervals in which fewer than half the years have sequences.
DataCleanedStats <- DataCleaned %>% ungroup() %>%
  mutate(Interval=REFYEAR + INTERVAL*floor((Year-REFYEAR)/INTERVAL)) %>%
  group_by(Gene, Syn, Interval) %>%
  filter(n_distinct(Year)>INTERVAL/2) %>%
  do(tidy(lm(DivPerSite ~ CollectionDate, .)))
DataCleanedStats <- DataCleanedStats %>%
  filter(term=="CollectionDate") %>%
  dplyr::select(-term, -statistic, -p.value) %>%
  dplyr::rename(MeanDivPerSitePerYear=estimate, SEDivPerSitePerYear=std.error) %>%
  mutate(MeanDivPerSitePerDay=MeanDivPerSitePerYear/365,
         SEDivPerSitePerDay=SEDivPerSitePerYear/365)

# Output the calculated rates of global evolution.
write.table(DataCleanedStats,
            file=paste0(outdir,"GlobalRates-H3N2-",REFYEAR,".data"),
            quote=FALSE, row.names=FALSE)
