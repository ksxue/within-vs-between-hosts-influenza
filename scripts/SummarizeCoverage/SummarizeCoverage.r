#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Summarizes sequencing coverage from BAM summary data produced by SummarizeBAM version 1.21
# and AnnotateVariants version 1.13.

library(tidyverse)

# Verify that the correct number of arguments are given.
if(length(args)!=3){
  stop("The following arguments must be supplied: input summary file, output variants file,
       and minimum coverage.", 
       call.=FALSE)
}

VariantFile <- args[1]
OutputFile <- args[2]
MinCoverage <- as.numeric(args[3])

# Read in variant file.
Data <- read.table(VariantFile, header=FALSE, stringsAsFactors = FALSE)
colnames(Data) <- c("Chr","Pos","Base","RefBase","GenomePos","Count","AvgQ","AvgReadPos",
                    "Gene","Codon","RefAA","AltAA","Syn","FourfoldSyn","RefCodon","AltCodon","CodonPos",
                    "Sample")
# Calculate sequencing coverage, taking care not to double coverage
# due to annotations of overlapping genes.
Data <- Data %>% 
  group_by(Sample, Chr, GenomePos, Base) %>% summarize(Count=mean(Count)) %>%
  ungroup() %>% group_by(Sample, Chr, GenomePos) %>% 
  summarize(Coverage=sum(Count))

# Summarize sequencing coverage at each chromosome.
DataSummary <- Data %>% 
  ungroup() %>% group_by(Sample, Chr) %>% 
  summarize(MeanCov=mean(Coverage), MedCov=median(Coverage),
            IQRCov=quantile(Coverage,0.75)-quantile(Coverage,0.25),
            MinCov=min(Coverage), MaxCov=max(Coverage),
            PercentSitesAboveMin=sum(Coverage>MinCoverage)/n())

# Export coverage summary.
write.table(DataSummary,
            file=OutputFile, quote=FALSE, sep='\t',
            row.names=FALSE, col.names=FALSE)
