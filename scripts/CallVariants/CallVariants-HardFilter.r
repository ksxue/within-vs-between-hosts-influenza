#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Calls variants from BAM summary data produced by SummarizeBAM version 1.2
# and AnnotateVariants version 1.12.
# Calls variants above a specified frequency threshold in the viral population
# at sites that have coverage above a defined threshold.

library(dplyr)

# Verify that the correct number of arguments are given.
if(length(args)!=4){
  stop("The following arguments must be supplied: input summary file, output variants file,
       threshold frequency, and minimum coverage.", 
       call.=FALSE)
}

VariantFile <- args[1]
OutputFile <- args[2]
ThresholdFreq <- as.numeric(args[3])
MinCoverage <- as.numeric(args[4])

# Read in variant file.
Data <- read.table(VariantFile, header=FALSE, stringsAsFactors = FALSE)
colnames(Data) <- c("Chr","Pos","Base","RefBase","GenomePos","Count","AvgQ","AvgReadPos",
                    "Gene","Codon","RefAA","AltAA","Syn","FourfoldSyn","RefCodon","AltCodon","CodonPos",
                    "Sample")
Data <- Data %>% filter(Gene != "none") %>% group_by(Gene, Pos) %>%
  mutate(Consensus=Base[which.max(Count)], Freq=Count/sum(Count),
         Coverage=sum(Count))

# Call minor alleles that are at or above the threshold frequency
# and that attain the minimum coverage at that site.
# Export called variants.
write.table(Data %>% filter(Base!=Consensus, Freq>ThresholdFreq,
                            Coverage>MinCoverage) %>% 
              dplyr::select(Sample,Chr,GenomePos,Gene,Pos,Consensus,Base,Codon,CodonPos,RefAA,AltAA,
                            RefCodon,AltCodon,Count,AvgQ,AvgReadPos,Coverage,Freq),
            file=OutputFile, quote=FALSE, sep='\t',
            row.names=FALSE, col.names=FALSE)
