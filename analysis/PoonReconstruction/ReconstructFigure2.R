require(tidyverse)
require(cowplot)
require(purrr)
require(readr)
require(foreach)

# Reconstruct Figure 2 from the published Poon et al. paper as closely as possible.
# To do this, first identify variable sites in the H3N2 samples that are identified
# at a frequency above 0.03 at sites with coverage above 200 in both replicates.
# Extract the full read count information for all H3N2 samples at those variable sites.
# Retain only information for samples with both sequencing replicates available.

# List output directory.
outdir <- "analysis/PoonReconstruction/out/"

# Read in the list of variants identified in both replicates of H3N2 influenza samples.
Variants <- read.table("analysis/PoonReconstruction/out/ReplicateVariants-H3N2.txt",
                       header=TRUE, stringsAsFactors = FALSE)

# Retain only HA variants that are observed in more than one
# distinct biological sample.
Variants <- Variants %>% group_by(GenomePos) %>%
  filter(n_distinct(Strain)>1) %>%
  filter(Gene=="4-HA", Codon)

# Read in all of the raw read count files for the H3N2 variants.
subtype <- "H3N2"
path <- "nobackup/DownloadDataCallVariants/flu-Poon/"
files <- dir(path, pattern=paste0(subtype,"-realigned-annotated.summary$"))

# Combine all of the variant files into a single dataframe.
Data <- files %>%
  map(~ read.table(file.path(path, .), 
                   header=FALSE, stringsAsFactors = FALSE)) %>%
  reduce(rbind)
colnames(Data) <- c("Chr","Pos","Base","RefBase","GenomePos","Count","AvgQ","AvgReadPos",
                    "Gene","Codon","RefAA","AltAA","Syn","FourfoldSyn",
                    "Sample")
Data <- Data %>% filter(Gene != "none") %>% group_by(Sample, Gene, Pos) %>%
  mutate(Consensus=Base[which.max(Count)], Freq=Count/sum(Count),
         Coverage=sum(Count))

# Retain only sites of HA variation in two or more samples,
# as determined previously.
DataVariable <- Data %>% filter(GenomePos %in% unique(Variants$GenomePos))

# Retain only samples for which both sequencing replicates are available.
# Import the list of samples that contain both sequencing replicates.
Replicates <- read.table("data/metadata/flu-Poon/Synapse-samplenames-replicates.txt",
                         header=FALSE, stringsAsFactors = FALSE)
colnames(Replicates) <- c("Sample")
DataVariable <- DataVariable %>% filter(Sample %in% Replicates$Sample)

# For each site and each strain, average the variant frequencies
# across the two sequencing replicates.
DataVariable <- DataVariable %>% 
  separate(Sample, into=c("Household","Visit","Replicate")) %>%
  mutate(Strain=paste(Household,Visit,sep="-")) %>%
  dplyr::select(-Household,-Visit)
DataVariable <- DataVariable %>% ungroup() %>%
  group_by(Strain, Chr, Pos, Base, GenomePos,
           Gene, Codon) %>% 
  summarize(Freq=mean(Freq))

# Export the averaged base frequencies at each variable HA site.
write.table(DataVariable, paste0(outdir,"Figure2Frequencies.txt"),
            col.names=TRUE, row.names=FALSE, quote=FALSE)
