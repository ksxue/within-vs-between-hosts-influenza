require(tidyverse)
require(cowplot)
require(gplots)
require(ggdendro)
source("analysis/PlotThemes.R")

# Read in pair data.
Data <- read.table("nobackup/ReadPairs/readIDrun_counts.txt",
                   header=FALSE, stringsAsFactors = FALSE)
colnames(Data) <- c("NumReads","Run1","Read1","Run2","Read2")
outdir <- "analysis/ReadPairs/out/"

# For now, ignore run information.
Data <- Data %>% dplyr::select(-Run1, -Run2)

# Process sample names to remove the read pair indicator
# (a 1 or 2 in front of the sample name)
Data <- Data %>% 
  separate(Read1, into=c("Read","Read1"), sep="-", extra="merge") %>%
  dplyr::select(-Read) %>%
  separate(Read2, into=c("Read","Read2"), sep="-", extra="merge") %>%
  dplyr::select(-Read)


# Convert the data from tidy format to a matrix for hierarchical clustering.
Matrix <- Data %>% spread(key=Read2, value=NumReads, fill=0) %>%
  dplyr::select(-Read1)
rownames(Matrix) <- colnames(Matrix)
Matrix <- as.matrix(Matrix)

# Create a dataframe that summarizes, for each pair of samples,
# the number of read pairs that are shared between them,
# regardless of whether the read pairs are 1-2 or 2-1.
DataMerged <- as.data.frame(Matrix + t(Matrix)) %>% rownames_to_column("Sample1") %>%
  gather(key=Sample2, value=NumReads, -Sample1)

# Summarize read pairs into homotypic and heterotypic.
# Homotypic means that both read pairs are found in the same sample file.
# Heterotypic means that the read pairs are split between different sample files.
DataMerged <- DataMerged %>%
  mutate(Type=ifelse(Sample1==Sample2,"homotypic","heterotypic")) %>%
  group_by(Sample1, Type) %>% summarize(NumReads=sum(NumReads)) %>%
  rename(Sample=Sample1)

# Import data on unpaired reads.
# Organize the data into the same format and append it to the dataframe.
Unpaired <- read.table("nobackup/ReadPairs/single_counts.txt.2",
                       header=FALSE, stringsAsFactors = FALSE)
colnames(Unpaired) <- c("NumReads","ReadSample")

# Combine unpaired read 1's and unpaired read 2's into a single count
# for each sample file.
Unpaired <- Unpaired %>% 
  separate(ReadSample, into=c("Read","Sample"), extra="merge") %>%
  dplyr::select(-Read) %>% group_by(Sample) %>%
  summarize(NumReads=sum(NumReads)) %>%
  mutate(Type="unpaired")

# Merge counts for paired and unpaired samples.
DataAll <- bind_rows(DataMerged, Unpaired) %>% arrange(Sample)

# Export a summary of unpaired reads, heterotypic pairs,
# and homotypic pairs.
write.table(DataAll, paste0(outdir,"ReadSummary.txt"),
            row.names=FALSE, quote=FALSE)
