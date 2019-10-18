library(tidyverse)
library(cowplot)

source("analysis/PlotThemes.R")

SRA <- read.table("data/metadata/flu-McCrone/SraRunTable.txt",
                   header=TRUE, stringsAsFactors = FALSE, sep="\t")
Metadata <- read.csv("data/metadata/flu-McCrone/metadata.data",
                       header=TRUE, stringsAsFactors = FALSE)

# Download data from the eLife-associated repository on longitudinal samples.
# This is Figure 2, source data 4.
download.file("https://raw.githubusercontent.com/lauringlab/Host_level_IAV_evolution/master/results/Figures/data/Figure_2-source_data_4.csv",
  "data/metadata/flu-McCrone/Fig2SourceData4.csv")
McCroneLongitudinal <- read.csv("data/metadata/flu-McCrone/Fig2SourceData4.csv",
                                header=TRUE, stringsAsFactors = FALSE)

# Using the raw metadata, identify samples that are collected from the same individual
# at different times during the infection.
# This method identifies 49 ENROLLIDs associated with longitudinal samples
# rather than the 43 in the source data downloaded above.
# Use the larger set of longitudinal samples for this meta-analysis.
LongitudinalSamples <-
  (Metadata %>% group_by(ENROLLID) %>%
     filter(n()>1) %>% arrange(desc(DPI)))
# Plot the distribution of collection times for longitudinal samples
# to verify that they match what has previously been reported.
# This plot loosely mimics Figure 2D in McCrone et al., eLife 2018.
p <- LongitudinalSamples %>% arrange(desc(DPI)) %>%
  ggplot() +
  geom_point(aes(x=DPI, y=fct_inorder(ENROLLID))) +
  geom_line(aes(x=DPI,y=fct_inorder(ENROLLID), group=factor(ENROLLID))) +
  ylab("ENROLLID") +
  THEME_ALL
save_plot("data/metadata/flu-McCrone/LongitudinalPairs.png",p)
# Output the list of longitudinal samples.
write.table(LongitudinalSamples,
            "data/metadata/flu-McCrone/LongitudinalSamples.data",
            row.names=FALSE, quote=FALSE)


# # Download data from the eLife-associated repository on transmission pairs.
# This is Figure 3, source data 2.
download.file("https://raw.githubusercontent.com/lauringlab/Host_level_IAV_evolution/master/results/Figures/data/Figure_3-source_data_2.csv",
              "data/metadata/flu-McCrone/Fig3SourceData2.csv")
McCroneTransmission <- read.csv("data/metadata/flu-McCrone/Fig3SourceData2.csv",
                                header=TRUE, stringsAsFactors = FALSE)
# Output the list of samples from donor individuals.
write.table(Metadata %>% filter(ENROLLID %in% McCroneTransmission$ENROLLID1),
            "data/metadata/flu-McCrone/TransmissionDonorSamples.data",
            quote=FALSE, row.names=FALSE)
# Output the list of samples from recipient individuals.
write.table(Metadata %>% filter(ENROLLID %in% McCroneTransmission$ENROLLID2),
            "data/metadata/flu-McCrone/TransmissionRecipientSamples.data",
            quote=FALSE, row.names=FALSE)
