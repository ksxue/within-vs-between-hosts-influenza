require(tidyverse)
require(cowplot)
require(ggdendro)
require(MASS)
require(scales)
require(seqinr)

outdir <- "analysis/figures/"
outdirpres <- "analysis/figures/presentation/"
# Import palettes and themes.
source("analysis/PlotThemes.R")

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}


# Determine chromosome lengths and annotations ----------------------------

# Determine the lengths of each gene segment sequenced for each subtype.
GenomeH3N2 <- as.data.frame(
  sapply(read.fasta("reference/flu-H3N2/H3N2-Victoria-2011.fasta"), length)) %>%
  rownames_to_column(var="Chr") %>%
  `colnames<-`(c("Chr","Length")) %>%
  mutate(End=cumsum(Length)) %>% mutate(Start = lag(End, n=1, default=0)) %>%
  mutate(Midpoint=(Start+End)/2, Subtype="H3N2")

GenomepdmH1N1 <- as.data.frame(
  sapply(read.fasta("reference/flu-pdmH1N1/pdmH1N1-California-2009.fasta"), length)) %>%
  rownames_to_column(var="Chr") %>%
  `colnames<-`(c("Chr","Length")) %>%
  mutate(End=cumsum(Length)) %>% mutate(Start = lag(End, n=1, default=0)) %>%
  mutate(Midpoint=(Start+End)/2, Subtype="pdmH1N1")

GenomeCoordinates <- rbind(GenomeH3N2, GenomepdmH1N1)
GenomeAverage <- full_join(GenomeH3N2, GenomepdmH1N1, by=c("Chr")) %>%
  mutate(Midpoint=(Midpoint.x+Midpoint.y)/2)

# Reconstruct Figure 2 ----------------------------------------------------

# Set palette and sample ordering for this figure.
DNAPALETTE <- c("darkgoldenrod2","firebrick1","forestgreen","dodgerblue")
H3N2SampleOrder <- c("707-V10","781-V10","671-V10","720-V10","720-V20","720-V31","720-V21",
                     "756-V10","662-V10","669-V10","739-V22","739-V23","739-V10","739-V33",
                     "769-V10","747-V10","747-V22","725-V10","710-V10","695-V10","729-V10",
                     "672-V10","714-V10","733-V10","755-V31","726-V10","734-V20","734-V32",
                     "734-V10","692-V10","737-V10","770-V10","689-V22","689-V20","689-V10",
                     "750-V10","741-V10","764-V10","752-V10","703-V10","763-V23","763-V10",
                     "724-V10","724-V20","736-V10","688-V10")
H3N2SampleColor <- c("white","white","white",
                     "cadetblue","cadetblue","cadetblue","cadetblue",
                     "white","white","white",
                     "darkorange1","darkorange1","darkorange1","darkorange1","white",
                     "darkorchid1","darkorchid1",
                     "white","white","white","white","white","white","white","white","white",
                     "violetred1","violetred1","violetred1",
                     "white","white","white",
                     "slateblue1","slateblue1","slateblue1",
                     "white","white","white","white","white",
                     "goldenrod1","goldenrod1",
                     "white","white","white","white")

# Import the raw variant frequency data for Figure 2.
# This data is based on single-end mapping of the Synapse raw data files.
# I identified all sites in H3N2 HA at which variants were present at a frequency
# of at least 0.03 and at sites with coverage at least 200
# in both sequencing replicates.
# I then extracted the raw read count data for those sites from all H3N2 samples,
# and averaged their frequencies in the two sequencing replicates.
# Note that samples without two sequencing replicates are excluded
# from this analysis.
Data <- read.table("analysis/PoonReconstruction/out/Figure2Frequencies.txt",
                   header=TRUE, stringsAsFactors = FALSE)
Data$Strain <- factor(Data$Strain, levels=H3N2SampleOrder)
Data$HouseholdColor <- sapply(Data$Strain, function(x) 
                              H3N2SampleColor[which(H3N2SampleOrder==x)][1])
Codons <- Data %>% dplyr::select(Codon, GenomePos) %>% distinct() %>% filter(Codon<336)
GenomePosToCodons <- as.character(Codons$Codon)
names(GenomePosToCodons) <- Codons$GenomePos

p <- ggplot(Data %>% filter(Freq>0.01, !is.na(Strain), Codon<336)) +
  geom_bar(aes(x=factor(Strain), y=1, fill=(HouseholdColor)), stat="identity",
           width=1, alpha=0.2) +
  geom_hline(yintercept=0) +
  geom_line(aes(x=Strain, y=Freq, color=factor(Base), group=factor(Base))) +
  geom_point(aes(x=Strain, y=Freq, color=factor(Base),
                 shape=factor(Freq>0.03)), size=0.8, fill="white") +
  facet_wrap(~GenomePos, ncol=2, dir="v", strip.position="right",
             labeller=labeller(GenomePos=GenomePosToCodons),
             scale="free_y") +
  scale_shape_manual(values=c(1,19), name="Frequency",
                     labels=c(">1% and <3%",">3%")) +
  scale_color_manual(values=DNAPALETTE, name="") +
  scale_y_continuous(expand=c(0,0), limits=c(0,1)) +
  theme(axis.text.x=element_text(size=4, angle=90, vjust=0.5),
        axis.text.y=element_text(size=4),
        axis.title=element_text(size=7),
        strip.text=element_text(size=5),
        legend.text=element_text(size=5),
        legend.position=c(0.65,-0.04),
        legend.direction="horizontal",
        legend.title=element_text(size=5),
        legend.spacing=unit(0,"cm"),
        legend.margin=margin(t = 10, r = 0, b = 0, l = 0, unit = "pt"),
        legend.box="vertical") +
  scale_fill_identity() +
  xlab("Study subject") + ylab("Frequency (%)")
p
save_plot(paste0(outdir,"ReconstructedFigure2.pdf"),p,
          base_width=7, base_height=5)
save_plot(paste0(outdir,"ReconstructedFigure2.jpg"),p,
          base_width=7, base_height=5)

p <- ggplot(Data %>% filter(Freq>0.01, !is.na(Strain), Codon==335)) +
  geom_bar(aes(x=factor(Strain), y=1, fill=(HouseholdColor)), stat="identity",
           width=1, alpha=0.2) +
  geom_hline(yintercept=0) +
  geom_line(aes(x=Strain, y=Freq, color=factor(Base), group=factor(Base))) +
  geom_point(aes(x=Strain, y=Freq, color=factor(Base),
                 shape=factor(Freq>0.03)), size=0.8, fill="white") +
  scale_shape_manual(values=c(1,19), name="Frequency",
                     labels=c(">1% and <3%",">3%")) +
  scale_color_manual(values=c(DNAPALETTE[1],DNAPALETTE[3]), name="") +
  theme(axis.text.x=element_text(size=4, angle=90, vjust=0.5),
        axis.text.y=element_text(size=4),
        axis.title=element_text(size=7),
        strip.text=element_text(size=5),
        legend.text=element_text(size=5),
        legend.position="right",
        legend.direction="vertical",
        legend.title=element_text(size=5),
        legend.spacing=unit(0,"cm"),
        legend.margin=margin(t = 0, r = 0, b = 10, l = 0, unit = "pt"),
        legend.box="vertical") +
  scale_fill_identity() +
  scale_y_continuous(expand=c(0,0), limits=c(0,1)) +
  xlab("Study subject") + ylab("Frequency (%)")
p
save_plot(paste0(outdir,"ReconstructedFigure2-codon335.pdf"),p,
          base_width=5, base_height=1.5)


# Analyze reconstructed read pairs ----------------------------------------

# Read in pair data.
DataPairs <- read.table("analysis/ReadPairs/out/readIDrun_counts.txt",
                   header=FALSE, stringsAsFactors = FALSE)
colnames(DataPairs) <- c("NumReads","Run1","Read1","Run2","Read2")

# For now, ignore run information.
Data <- DataPairs %>% dplyr::select(-Run1, -Run2)

# Process sample names to remove the read pair indicator
# (a 1 or 2 in front of the sample name)
Data <- Data %>% 
  separate(Read1, into=c("Read","Read1"), sep="-", extra="merge") %>%
  dplyr::select(-Read) %>%
  separate(Read2, into=c("Read","Read2"), sep="-", extra="merge") %>%
  dplyr::select(-Read)

# Create a square data frame of read pair values,
# filling in missing values with zeroes.
Data <- Data %>% spread(key=Read2, value=NumReads, fill=0) %>%
  gather(key=Read2, value=NumReads, -Read1)

# Convert the data from tidy format to a matrix for hierarchical clustering.
Matrix <- Data %>% spread(key=Read2, value=NumReads, fill=0) %>%
  dplyr::select(-Read1)
rownames(Matrix) <- colnames(Matrix)
Matrix <- as.matrix(Matrix)

# Perform hierarchical clustering and plot the resulting heatmap of read pairs
# using the ordering of samples based on similarity from the dendrogram.
Hclust <- hclust(dist(scale(Matrix)))
HclustDendro <- as.dendrogram(Hclust)
Matrix <- as.data.frame(Matrix[order.dendrogram(HclustDendro),
                               order.dendrogram(HclustDendro)],
                        stringsAsFactors = FALSE)
DataMatrix <- Data %>%
  mutate(Read1=factor(Read1, levels=colnames(Matrix)),
         Read2=factor(Read2, levels=colnames(Matrix)))
# Plot read pair matches.
p <- ggplot(DataMatrix) + geom_tile(aes(x=Read1, y=Read2, fill=log10(NumReads+1))) +
  scale_fill_gradient(low = "white", high = PALETTE[2],
                      name="Log 10 number\nof reads") +
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  THEME_ALL +
  theme(axis.ticks=element_blank(),
        axis.text=element_blank(),
        legend.position = c(0.5,0.1),
        legend.key.size=unit(0.3,"cm"),
        legend.direction="horizontal") +
  xlab("Sample containing read 1") + ylab("Sample containing read 2") +
  annotate(geom="text", x=28, y=60, label="@SOLEXA4_0078:8", size=1.8) +
  annotate(geom="text", x=85, y=120, label="@SOLEXA4_0078:7", size=1.8) +
  annotate(geom="text", x=55, y=140, label="pandemic H1N1", size=2, fontface="bold") +
  annotate(geom="segment", x=2, xend=115, y=135, yend=135) +
  annotate(geom="text", x=28, y=65, label="replicate 1", size=1.8) +
  annotate(geom="text", x=85, y=125, label="replicate 2", size=1.8) +
  annotate(geom="text", x=140, y=113, label="@SOLEXA4_0079:3", size=1.8) +
  annotate(geom="text", x=190, y=160, label="@SOLEXA4_0078:1", size=1.8) +
  annotate(geom="text", x=165, y=93, label="H3N2", size=2, fontface="bold") +
  annotate(geom="segment", x=118, xend=220, y=98, yend=98) +
  annotate(geom="text", x=140, y=108, label="replicate 1", size=1.8) +
  annotate(geom="text", x=190, y=155, label="replicate 2", size=1.8)
p
save_plot(paste0(outdir,"ReadPairs-heatmap.pdf"),p,base_width=3.75, base_height=3.75)
save_plot(paste0(outdirpres,"ReadPairs-heatmap.pdf"),p + THEME_PRES,
          base_width=7, base_height=6)
# Output heatmap with labels for easy reading of sample names.
save_plot(paste0(outdirpres,"ReadPairs-heatmap-key.pdf"),
          p+theme(axis.text=element_text(size=3),
                  axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)),
          base_width=10, base_height=10)

# Read in summary of unpaired and paired reads.
Data <- read.table("analysis/ReadPairs/out/ReadSummary.txt",
                   header=TRUE, stringsAsFactors = FALSE)
Data$Type <- factor(Data$Type, levels=c("unpaired","heterotypic","homotypic"))

# Output the total number of unpaired, heterotypic, and homotypic reads
# across all samples.
write.table(Data %>% group_by(Type) %>% summarize(NumReads=sum(NumReads)) %>%
              ungroup() %>% mutate(PercentReads=NumReads/sum(NumReads)),
            paste0(outdir, "TotalReadCounts.txt"),
            quote=FALSE, row.names=FALSE)

# Order samples based on total number of reads.
Data$Sample <- factor(Data$Sample,
                      levels=(Data %>% group_by(Sample) %>% 
                                summarize(NumReads=sum(NumReads)) %>%
                                arrange(NumReads))$Sample)

# Plot a bar plot showing what proportion of each sample
# consists of homotypic paired reads.
p <- ggplot(Data %>% group_by(Sample)) + 
  geom_bar(aes(x=Sample, y=NumReads, fill=factor(Type)), 
           stat="identity", width=1) +
  scale_fill_manual(values=PALETTE[-1], name="Read type",
                    labels=c("unpaired","paired, different file",
                             "paired, same file")) +
  THEME_ALL +
  theme(axis.text.x=element_blank(),
        legend.position=c(0.1,0.8),
        axis.ticks.x=element_blank(),
        legend.key.size=unit(0.3,"cm")) +
  scale_y_continuous(expand=c(0,0), labels =comma) +
  xlab("Sample") + ylab("Number of reads")
p
save_plot(paste0(outdir,"ReadPairs-summary.pdf"),p,
          base_width=3.25, base_height=2)
save_plot(paste0(outdirpres,"ReadPairs-summary.pdf"),p+THEME_PRES,
          base_width=7)



# Compare within-host diversity across different datasets -----------------

# Read in the number of variants in each sample analyzed in this study.
# Note that this differs from the list of variants analyzed below
# because this analysis accounts for samples in which no variants
# were identified.

Data <- read.table("analysis/CompareDiversity/out/variants-per-sample-all.txt",
                   header=TRUE, stringsAsFactors = FALSE) %>%
  mutate(StudySubtype=paste(Study,Subtype,sep="\n"),
         StudySample=paste(Study,Sample,sep="\n"))
Data$StudySubtype <- factor(Data$StudySubtype,
                            levels=(c("flu-Dinis\nH3N2","flu-Dinis\npdmH1N1","flu-Debbink\nH3N2",
                                         "flu-McCrone\nH3N2","flu-McCrone\npdmH1N1",
                                         "flu-Poon\nH3N2","flu-Poon\npdmH1N1")))
# Format study names to display the publication information.
PubYears <- c("2016","2016","2017","2018")
names(PubYears) <- c("flu-Poon","flu-Dinis","flu-Debbink","flu-McCrone")
Data <- Data %>% 
  separate(Study, into=c("Virus","Author"), sep="-", remove=FALSE) %>%
  mutate(PubYear=PubYears[Study]) %>%
  mutate(DisplayStudySubtype=paste0(Author," et al.\n",PubYear,"\n",Subtype))
Data$DisplayStudySubtype <- factor(Data$DisplayStudySubtype,
                                   levels=(c("Dinis et al.\n2016\nH3N2", "Dinis et al.\n2016\npdmH1N1",
                                            "Debbink et al.\n2017\nH3N2", "McCrone et al.\n2018\nH3N2",
                                            "McCrone et al.\n2018\npdmH1N1","Poon et al.\n2016\nH3N2",
                                            "Poon et al.\n2016\npdmH1N1")))


# Plot the number of variants per sample identified in each study.
p_NumVariants <- ggplot(Data) + 
  geom_boxplot(aes(x=DisplayStudySubtype, y=NumVariants/(BpSequenced/1000)),
               outlier.size=0.5, size=0.5) +
  xlab("Study") + ylab("Number of within-host variants\nper kB sequenced") +
  THEME_ALL #+ theme(axis.text.x=element_text(angle=90), vjust=0.5, hjust=0.5)
p_NumVariants
save_plot(paste0(outdir, "CompareDiversity-NumVariants.pdf"),p_NumVariants,
          base_width=4.75, base_height=3.5)
save_plot(paste0(outdir, "CompareDiversity-NumVariants.jpg"),p_NumVariants,
          base_width=4.75, base_height=3.5)
save_plot(paste0(outdir, "CompareDiversity-NumVariants-zoom.pdf"),p_NumVariants + ylim(0,10))

# Read in variants identified in all datasets.
Data <- read.table("analysis/CompareDiversity/out/variants-all.txt",
                   header=TRUE, stringsAsFactors = FALSE) %>%
  mutate(StudySubtype=paste(Study,Subtype,sep="\n"),
         StudySample=paste(Study,Sample,sep="\n"))
Data$StudySubtype <- factor(Data$StudySubtype,
                            levels=rev(c("flu-Dinis\nH3N2","flu-Dinis\npdmH1N1","flu-Debbink\nH3N2",
                                         "flu-McCrone\nH3N2","flu-McCrone\npdmH1N1",
                                         "flu-Poon\nH3N2","flu-Poon\npdmH1N1")))

# Format study names to display the publication information.
PubYears <- c("2016","2016","2017","2018")
names(PubYears) <- c("flu-Poon","flu-Dinis","flu-Debbink","flu-McCrone")
Data <- Data %>% 
  separate(Study, into=c("Virus","Author"), sep="-", remove=FALSE) %>%
  mutate(PubYear=PubYears[Study]) %>%
  mutate(DisplayStudySubtype=paste0(Author," et al. ",PubYear,", ",Subtype))
Data$DisplayStudySubtype <- factor(Data$DisplayStudySubtype,
                                   levels=c("Dinis et al. 2016, H3N2", "Dinis et al. 2016, pdmH1N1",
                                            "Debbink et al. 2017, H3N2","McCrone et al. 2018, H3N2",
                                            "McCrone et al. 2018, pdmH1N1","Poon et al. 2016, H3N2",
                                            "Poon et al. 2016, pdmH1N1"))

# Output the number of genome sites in each study and subtype
# at which more than half of samples display within-host variation.
SharedVariation <- Data %>% group_by(StudySubtype, NumSamplesSequenced, GenomePos) %>%
  summarize(NumVariableSamples=n_distinct(Sample)) %>% ungroup() %>%
  group_by(StudySubtype, NumSamplesSequenced) %>%
  summarize(NumSharedVariableSites=sum(NumVariableSamples > (NumSamplesSequenced)/2)) %>%
  ungroup() %>% separate(StudySubtype, into=c("Study", "Subtype"), sep="\n")
write.table(SharedVariation, paste0(outdir, "SitesOfSharedVariation.txt"),
            quote=FALSE, row.names=FALSE)

# Plot shared variation across the genome for each study.
p_Shared <- Data %>% group_by(DisplayStudySubtype, NumSamplesSequenced, GenomePos, Chr) %>% 
  summarize(NumSamples=n_distinct(Sample)) %>%
  ggplot() + geom_bar(aes(x=GenomePos, y=NumSamples/NumSamplesSequenced,
                          color=factor(Chr)), stat="identity") +
  facet_wrap(~DisplayStudySubtype, ncol=1, scales="free") +
  ylim(0,1) + xlab("Genome position (bp)") + 
  ylab("Proportion of samples in study\nwith within-host variation") +
  scale_color_manual(values=rep(c(PALETTE[1],PALETTE[2]),4), name="Gene") +
  guides(color=FALSE) +
  scale_x_continuous(limits=c(0,13600), expand=c(0,0),
                     breaks=GenomeAverage$Midpoint,
                     labels=c("PB2","PB1","PA","HA","NP","NA","M","NS")) +
  scale_y_continuous(limits=c(0,1), breaks=c(0,0.5,1)) +
  THEME_ALL
p_Shared
save_plot(paste0(outdir,"CompareDiversity-SharedVariants.pdf"), p_Shared,
          base_width=7, base_height=6)
save_plot(paste0(outdirpres,"CompareDiversity-SharedVariants.pdf"), 
          p_Shared+THEME_PRES + theme(axis.text=element_text(size=7),
                               strip.text=element_text(size=11, face="bold")),
          base_width=7, base_height=6)

p_Compare <- plot_grid(p_NumVariants + coord_flip(), p_Shared, 
                       ncol=1, rel_heights=c(2,5,5), labels=c("A","B"),
                       label_size=8)
p_Compare
save_plot(paste0(outdir, "CompareDiversity.pdf"), p_Compare,
          base_width=7, base_height=7)


# Analyze read 1 and read 2 separately ------------------------------------

# Import the raw variant frequency data for Figure 2.
# This particular figure analyzes read 1 and read 2 separately..
# I identified all sites in H3N2 HA at which variants were present at a frequency
# of at least 0.03 and at sites with coverage at least 200
# in both sequencing replicates of the ORIGINAL dataset
# (containing both read 1 and read 2 reads).
# I then extracted the raw read count data for those sites from all H3N2 samples,
# and averaged their frequencies in the two sequencing replicates.
# Note that samples without two sequencing replicates are excluded
# from this analysis.
Data <- read.table("analysis/CompareR1R2/out/Figure2Frequencies.txt",
                   header=TRUE, stringsAsFactors = FALSE)
Data$Strain <- factor(Data$Strain, levels=H3N2SampleOrder)
Data$HouseholdColor <- sapply(Data$Strain, function(x) 
  H3N2SampleColor[which(H3N2SampleOrder==x)][1])
Codons <- Data %>% dplyr::select(Codon, GenomePos) %>% distinct() %>% filter(Codon<336)
GenomePosToCodons <- as.character(Codons$Codon)
names(GenomePosToCodons) <- Codons$GenomePos

PlotVariantFreqs <- function(read){
  p <- ggplot(Data %>% filter(Freq>0.01, !is.na(Strain), Codon<336, Read==read)) +
    geom_bar(aes(x=factor(Strain), y=1, fill=(HouseholdColor)), stat="identity",
             width=1, alpha=0.2) +
    geom_hline(yintercept=0) +
    geom_line(aes(x=Strain, y=Freq, color=factor(Base), group=factor(Base))) +
    geom_point(aes(x=Strain, y=Freq, color=factor(Base),
                   shape=factor(Freq>0.03)), size=0.8, fill="white") +
    facet_wrap(~GenomePos, ncol=2, dir="v", strip.position="right",
               labeller=labeller(GenomePos=GenomePosToCodons),
               scale="free_y") +
    scale_shape_manual(values=c(1,19), name="Frequency",
                       labels=c(">1% and <3%",">3%")) +
    scale_color_manual(values=DNAPALETTE, name="") +
    scale_y_continuous(expand=c(0,0), limits=c(0,1)) +
    theme(axis.text.x=element_text(size=4, angle=90, vjust=0.5),
          axis.text.y=element_text(size=4),
          axis.title=element_text(size=7),
          strip.text=element_text(size=5),
          legend.text=element_text(size=5),
          legend.position=c(0.65,-0.04),
          legend.direction="horizontal",
          legend.title=element_text(size=5),
          legend.spacing=unit(0,"cm"),
          legend.margin=margin(t = 10, r = 0, b = 0, l = 0, unit = "pt"),
          legend.box="vertical") +
    scale_fill_identity() +
    xlab("Study subject") + ylab("Frequency (%)")
  p
  save_plot(paste0(outdir,"ReconstructedFigure2-",read,".pdf"),p,
            base_width=7, base_height=5)
}

sapply(c("R1","R2"), PlotVariantFreqs)


# Analyze diversity on read 1 and read 2 separately -----------------------

# Read in the number of variants in each sample.
# This number is calculated separately for read 1 alone, read 2 alone,
# and both read 1 and read 2 together.

Data <- read.table("analysis/CompareR1R2/out/variants-per-sample-all.txt",
                   header=TRUE, stringsAsFactors = FALSE) %>%
  mutate(StudySubtype=paste(Study,Subtype,sep="\n"),
         StudySample=paste(Study,Sample,sep="\n"))
# Format analysis names for display.
Reads <- c("all reads", "read 1\nonly","read 2\nonly")
names(Reads) <- c("R1R2","R1","R2")
Data <- Data %>% mutate(DisplayRead=Reads[Read])
Data$DisplayRead <- factor(Data$DisplayRead,
                           levels=Reads)

# Plot the number of variants per sample identified in each study.
p_NumVariants <- ggplot(Data %>% filter(Subtype=="H3N2")) + 
  geom_boxplot(aes(x=DisplayRead, y=NumVariants/(BpSequenced/1000)),
               outlier.size=0.5, size=0.5) +
  xlab("Poon et al. 2016, H3N2") + 
  ylab("Number of within-host variants\nper kB sequenced") +
  THEME_ALL
p_NumVariants
save_plot(paste0(outdir, "CompareR1R2-NumVariants.pdf"),p_NumVariants,
          base_width=2.25, base_height=3.5)
save_plot(paste0(outdir, "CompareR1R2-NumVariants.jpg"),p_NumVariants,
          base_width=2.25, base_height=3.5)

# Read in variants identified in all datasets.
Data <- read.table("analysis/CompareR1R2/out/variants-all.txt",
                   header=TRUE, stringsAsFactors = FALSE) %>%
  mutate(StudySubtype=paste(Study,Subtype,sep="\n"),
         StudySample=paste(Study,Sample,sep="\n"))
# Format analysis names for display.
Reads <- c("Poon et al. 2016, H3N2, all reads", 
           "Poon et al. 2016, H3N2, read 1 only",
           "Poon et al. 2016, H3N2, read 2 only")
names(Reads) <- c("R1R2","R1","R2")
Data <- Data %>% mutate(DisplayRead=Reads[Read])
Data$DisplayRead <- factor(Data$DisplayRead,
                           levels=Reads)

# Output the number of genome sites in each study and subtype
# at which more than half of samples display within-host variation.
SharedVariation <- Data %>% group_by(StudySubtype, Read, NumSamplesSequenced, GenomePos) %>%
  summarize(NumVariableSamples=n_distinct(Sample)) %>% ungroup() %>%
  group_by(StudySubtype, Read, NumSamplesSequenced) %>%
  summarize(NumSharedVariableSites=sum(NumVariableSamples > (NumSamplesSequenced)/2)) %>%
  ungroup() %>% separate(StudySubtype, into=c("Study", "Subtype"), sep="\n")
write.table(SharedVariation, paste0(outdir, "SitesOfSharedVariation-R1R2.txt"),
            quote=FALSE, row.names=FALSE)

# Plot shared variation across the genome for each study.
p_Shared <- Data %>% filter(Subtype=="H3N2") %>%
  group_by(DisplayRead, NumSamplesSequenced, GenomePos, Chr) %>% 
  summarize(NumSamples=n_distinct(Sample)) %>%
  ggplot() + geom_bar(aes(x=GenomePos, y=NumSamples/NumSamplesSequenced,
                          color=factor(Chr)), stat="identity") +
  facet_wrap(~DisplayRead, ncol=1, scales="free") +
  ylim(0,1) + xlab("Genome position (bp)") + 
  ylab("Proportion of samples in study\nwith within-host variation") +
  scale_color_manual(values=rep(c(PALETTE[1],PALETTE[2]),4), name="Gene") +
  guides(color=FALSE) +
  scale_x_continuous(limits=c(0,13600), expand=c(0,0),
                     breaks=GenomeAverage$Midpoint,
                     labels=c("PB2","PB1","PA","HA","NP","NA","M","NS")) +
  scale_y_continuous(limits=c(0,1), breaks=c(0,0.5,1)) +
  THEME_ALL
p_Shared
save_plot(paste0(outdir,"CompareR1R2-SharedVariants.pdf"), p_Shared,
          base_width=7, base_height=3)
save_plot(paste0(outdirpres,"CompareR1R2-SharedVariants.pdf"), 
          p_Shared+THEME_PRES + theme(axis.text=element_text(size=7),
                                      strip.text=element_text(size=11, face="bold")),
          base_width=9, base_height=6)

p_Compare <- plot_grid(p_NumVariants, p_Shared, 
                       ncol=2, rel_widths=c(2,5), labels=c("A","B"),
                       label_size=8)
p_Compare
save_plot(paste0(outdir, "CompareR1R2-Diversity.pdf"), p_Compare,
          base_width=7, base_height=3)

