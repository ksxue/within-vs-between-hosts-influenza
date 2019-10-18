library(ggplot2)
library(tidyverse)
library(cowplot)
library(foreach)
library(Cairo)
library(scales)

options(bitmapType='cairo')

source("analysis/PlotThemes.R")
source("analysis/FileFormats.R")

# Import sample exclusions.
source("analysis/CalculateAcuteRates/IdentifySampleExclusions.R")

outdir <- "analysis/figures-WvB/"
outpres <- "analysis/figures-WvB/pres/"

# Render values in scientific notation.
# Function modified from 
# https://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales
scientific_10 <- function(x) {
  ifelse(x==0,0,parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x))))
}

# Function to rename Project variable (flu-Dinis, flu-Debbink, flu-McCrone)
# to display manuscript names instead.
Projects <- list("Dinis et al., 2016","Debbink et al., 2017","McCrone et al., 2018")
names(Projects) <- c("flu-Dinis","flu-Debbink","flu-McCrone")
ProjectOrder <- unlist(Projects)

# Function to rename Gene variable to remove numeric identifier.
Genes <- list("PB2","PB1","PA","HA","NP","NA","M1","M2","NS1","NEP")
names(Genes) <- FluGenes

# Mutation type ordering.
MutationTypes <- c("S","NS","Stop")
MUTATIONTYPESPALETTE <- c(DARK2PALETTE[2],DARK2PALETTE[1],DARK2PALETTE[3])

# Figure 2. Plot an overview of variant-calling thresholds --------------------------

# Import summary of number of variants per sample.
DataNumVariants <- read.table("analysis/CompareVariantCallingThresholds/out/NumVariantsPerSample.data",
                              header=TRUE, stringsAsFactors = FALSE)

# Plot a distribution of the number of variants per sample
# at different variant-calling thresholds.
# Plot for only H3N2 samples.
# Exclude samples from plasmid controls.
# Note that each variant represents a distinct genome position,
# so variants in overlapping regions of the genome are counted once.
p <- DataNumVariants %>%
  filter(!(Sample %in% ControlSamples$Sample),
         Subtype=="H3N2") %>%
  ggplot() +
  geom_boxplot(aes(x=(MinFreq), y=NumVariants, 
                   fill=factor(MinFreq==0.005), group=factor(MinFreq)), 
               size=0.5, outlier.size=0.1) +
  scale_y_log10() +
  scale_fill_manual(values=c("gray60","firebrick3")) +
  guides(fill=FALSE) +
  xlab("Mutation-frequency threshold") + ylab("Number of mutations per sample") +
  THEME_ALL
save_plot(paste0(outdir,"NumVariantsPerSample-ByMinFreq-H3N2.png"),p,
          base_height=3, base_width=3)
save_plot(paste0(outdir,"NumVariantsPerSample-ByMinFreq-H3N2.pdf"),p,
          base_height=3, base_width=3)

# Import variants above a frequency of 0.5%.
DataVariants <- read.table("analysis/CallWithinHostVariants/variants-annotated-0.005.data.gz",
                           header=TRUE, stringsAsFactors = FALSE)
# Plot a CDF of variant frequencies for variants above 0.5% frequency.
# Exclude samples according to pre-determined criteria.
p <- DataVariants %>%
  filter(Subtype=="H3N2", Project!="flu-Xue-chronic") %>%
  filter(!Sample %in% SampleExclusions$Sample) %>%
  ggplot() +
  stat_ecdf(aes(x=Freq), size=1) +
  xlab("Mutation frequency") + ylab("Proportion of mutations") +
  THEME_ALL
save_plot(paste0(outdir,"VariantFrequencyCDF-0.005-H3N2.png"),p,
          base_width=2, base_height=1.5)
save_plot(paste0(outdir,"VariantFrequencyCDF-0.005-H3N2.pdf"),p,
          base_width=2, base_height=1.5)

# Import proportion of samples at each codon position.
DataCodonPosition <- read.table("analysis/CompareVariantCallingThresholds/out/VariantCodonPositions.data",
                                header=TRUE, stringsAsFactors = FALSE)
# Plot the proportion of variants in each frequency bin
# at each codon position.
# Variants from the M1/M2 and NS1/NEP genes have been excluded,
# since the interpretation of variant position is not clear.
# Plot only variants from H3N2 samples, and require 200 variants in each bin.
p <- DataCodonPosition %>%
  filter(Subtype=="H3N2") %>%
  group_by(FreqBin) %>% filter(sum(NumVariants)>200) %>%
  mutate(Freq=10^FreqBin) %>%
  ggplot() +
  geom_vline(xintercept=0.005, linetype="dashed") +
  geom_line(aes(x=Freq, y=PercentVariants, color=factor(CodonPos)), size=1) +
  scale_color_brewer(palette="Dark2", name="Codon\nposition") +
  scale_x_log10(
    limits=c(1e-05,10^(-1.5)),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) + 
  xlab("Mutation frequency") +
  ylab("Proportion of mutations") + ylim(0,1) +
  theme(legend.position=c(0.75,0.8),
        legend.key.size=unit(0.25,"cm"),
        legend.text=element_text(size=6),
        legend.title = element_text(size=6)) +
  THEME_ALL
save_plot(paste0(outdir,"VariantProportion-CodonPosition-H3N2.png"),p,
          base_width=2, base_height=1.5)
save_plot(paste0(outdir,"VariantProportion-CodonPosition-H3N2.pdf"),p,
          base_width=2, base_height=1.5)


# Figure 3. Plot evolutionary rate comparisons --------------------------------------

# Import the aggregated list of all evolutionary rates.
DataRatesAll <- read.table("analysis/CalculateGlobalRates/out/CombinedRates.data",
                           header=TRUE, stringsAsFactors = FALSE)

# Add in global rates of stop codon evolution.
# These are set to zero rather than being calculated directly from the data
# due to issues of definitional circularity.
DataRatesAll <- 
  union(DataRatesAll,
        DataRatesAll %>%
  filter(Scale=="global") %>% dplyr::select(Subtype,MinFreq,Gene,Type,Scale) %>%
  unique() %>%
  mutate(Syn="Stop", MeanDivPerSitePerDay=0, SEDivPerSitePerDay=0))

# Compare acute and global rates of evolution
# using different variant-call thresholds for the acute point estimate.
p <- DataRatesAll %>%
  filter(Type %in% c("point-all-0.005","point-all-0.01", "point-all-0.02","global")) %>%
  filter(Gene %in% FluGenes) %>%
  ggplot() +
  geom_errorbar(aes(x=Gene, ymin=MeanDivPerSitePerDay-SEDivPerSitePerDay,
                    ymax=MeanDivPerSitePerDay+SEDivPerSitePerDay, color=factor(Syn)),
                width=0.3) +
  geom_point(aes(x=Gene, y=MeanDivPerSitePerDay, color=factor(Syn),
                 shape=factor(Scale))) +
  facet_wrap(~MinFreq, ncol=4) +
  scale_color_brewer(palette="Dark2", name="Mutation\ntype") +
  scale_shape_manual(values=c(1,16), name="Scale") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  ylab("divergence / site / day") +
  coord_cartesian(ylim=c(0,4e-05)) +
  THEME_ALL
save_plot(paste0(outdir,"CompareAcuteGlobalRates-ByMinFreq.png"), p, ncol=2)  

# Compare rates of acute and global rates of evolution
# using different variant-call thresholds for the acute point estimate
# at the scale of each gene.
p <- DataRatesAll %>%
  filter(Type %in% c("point-all-0.005","point-all-0.01", "point-all-0.02","global")) %>%
  filter(Gene %in% FluGenes, !(Gene %in% c("7-M2", "8-NEP"))) %>%
  ggplot() +
  geom_errorbar(aes(x=MinFreq, ymin=MeanDivPerSitePerDay-SEDivPerSitePerDay,
                    ymax=MeanDivPerSitePerDay+SEDivPerSitePerDay, color=factor(Syn)),
                width=0.3) +
  geom_point(aes(x=MinFreq, y=MeanDivPerSitePerDay, color=factor(Syn),
                 shape=factor(Scale))) +
  facet_wrap(~Gene, ncol=4, scales="free_x") +
  scale_color_brewer(palette="Dark2", name="Mutation\ntype") +
  scale_shape_manual(values=c(1,16), name="Scale") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  xlab("Mutation-frequency threshold") +
  ylab("divergence / site / day") +
  coord_cartesian(ylim=c(0,4e-05)) +
  THEME_ALL
save_plot(paste0(outdir,"CompareAcuteGlobalRates-ByMinFreqByGene.png"),p, ncol=2)

# Compare rates of S and NS evolution within hosts.
DataWithinRatesRatios <-
  DataRatesAll %>%
  filter(Type %in% c("point-all-0.005")) %>%
  filter(Gene %in% FluGenesNonOverlapping) %>%
  dplyr::select(-SEDivPerSitePerDay) %>%
  spread(Syn,MeanDivPerSitePerDay) %>%
  mutate(StoNS=S/NS, StoStop=S/Stop)

# Compare acute and global rates of evolution
# using a single variant-call threshold for the acute point estimate.
p <- DataRatesAll %>%
  filter(Type %in% c("point-all-0.005","global")) %>%
  filter(Gene %in% FluGenesNonOverlapping) %>%
  mutate(Scale=ifelse(Scale=="acute","Within-host","Global")) %>%
  mutate(Gene=unlist(Genes[Gene]), Gene=fct_relevel(Gene,unlist(Genes))) %>%
  ggplot() +
  geom_errorbar(aes(x=Gene, ymin=MeanDivPerSitePerDay-SEDivPerSitePerDay,
                    ymax=MeanDivPerSitePerDay+SEDivPerSitePerDay, 
                    color=fct_relevel(Syn, MutationTypes)),
                width=0.1) +
  geom_point(aes(x=Gene, y=MeanDivPerSitePerDay, 
                 color=fct_relevel(Syn, MutationTypes),
                 shape=fct_relevel(Syn, MutationTypes)), size=1.8) +
  facet_wrap(~fct_rev(Scale), ncol=2, scales="free") +
  scale_color_manual(values=MUTATIONTYPESPALETTE, name="Mutation\ntype") +
  scale_shape_manual(values=SHAPES, name="Mutation\ntype") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  ylab("Divergence / site / day") +
  coord_cartesian(ylim=c(0,2.3e-05)) +
  scale_y_continuous(labels=scientific_10) +
  THEME_ALL
save_plot(paste0(outdir,"CompareAcuteGlobalRates.png"), p,
          base_height=2.5, base_width=5)
save_plot(paste0(outdir,"CompareAcuteGlobalRates.pdf"), p,
          base_height=2.5, base_width=5)
save_plot(paste0(outdir,"CompareAcuteGlobalRates-resize.pdf"), p,
          base_height=3, base_width=4.5)

# Figure 4. Plot evolutionary rates in HA and NA special sites. ------------------

# Import the aggregated list of all evolutionary rates.
DataRatesAll <- read.table("analysis/CalculateGlobalRates/out/CombinedRates.data",
                           header=TRUE, stringsAsFactors = FALSE)
# Analyze only the point-estimate rates calculated at a 0.5% threshold
# across all studies and in the global flu population.
DataRatesAll <- DataRatesAll %>%
  filter(Type %in% c("point-all-0.005","global"))
# Add in global rates of stop codon evolution.
# These are set to zero rather than being calculated directly from the data
# due to issues of definitional circularity.
DataRatesAll <-
  union(DataRatesAll,
        DataRatesAll %>%
          filter(Scale=="global") %>% dplyr::select(Subtype,MinFreq,Gene,Type,Scale) %>%
          unique() %>%
          mutate(Syn="Stop", MeanDivPerSitePerDay=0, SEDivPerSitePerDay=0))

# Plot all rates calculated in the HA gene.
p <- DataRatesAll %>%
  filter(Gene %in% c("4-HA-antigenic-Wolf","4-HA-nonantigenic-Wolf")) %>%
  mutate(Gene=ifelse(Gene=="4-HA-antigenic-Wolf","antigenic","non-antigenic")) %>%
  mutate(Type=ifelse(Type=="point-all-0.005","Within-host","Global")) %>%
  ggplot() +
  geom_errorbar(aes(x=Gene, ymin=MeanDivPerSitePerDay-SEDivPerSitePerDay,
                    ymax=MeanDivPerSitePerDay+SEDivPerSitePerDay, 
                    color=fct_relevel(Syn, MutationTypes)),
                width=0.1) +
  geom_point(aes(x=Gene, y=MeanDivPerSitePerDay, 
                 color=fct_relevel(Syn, MutationTypes),
                 shape=fct_relevel(Syn, MutationTypes)), size=1.8) +
  scale_color_manual(values=MUTATIONTYPESPALETTE, name="Mutation\ntype") +
  scale_shape_manual(values=SHAPES, name="Mutation\ntype") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  coord_cartesian(ylim=c(0e-05, 2.3e-05)) +
  scale_y_continuous(labels=scientific_10) +
  facet_wrap(~fct_rev(Type)) +
  ylab("Divergence / site / day") + xlab("Sites") +
  THEME_ALL
save_plot(paste0(outdir,"Figure4-AcuteRates-AntigenicSites.png"), p,
          base_height=3, base_width=3.5)
save_plot(paste0(outdir,"Figure4-AcuteRates-AntigenicSites.pdf"), p,
          base_height=3, base_width=3.5)

# Plot all rates calculated in the NA gene.
p <- DataRatesAll %>%
  filter(Gene %in% c("6-NA-surface-Bhatt","6-NA-nonsurface-Bhatt")) %>%
  mutate(Gene=ifelse(Gene=="6-NA-surface-Bhatt","exposed","buried")) %>%
  mutate(Type=ifelse(Type=="point-all-0.005","Within-host","Global")) %>%
  ggplot() +
  geom_errorbar(aes(x=fct_rev(Gene), ymin=MeanDivPerSitePerDay-SEDivPerSitePerDay,
                    ymax=MeanDivPerSitePerDay+SEDivPerSitePerDay, color=factor(Syn)),
                width=0.1) +
  geom_point(aes(x=fct_rev(Gene), y=MeanDivPerSitePerDay, color=factor(Syn)), size=1.8) +
  scale_color_brewer(palette="Dark2", name="Mutation\ntype") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  coord_cartesian(ylim=c(0e-05, 2.3e-05)) +
  scale_y_continuous(labels=scientific_10) +
  facet_wrap(~fct_rev(Type)) +
  ylab("Divergence / site / day") + xlab("Sites") +
  THEME_ALL
# save_plot(paste0(outdir,"Figure4-AcuteRates-SurfaceSites.png"), p,
#          base_height=3, base_width=3.5)
# save_plot(paste0(outdir,"Figure4-AcuteRates-SurfaceSites.pdf"), p,
#          base_height=3, base_width=3.5)


# Figure 5. Plot SFS comparisons ----------------------------------------

# Import NS proportions within and between hosts.
DataNSProportions <- read.table("analysis/InferGlobalVariants/out/H3N2-unpassaged-NSproportions.data",
                                header=TRUE, stringsAsFactors = FALSE)

# Import proportion of available NS sites.
AvailableSites <- read.table("analysis/CalculateAcuteRates/out/AvailableSitesByType.data",
                             header=TRUE, stringsAsFactors = FALSE)
# Classify stop sites as NS sites.
AvailableNS <- (AvailableSites %>% 
                  spread(Syn,PercentSites) %>%
                  mutate(NS=NS+Stop))$NS[1]

# Calculate bootstrap intervals for display.
# When there are no NS and S variants (i.e. ProportionNS cannot be calculated),
# set the NS proportion to zero.
DataNSProportionsAvg <- 
  DataNSProportions %>%
  mutate(Type=ifelse(Year=="within-host","Within-host","Global")) %>%
  group_by(Type, Gene, BinFreq) %>%
  replace_na(list(ProportionNS=0)) %>%
  # Calculate within-host bootstrap intervals.
  # Also calculate minimum and maximum between-host NS proportions
  # from the three seasons in which values were calculated.
  # Also calculate the mean of the three seasons.
  summarize(MeanProportionNS=mean(ProportionNS),
            Lower=quantile(ProportionNS,0.025, type=3),
            Upper=quantile(ProportionNS,0.975), type=3) %>%
  filter(Gene %in% FluGenesNonOverlapping,
         BinFreq<0.3, !(Type=="Within-host" & BinFreq==0.1)) %>%
  ungroup()

PlotNSProportions <- function(gene){
  p <- DataNSProportionsAvg %>% filter(Gene==gene) %>%
    mutate(Gene=unlist(Genes[Gene]), Gene=fct_relevel(Gene,unlist(Genes))) %>%
    ggplot() +
    geom_ribbon(aes(x=BinFreq, ymin=Lower,
                    ymax=Upper, fill=factor(fct_rev(Type)))) +
    geom_point(aes(x=BinFreq, y=MeanProportionNS)) +
    geom_line(aes(x=BinFreq, y=MeanProportionNS)) +
    geom_hline(aes(yintercept=AvailableNS), linetype="dashed") +
    facet_wrap(~fct_rev(Type), ncol=2, scales="free_x") +
    scale_linetype_discrete(name="Type") +
    scale_fill_manual(values=c("firebrick3", #red
                                "dodgerblue3")) + #light blue
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x, n=3),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) + ylim(0,1) +
    xlab("Mutation frequency") + ylab("Fraction of mutations\nthat are nonsynonymous") +
    guides(fill=FALSE) +
    ggtitle(unlist(Genes[gene])) +
    THEME_ALL + theme(plot.title=element_text(size=7),
                      plot.margin=unit(c(2,2,2,2),"mm"))
  return(p)  
}
p <- plot_grid(plotlist=lapply(FluGenesNonOverlapping, PlotNSProportions), ncol=2)
save_plot(paste0(outdir,"WvB-SFScomparison-NSproportions.png"),p,
          base_width=7,base_height=5.5)  
save_plot(paste0(outdir,"WvB-SFScomparison-NSproportions.pdf"),p,
          base_width=7,base_height=5.5)  



# Supplementary figure. Plot chronic evolutionary rates -----------------------------------------

# Import chronic evolutionary rates
# based on linear regression across all patients.
DataRates <- read.table("analysis/CalculateGlobalRates/out/CombinedRates.data",
                               header=TRUE, stringsAsFactors = FALSE)


# Import data on evolutionary rates in chronic infections,
# calculated separately for each patient.
DataChronicRatesByPatient <- DataRates %>%
  filter(Scale=="chronic", !grepl("chronic",Type))
# Plot the evolutionary rates in chronic infections.
# Plot evolutionary rates separately for each patient. Facet by patient.
p_ChronicByPatient <- DataChronicRatesByPatient %>% 
  filter(MinFreq==0.005) %>%
  filter(Gene %in% FluGenesNonOverlapping) %>%
  mutate(Gene=unlist(Genes[Gene]), Gene=fct_relevel(Gene,unlist(Genes))) %>%
  mutate(Patient=substr(Type,1,1)) %>%
  ggplot() +
  geom_errorbar(aes(x=Gene, ymin=MeanDivPerSitePerDay-SEDivPerSitePerDay,
                    ymax=MeanDivPerSitePerDay+SEDivPerSitePerDay, 
                    color=fct_relevel(Syn, MutationTypes)),
                width=0.1, alpha=0.6) +
  geom_point(aes(x=Gene, y=MeanDivPerSitePerDay, 
                 color=fct_relevel(Syn, MutationTypes),
                 shape=fct_relevel(Syn, MutationTypes)), size=1.8, alpha=0.6) +
  scale_color_manual(values=MUTATIONTYPESPALETTE, name="Mutation\ntype") +
  scale_shape_manual(values=SHAPES, name="Mutation\ntype") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  facet_wrap(~Patient, ncol=2, scales="free") + 
  coord_cartesian(ylim=c(-0.1e-05, 9e-05)) +
  scale_y_continuous(labels=scientific_10) +
  ylab("Divergence / site / day") +
  THEME_ALL
save_plot(paste0(outdir,"ChronicRates-ByPatient-0.005.png"),p_ChronicByPatient, 
          base_width=3, base_height=3)
save_plot(paste0(outdir,"ChronicRates-ByPatient-0.005.pdf"),p_ChronicByPatient, 
          base_width=3, base_height=3)

# Average the evolutionary rates across patients
# (note that here, the rates calculated for each patient through regression are being averaged).
# Plot a single rate for each gene and site class.
p_CompareAcuteChronic <- DataRates %>% 
  filter(((Scale=="chronic" & grepl("chronic",Type)) | 
            (Scale=="acute" & grepl("point-all",Type)))) %>%
  filter(MinFreq==0.005, Gene %in% FluGenesNonOverlapping) %>% ungroup() %>%
  mutate(Gene=unlist(Genes[Gene]), Gene=fct_relevel(Gene,unlist(Genes))) %>%
  mutate(Type=ifelse(grepl("chronic",Type),"Chronic","Acute")) %>%
  ggplot() +
  geom_errorbar(aes(x=Gene, ymin=MeanDivPerSitePerDay-SEDivPerSitePerDay,
                    ymax=MeanDivPerSitePerDay+SEDivPerSitePerDay, 
                    color=fct_relevel(Syn, MutationTypes)),
                width=0.1, alpha=0.6) +
  geom_point(aes(x=Gene, y=MeanDivPerSitePerDay, 
                 color=fct_relevel(Syn, MutationTypes),
                 shape=fct_relevel(Syn, MutationTypes)), size=1.8,alpha=0.6) +
  scale_color_manual(values=MUTATIONTYPESPALETTE, name="Mutation\ntype") +
  scale_shape_manual(values=SHAPES, name="Mutation\ntype") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  coord_cartesian(ylim=c(-0.1e-05, 4.3e-05)) +
  scale_y_continuous(labels=scientific_10) +
  facet_wrap(~Type, scales="free") +
  ylab("Divergence / site / day") +
  guides(color=FALSE, shape=FALSE) +
  THEME_ALL 
save_plot(paste0(outdir,"ChronicRates-AllPatients-0.005.png"),p_CompareAcuteChronic,
          base_width=3.5, base_height=3)
save_plot(paste0(outdir,"ChronicRates-AllPatients-0.005.pdf"),p_CompareAcuteChronic,
          base_width=3.5, base_height=3)

# Assemble plot components.
p <- plot_grid(p_CompareAcuteChronic, p_ChronicByPatient, 
          labels=c("A","B"), label_size=11, rel_widths=c(1,1.2))
save_plot(paste0(outdir,"FigureS6-ChronicRates.pdf"), p,
          base_width=7, base_height=2.8)
save_plot(paste0(outdir,"FigureS6-ChronicRates.png"), p,
          base_width=7, base_height=2.8)

# Supplemental figure. Plot sample DPI metadata. -----------------------------

# Import parsed sample DPI metadata for the Debbink and McCrone studies.
source("analysis/CalculateAcuteRates/ParseSampleMetadata.R")
# Import list of sample exclusions.
SampleExclusions <- read.table("analysis/CalculateAcuteRates/out/SampleExclusions.data",
                               header=TRUE, stringsAsFactors = FALSE)
DataDPI <- rbind(
  MetadataDPI %>%
    filter(!(Sample %in% SampleExclusions$Sample),
         Subtype=="H3N2") %>%
    group_by(Project,DPI) %>% summarize(NumPatients=n()) %>%
    ungroup(),
  read.table("data/metadata/flu-Dinis/metadata-H3N2-2012-2013.txt",
           header=TRUE, stringsAsFactors = FALSE) %>%
    mutate(Project="flu-Dinis") %>%
    dplyr::select(Project,DPI,NumPatients))

# Plot distribution of sample DPIs for each study.
p <- DataDPI %>%
  mutate(Project=unlist(Projects[Project])) %>%
  ggplot() +
  geom_bar(aes(x=DPI, y=NumPatients), stat="identity") +
  facet_wrap(~fct_relevel(Project,ProjectOrder), ncol=3) +
  xlab("Days post-symptom-onset") + ylab("Number of samples") +
  scale_x_continuous(limits=c(-1.5,8), breaks=seq(0,8,2)) +
  THEME_ALL
save_plot(paste0(outdir,"MetadataDPIs.pdf"),p,
          base_width=3.5, base_height=1.8)
save_plot(paste0(outdir,"FigureS1-MetadataDPIs.pdf"),p,
          base_width=3.5, base_height=1.8)
save_plot(paste0(outdir,"FigureS1-MetadataDPIs.png"),p,
          base_width=3.5, base_height=1.8)

# Supplemental figure. Plot within-host rates by min freq -------------------------------------

# Read in average evolutionary rates in acute infections,
# averaged across all studies.
DataWithinRates <- read.table("analysis/CalculateAcuteRates/out/AcuteRatesAllStudies.data",
                              header=TRUE, stringsAsFactors = FALSE)

# Plot the average evolutionary rates in acute infections for each influenza gene
# at different variant-frequency thresholds.
p <- DataWithinRates %>%
  filter(Gene %in% FluGenesNonOverlapping) %>%
  mutate(Gene=unlist(Genes[Gene]), Gene=fct_relevel(Gene,unlist(Genes))) %>%
  ggplot() +
  geom_errorbar(aes(x=factor(MinFreq), ymin=MeanDivPerSitePerDay-SEDivPerSitePerDay,
                    ymax=MeanDivPerSitePerDay+SEDivPerSitePerDay, 
                    color=fct_relevel(Syn, MutationTypes)),
                width=0.1) +
  geom_point(aes(x=factor(MinFreq), y=MeanDivPerSitePerDay, 
                 color=fct_relevel(Syn, MutationTypes),
                 shape=fct_relevel(Syn, MutationTypes)), size=1.8) +
  geom_line(aes(x=factor(MinFreq), y=MeanDivPerSitePerDay, 
                color=fct_relevel(Syn, MutationTypes),
                group=factor(Syn))) +
  scale_color_manual(values=MUTATIONTYPESPALETTE, name="Mutation\ntype") +
  scale_shape_manual(values=SHAPES, name="Mutation\ntype") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  facet_wrap(~Gene, ncol=6) +
  coord_cartesian(ylim=c(0,2.3e-05)) +
  scale_y_continuous(labels=scientific_10) +
  xlab("Mutation-frequency threshold") +
  ylab("Divergence / site / day") +
  THEME_ALL
save_plot(paste0(outdir,"AcuteRates-ByMinFreq.png"),p,
          base_height=1.8, base_width=7)
save_plot(paste0(outdir,"FigureS2-AcuteRates-ByMinFreq.pdf"),p,
          base_height=1.8, base_width=7)
save_plot(paste0(outdir,"FigureS2-AcuteRates-ByMinFreq.png"),p,
          base_height=1.8, base_width=7)


# Supplemental figure. Plot within-host rates by study. -------------------

# Read in average evolutionary rates in each gene,
# calculated separately for each study.
DataWithinRatesByStudy <- 
  read.table("analysis/CalculateAcuteRates/out/AcuteRatesByStudy.data",
             header=TRUE, stringsAsFactors = FALSE)

# Plot the average evolutionary rates in acute infections for each influenza gene,
# calculated separately for each study.
# Facet by study.
p <- DataWithinRatesByStudy %>% filter(MinFreq==0.005) %>%
  filter(Gene %in% FluGenesNonOverlapping) %>%
  mutate(Project=unlist(Projects[Project])) %>%
  mutate(Gene=unlist(Genes[Gene]), Gene=fct_relevel(Gene,unlist(Genes))) %>%
  ggplot() +
  geom_errorbar(aes(x=Gene, ymin=MeanDivPerSitePerDay-SEDivPerSitePerDay,
                    ymax=MeanDivPerSitePerDay+SEDivPerSitePerDay, 
                    color=fct_relevel(Syn, MutationTypes)),
                width=0.1) +
  geom_point(aes(x=Gene, y=MeanDivPerSitePerDay, 
                 color=factor(Syn, MutationTypes),
                 shape=factor(Syn, MutationTypes)), size=1.8) +
  scale_color_manual(values=MUTATIONTYPESPALETTE, name="Mutation\ntype") +
  scale_shape_manual(values=SHAPES, name="Mutation\ntype") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  facet_wrap(~fct_relevel(Project,ProjectOrder),
             ncol=3, scales="free_y") +
  ylab("Divergence / site / day") +
  scale_y_continuous(labels=scientific_10,lim=c(0,3.2e-05)) +
  THEME_ALL
save_plot(paste0(outdir,"AcuteRates-ByStudy-0.005.png"),p,
          base_width=7, base_height=2.5)
save_plot(paste0(outdir,"AcuteRates-ByStudy-0.005.pdf"),p,
          base_width=7, base_height=2.5)
save_plot(paste0(outdir,"FigureS4-AcuteRates-ByStudy-0.005.pdf"),p,
          base_width=7, base_height=2.5)
save_plot(paste0(outdir,"FigureS4-AcuteRates-ByStudy-0.005.png"),p,
          base_width=7, base_height=2.5)

# Supplemental figure. Within-host rates by regression. -------------------

# Read in average evolutionary rates in acute infections,
# calculated through linear regression across all studies.
DataWithinRatesRegression <- 
  read.table("analysis/CalculateAcuteRates/out/AcuteRatesAllStudies-Regression.data",
             header=TRUE, stringsAsFactors = FALSE)

# Plot a comparison of evolutionary rates in acute infections
# as calculated through the point estimate and regression methods.
p <- rbind(DataWithinRates %>% mutate(Type="point"),
           DataWithinRatesRegression %>% mutate(Type="regression")) %>%
  filter(MinFreq==0.005, Gene %in% FluGenesNonOverlapping) %>%
  mutate(Gene=unlist(Genes[Gene]), Gene=fct_relevel(Gene,unlist(Genes))) %>%
  ggplot() +
  geom_errorbar(aes(x=Gene, ymin=MeanDivPerSitePerDay-SEDivPerSitePerDay,
                    ymax=MeanDivPerSitePerDay+SEDivPerSitePerDay, 
                    color=fct_relevel(Syn, MutationTypes)),
                width=0.1) +
  geom_point(aes(x=Gene, y=MeanDivPerSitePerDay, 
                 color=fct_relevel(Syn, MutationTypes),
                 shape=fct_relevel(Syn, MutationTypes)), size=1.8) +
  facet_wrap(~Type, ncol=2) +
  scale_color_manual(values=MUTATIONTYPESPALETTE, name="Mutation\ntype") +
  scale_shape_manual(values=SHAPES, name="Mutation\ntype") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  ylab("Divergence / site / day") +
  coord_cartesian(ylim=c(0,2.3e-05)) +
  scale_y_continuous(labels=scientific_10) +
  THEME_ALL
save_plot(paste0(outdir,"AcuteRates-CompareRegression-0.005.png"),p,
          base_width=3.5, base_height=2)  
save_plot(paste0(outdir,"AcuteRates-CompareRegression-0.005.pdf"),p,
          base_width=3.5, base_height=2)  
save_plot(paste0(outdir,"FigureS5-AcuteRates-CompareRegression-0.005.pdf"),p,
          base_width=3.5, base_height=2)  
save_plot(paste0(outdir,"FigureS5-AcuteRates-CompareRegression-0.005.png"),p,
          base_width=3.5, base_height=2)



# Supplemental figure. Plot global evolutionary rates ------------------------------------------

# Import global sequence distances from reference sequence,
# after removing outliers.
DataGlobalSequenceDistances <-
  rbind(read.table(paste0("analysis/CalculateGlobalRates/out/1999-SequenceDistances.data"),
                   header=TRUE, stringsAsFactors = FALSE) %>% mutate(Interval=1999),
        read.table(paste0("analysis/CalculateGlobalRates/out/2007-SequenceDistances.data"),
                   header=TRUE, stringsAsFactors = FALSE) %>% mutate(Interval=2007))

# Plot the distance of each sequence from the reference over time.
p <- DataGlobalSequenceDistances %>%
  filter(Gene %in% FluGenesNonOverlapping) %>%
  filter(Interval==ifelse(Gene %in% c("4-HA","6-NA","8-NS1"),2007,1999)) %>%
  filter(CollectionDate>Interval,
         CollectionDate<(ifelse(Interval==1999,Interval+20,Interval+10))) %>%
  mutate(Gene=unlist(Genes[Gene]), Gene=fct_relevel(Gene,unlist(Genes))) %>%
  ggplot() +
  geom_point(aes(x=CollectionDate, y=Distance, 
                 color=fct_relevel(Syn, MutationTypes),
                 shape=fct_relevel(Syn, MutationTypes)), alpha=0.2) +
  facet_wrap(~Gene, scales="free_y", ncol=6) +
  scale_color_manual(values=MUTATIONTYPESPALETTE, name="Mutation\ntype") +
  scale_shape_manual(values=SHAPES, name="Mutation\ntype") +
  scale_x_continuous(breaks=pretty_breaks(3)) + ylim(0,80) +
  xlab("Year") + ylab("Divergence from\nreference sequence") +
  THEME_ALL + theme(axis.text=element_text(size=4.5))
save_plot(paste0(outdir,"FigureS7-GlobalSequenceDistances.pdf"),p,
          base_width=7, base_height=1.8)
save_plot(paste0(outdir,"FigureS7-GlobalSequenceDistances.png"),p,
          base_width=7, base_height=1.8)


# Supplemental figure. Plot sample exclusions. ----------------------------


# Import variants called at a frequency of at least 0.5%.
DataWithin <- read.table("analysis/CallWithinHostVariants/variants-annotated-0.005.data.gz",
                         header=TRUE, stringsAsFactors = FALSE) %>%
  filter(Project %in% c("flu-Dinis","flu-Debbink","flu-McCrone"))

# Summarize the number of variants identified per sample.
DataWithinNumVariants <- DataWithin %>% filter(Subtype=="H3N2") %>%
  group_by(Project, Sample) %>% summarize(NumVariants=n())

# Ensure that samples with zero variants are included in this analysis.
DataWithinNumVariants <-
  rbind(DataWithinNumVariants %>% ungroup(),
        Metadata %>% dplyr::select(Project, Subtype, Sample) %>%
          filter(Project %in% c("flu-Dinis","flu-Debbink","flu-McCrone"),
                 Subtype=="H3N2", !(Sample %in% DataWithinNumVariants$Sample)) %>%
          dplyr::select(-Subtype) %>% mutate(NumVariants=0))

# Import annotated list of sample exclusions.
AnnotatedSampleExclusions <- 
  read.table("analysis/CalculateAcuteRates/out/SampleExclusions-annotated.data",
             header=TRUE, stringsAsFactors = FALSE) %>%
  filter(Project!="flu-Xue-chronic")
# For the three samples that are excluded based on two criteria,
# select a single criterion to display the reason for exclusion.
# In this case, value TopNumbersOfVariants as a reason above the others.
AnnotatedSampleExclusions <- AnnotatedSampleExclusions %>%
  group_by(Sample) %>% dplyr::arrange(desc(Reason)) %>% top_n(1, Reason)

# Add the sample exclusion information to the data frame listing numbers of variants.
DataWithinNumVariants <- 
  left_join(DataWithinNumVariants,
            AnnotatedSampleExclusions,
            by=c("Project","Sample")) %>%
  replace_na(list(Reason="None"))

# Calculate the total number of samples per project and the number of samples analyzed
# after applying the sample exclusion criteria.
NumSamples <- DataWithinNumVariants %>%
  group_by(Project) %>% 
  summarize(SamplesAnalyzed=sum(Reason=="None"),
            SamplesTotal=n()) %>%
  mutate(RegionSequenced=ifelse(Project=="flu-Dinis","HA only","full genome"))

# Function to rename Project variable (flu-Dinis, flu-Debbink, flu-McCrone)
# to display manuscript names instead.
Projects <- list("Dinis et al., 2016","Debbink et al., 2017","McCrone et al., 2018")
names(Projects) <- c("flu-Dinis","flu-Debbink","flu-McCrone")
ProjectOrder <- unlist(Projects)
# Function to rename Reason variable to display a more expanded explanation.
Reasons <- list("Top 10% of samples,\ntotal number of variants", "Low sequencing\ncoverage",
                "Longitudinal pair,\nfirst timepoint", "Plasmid control", "None")
names(Reasons) <- c("TopNumbersOfVariants","LowCoverage","LongitudinalFirstTimepoint",
                    "PlasmidControl","None")
ReasonsOrder <- unlist(Reasons)

# Plot sample exclusions using a dotplot.
p <- DataWithinNumVariants %>%
  mutate(Project=unlist(Projects[Project])) %>%
  mutate(Reason=unlist(Reasons[Reason])) %>%
  ggplot() +
  geom_dotplot(aes(x=NumVariants, 
                   fill=fct_relevel(Reason, ReasonsOrder), 
                   color=fct_relevel(Reason, ReasonsOrder)),
               stackgroups=TRUE, method="histodot", dotsize=1) +
  geom_text(data=NumSamples %>% mutate(Project=unlist(Projects[Project])), 
            aes(x=25, y=0.97, 
                label=paste0("total samples: ", SamplesTotal,
                             "\nsamples analyzed: ", SamplesAnalyzed,
                             "\nregion sequenced: ", RegionSequenced)),
            size=2) +
  facet_wrap(~fct_relevel(Project,ProjectOrder)) + scale_x_log10() +
  scale_color_manual(values=c("#7570B3","#E6AB02","#66A61E","#E7298A","grey30"),
                     name="Reason for\nsample exclusion") +
  scale_fill_manual(values=c("#7570B3","#E6AB02","#66A61E","#E7298A","grey30"),
                    name="Reason for\nsample exclusion") +
  xlab("Number of mutations") + 
  theme(axis.line.y=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  THEME_ALL
save_plot(paste0(outdir,"FigureS3-SampleExclusions.pdf"),p,base_width=7, base_height=3)
save_plot(paste0(outdir,"FigureS3-SampleExclusions.png"),p,base_width=7, base_height=3)
