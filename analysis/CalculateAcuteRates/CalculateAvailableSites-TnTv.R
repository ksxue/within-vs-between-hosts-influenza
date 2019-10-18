library(ggplot2)
library(tidyverse)
library(cowplot)
library(foreach)
library(Cairo)
library(Biostrings)
library(seqinr)

# Output directory for generated data.
outdir<-"analysis/CalculateAcuteRates/out/"

# Transition-transversion ratio: Tn/Tv
# Bloom, MBE 2014 estimates this ratio as 3.
# Pauly et al., eLife 2017 estimates this ratio as 2.7-3.6.
TnTv <- 3
# There are three possible mutations at a site,
# one of which is a transition and two of which are transversions.
# Given the transition-transversion ratio given above,
# assign weights to the transition and transversion mutations
# such that the total mutational "weight" among the three mutations is 3,
# but the transition has more weight proportional to its
# increased likelihood of occurring.
weightTn <- 3*TnTv/(1+TnTv)
weightTv <- (3-weightTn)/2

# Write a function to parse a reference coding sequence
# and count the available sites for stop codons.
CountStopCodonSites <- function(refseq){
  refseq <- s2c(refseq)
  # Split the reference sequence into a list of codons.
  refcodons <- sapply(seq(1,length(refseq)-1,3), function(x) c2s(refseq[x:(x+2)]))
  # Score each codon to see how many stop codons can be produced,
  # and sum the results for the entire sequence.
  return(sum(sapply(refcodons, ScoreCodonForStops)))
}
# Write a function to parse a reference coding sequence
# and count the available sites for synonymous mutations.
CountSMutationSites <- function(refseq){
  refseq <- s2c(refseq)
  # Split the reference sequence into a list of codons.
  refcodons <- sapply(seq(1,length(refseq)-1,3), function(x) c2s(refseq[x:(x+2)]))
  # Score each codon to see how many synonymous mutations can be produced,
  # and sum the results for the entire sequence.
  return(sum(sapply(refcodons, ScoreCodonForSMutations)))
}

# Score each codon to see how many stop codons can be produced.
# Treat each type of mutation equivalently.
ScoreCodonForStops <- function(codon){
  return(sum(sapply(seq(1:3), 
                    function(x) ScoreCodonPosition(codon, x, "*"))))
}
# Score each codon to see how many S mutations can be produced.
# Treat each type of mutation equivalently.
ScoreCodonForSMutations <- function(codon){
  return(sum(sapply(seq(1:3), 
                    function(x) ScoreCodonPosition(codon, x, translate(s2c(codon))))))
}
# Score each codon position for possible mutations that
# result in the codonAA amino acid.
ScoreCodonPosition <- function(codon, position, codonAA){
  codonchars <- s2c(codon)
  bases <- c("A","C","G","T")
  # for all bases that are not the reference base
  refbase <- codonchars[position]
  return(sum(sapply(bases[which(bases!=refbase)],
                    function(x)
                    {
                      testcodonchars <- codonchars
                      # replace the position of interest with an alternate base
                      testcodonchars[position] <- x
                      # determine whether the replacement was a transition or transversion
                      # and assign the appropriate corresponding weight
                      weight <- weightTv
                      if(paste0(refbase,x) %in% c("AG","GA","TC","CT")){
                        weight <- weightTn 
                      }
                      # return 1/3 if the base is changed to the base of interest,
                      # and 0 otherwise
                      # weight by whether the mutation is a transition or transversion
                      if(translate(testcodonchars) == codonAA){
                        return(weight*1/3)
                      }
                      else{
                        return(0)
                      }
                    })))
}


# Calculate the number of positions in a reference genome
# that, when mutated, can create a stop codon.
# These are the available sites for stop codon mutations.
# Each base can be mutated to three other bases,
# and those bases each count for 1/3 of an available site.
# Read in reference FASTA file containing coding sequences
# for each gene.
fastaFile <- readDNAStringSet("reference/flu-H3N2/H3N2-Victoria-2011-coding.fasta")
RefSeqs <- data.frame(names(fastaFile), paste(fastaFile)) %>%
  `colnames<-`(c("Header","Sequence"))
rm(fastaFile)
RefSeqs$Sequence <- as.character(RefSeqs$Sequence)
RefSeqs <- RefSeqs %>% separate(Header, into=c("Chr"), sep=" ")
# Remap chromosome names onto gene names.
RefSeqs <- RefSeqs %>% 
  mutate(Gene=ifelse(Chr=="7-M","7-M1",
              ifelse(Chr=="8-NS","8-NS1",Chr)))

# Calculate the number of S and stop codon sites for each sequence.
RefSeqs$StopSites <- sapply(RefSeqs$Sequence, CountStopCodonSites)
RefSeqs$SSites <- sapply(RefSeqs$Sequence, CountSMutationSites)
RefSeqs <- RefSeqs %>%
  mutate(NSSites=nchar(Sequence)-SSites,
         Length=nchar(Sequence),
         NSProportion=NSSites/Length,
         StopProportion=StopSites/Length)

# Average the proportion of S, NS, and stop sites between genes.
SiteProportions <- RefSeqs %>%
  ungroup() %>% 
  summarize(NSProportion=mean(NSProportion),
            StopProportion=mean(StopProportion)) %>%
  mutate(SProportion=1-NSProportion-StopProportion) %>%
  gather(key="Syn", value="PercentSites") %>%
  mutate(Syn=gsub("Proportion","",Syn))

# Export the average genome-wide proportions of each site.
write.table(SiteProportions, paste0(outdir,"AvailableSitesByType.data"),
            quote=FALSE, row.names=FALSE)
