#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

require(seqinr)
require(stringr)

# Verify that the correct number of arguments are given.
if(length(args)!=3){
  stop("The following arguments must be supplied: original reference sequence (complete),
       new reference sequence (possible missing segments), output path.", 
       call.=FALSE)
}

RefOldFile <- args[1]
RefNewFile <- args[2]
RefOutFile <- args[3]

# Read in original reference file and the newly generated reference.
RefOld <- read.fasta(RefOldFile, as.string = TRUE, forceDNAtolower = FALSE)
RefNew <- read.fasta(RefNewFile, as.string = TRUE, forceDNAtolower = FALSE)

# Compare each gene segment for the old and new reference sequences.
# If a gene segment in the new reference consists solely of N characters,
# then replace it with the segment from the old reference.
CompareSegment <- function(x){
  out <- RefNew[x]
  if(str_count(RefNew[[x]], LETTERS)[14]==nchar(RefNew[[x]])){
    out <- RefOld[x]
  }
  return(out)
}
RefOut <- sapply(seq(length(RefOld)), CompareSegment)

write.fasta(RefOut, attr(RefOut,"name"), RefOutFile, as.string=TRUE, nbchar=70)
