require(tidyverse)
require(Biostrings)
require(seqinr)

# Write a function that, given a reference and test sequence in string form,
# returns the number of synonymous and nonsynonymous differences
# between the reference and test sequences.
CountSeqDifferences <- function(refseq, testseq){
  # Verify that the two sequences are the same length.
  if(nchar(refseq) != nchar(testseq)){
    return(c(NA, NA))
  } else{
    refseq <- s2c(refseq)
    testseq <- s2c(testseq)
    # Split the reference and test sequences into lists of codons.
    refcodons <- sapply(seq(1,length(refseq)-1,3), function(x) c2s(refseq[x:(x+2)]))
    testcodons <- sapply(seq(1,length(testseq)-1,3), function(x) c2s(testseq[x:(x+2)]))
    # Count the number of S and NS differences at each codon.
    # This is returned as a matrix with two rows, S and NS differences respectively.
    Diffs <- (base::sapply(seq(1:length(refcodons)), 
                           function(x) CountCodonDifferences(refcodons[x], testcodons[x])))
    return(c(sum(Diffs[1,]), sum(Diffs[2,]))) 
  }
}

# Given two codons, this function returns
# the number of synonymous and nonsynonymous differences
# between them, respectively.
CountCodonDifferences <- function(refcodon, testcodon){
  if(refcodon==testcodon){
    return(c(0,0))
  } else if (translate(s2c(refcodon))==translate(s2c(testcodon))) {
    return(c(1,0))
  } else {
    return(c(0,1))
  }
}