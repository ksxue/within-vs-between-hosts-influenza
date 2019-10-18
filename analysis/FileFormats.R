# List some standardized tab-delimited file formats
# for within-host work.

# The output of the sequencing pipeline is an aligned, annotated summary file
# of base counts at each site. The columns are ordered as below:
RealignedAnnotatedSummaryCols <- 
  c("Chr","Pos","Base","RefBase","GenomePos","Count","AvgQ","AvgReadPos",
    "Gene","Codon","RefAA","AltAA","Syn","FourfoldSyn","RefCodon","AltCodon","CodonPos",
    "Sample")

# The hard-filter variant-calling script produces variant call files.
# The columns are ordered as below:
VariantCols <-
  c("Sample","Chr","GenomePos","Gene","Pos","Consensus","Base","Codon","CodonPos","RefAA","AltAA",
    "RefCodon","AltCodon","Count","AvgQ","AvgReadPos","Coverage","Freq")

# The coverage summary script produces a table with columns ordered as below:
CovSummaryCols <-
  c("Sample","Chr","MeanCov","MedCov","IQRCov","MinCov","MaxCov","PercentSitesAboveMin")

# Metadata files for each project contain file paths and information
# about each sample, i.e. subtype.
# The columns are ordered as below:
MetadataCols <- c("Fastq1","Fastq2","Trimmed1","Trimmed2","Sample","Project",
              "Ref","Subtype")

# Calculate the length of the coding sequence for each gene.
# Use the BED file for each gene to do this calculation.
BEDCols <- c("Chr","ChrStart","ChrEnd","Gene","Score","Strand","ThickStart","ThickEnd",
             "ItemRGB","BlockCount","BlockSizes","BlockStarts")

# The information contained in the FASTA headers downloaded from GISAID
# are listed in the fields below.
# The fields are delimited by the | character.
GISAIDFASTAHeader <- c("IsolateName","DNAAccession","Type","Lineage","CollectionDate",
                 "Segment","PassageHistory","OriginatingLab","SubmittingLab")

# List influenza full-length genes.
FluGenes <- c("1-PB2","2-PB1","3-PA","4-HA","5-NP",
              "6-NA","7-M1","7-M2","8-NS1","8-NEP")
# List non-overlapping influenza genes.
FluGenesNonOverlapping <- c("1-PB2","2-PB1","3-PA","4-HA","5-NP","6-NA")


