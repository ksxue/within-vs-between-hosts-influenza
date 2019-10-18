#!/bin/bash

# This script takes trimmed sequencing reads,
# aligns them to the specified reference,
# removes sequencing duplicates,
# and generates a pileup and gVCF file.
# Where possible or applicable, it follows the best practices described by GATK:
# https://software.broadinstitute.org/gatk/best-practices/bp_3step.php?case=GermShortWGS
# Note that this pipeline requires that the reference has already been indexed by
# both bowtie, samtools, and GATK.
# This pipeline for analysis of longitudinal samples remaps reads to the consensus reference
# of the initial timepoint sample, which must be specified in the input parameters.

# Load dependencies.
module load modules modules-init modules-gs
module load bowtie2/2.2.3
module load samtools/1.3
module load java/8u25
module load GATK/3.7
module load picard/1.43

# Input parameters.
fastq1="$1" # Give trimmed read 1 file.
fastq2="$2" # Give trimmed read 2 file.
sample="$3" # Give the sample name.
dir="$4" # Give the desired location of intermediate output files.
reference="$5" # Give reference path WITHOUT .fasta or .fa extension.
refname="$6" # Give a short nickname for the reference sequence.
bedfile="$7" # Give the BED format annotation file for the reference sequence.
initialtimepoint="$8" # Give the name of the initial timepoint sample.
clean="$9" # If 0, reruns entire pipeline and regenerates all intermediates.

# Other paths.
SummarizeBAM="bin/SummarizeBAM-1.21"
AnnotateVariants="bin/AnnotateVariants-1.13"
FillGapsRefSeq="scripts/FillGapsRefSeq.R"
hgreference="/net/shendure/vol10/nobackup/shared/alignments/bowtie2-2.0.2/human_g1k_hs37d5/hs37d5"

# Concatenate the sample and reference names
# so that all files reflect both the strain and reference.
sampleshort="$sample"
sample="$sample-$refname"

# Notifies user when all intermediate files are being regenerated.
if [ ${clean} -eq "0" ]; then
  echo "Pipeline will run in its entirety, overwriting all existing intermediates."
fi

:<<END
# Some analyses require filtering of reads that map to the human genome
# in order to remove all human genetic information.
# Since these reads were downloaded from the SRA, skip the filtering step.

# Filters out reads that map to the human genome.
# Retain unmapped reads for further analysis.
# Use a bowtie2 index built for human_g1k_hs37d5, a version of the hg19 reference genome.
echo "Filter out reads that map to the human genome."
if [ ! -f ${dir}/${sampleshort}-filtered.1.fastq.gz ] || [ ! -f ${dir}/${sampleshort}-filtered.2.fastq.gz ] \
  || [ ${clean} -eq "0" ];
then
  bowtie2 --very-fast-local --un-conc-gz ${dir}/${sampleshort}-filtered \
  -X 1000 \
  -x ${hgreference} \
  -1 ${fastq1} \
  -2 ${fastq2} \
  -S ${dir}/${projectdir}/${sampleshort}-hg.sam \
  2> ${dir}/${projectdir}/${sampleshort}-hg.bt2.log
  
  # Rename output files to have a .fastq.gz extension
  mv ${dir}/${sampleshort}-filtered.1 ${dir}/${sampleshort}-filtered.1.fastq.gz
  mv ${dir}/${sampleshort}-filtered.2 ${dir}/${sampleshort}-filtered.2.fastq.gz
fi
END

# Align reads with bowtie2.
# Convert SAM files to BAM files, sort, and delete SAM intermediates.
# Note that the GATK best practices suggest using bwa,
# but I have been having some trouble with picard running on bwa output.
echo "Align reads with bowtie2."
if [ ! -f ${dir}/${sample}.bt2.log ] || [ ${clean} -eq "0" ];
then
  bowtie2 --very-sensitive-local --un-conc-gz ${dir}/${sample}-unmapped \
	-X 2300 \
	-x ${reference} \
	-1 ${fastq1} \
	-2 ${fastq2} \
	-S ${dir}/${sample}.sam \
	2> ${dir}/${sample}.bt2.log
  samtools view -bS ${dir}/${sample}.sam -o ${dir}/${sample}.bam
  # Remove intermediate files.
  rm -f ${dir}/${sample}.sam
  rm -f ${dir}/${sample}-unmapped.1
  rm -f ${dir}/${sample}-unmapped.2
fi

# Tally reads at each position in the reference sequence.
# Use the custom script SummarizeBAM.
# Output a FASTA sequence corresponding to the consensus sequence for the population.
# Note that this new FASTA sequence will have the same coordinates
# as the original reference sequence.
echo "Summarize aligned reads using custom scripts."
if [ ! -f ${dir}/${sample}.fasta ] || [ ${clean} -eq "0" ];
then
  ${SummarizeBAM} -i <(samtools view ${dir}/${sample}.bam) \
    -f ${reference}.fasta -o ${dir}/${sample}.summary -s ${dir}/${sample}.fasta
  rm -f ${dir}/${sample}.bam
fi

# Some datasets contain sequence information for only a few gene segments.
# Segments that consist solely of N's cause problems for bowtie.
# Therefore, replace any unsequenced segments with the reference segment.
# This maintains the same genome coordinates.
echo "Replace missing segments in the new reference sequence."
if [ ! -f ${dir}/${sample}-clean.fasta ] || [ ${clean} -eq "0" ];
then
  Rscript ${FillGapsRefSeq} ${reference}.fasta ${dir}/${sample}.fasta \
    ${dir}/${sample}-clean.fasta
  # Build a reference from a FASTA sequence.
  bowtie2-build ${dir}/${sample}-clean.fasta ${dir}/${sample}-clean
fi

# Align all of the original, trimmed reads to the FASTA consensus sequence
# for the viral population at the initial timepoint.
echo "Re-align reads to FASTA consensus sequence for sample."
if [ ! -f ${dir}/${sample}-realigned.summary ] || [ ${clean} -eq "0" ];
then
  # Align reads.
  bowtie2 --very-sensitive-local --un-conc-gz ${dir}/${sample}-realigned-unmapped \
	-X 2300 \
	-x ${dir}/${initialtimepoint}-clean \
	-1 ${fastq1} \
	-2 ${fastq2} \
	-S ${dir}/${sample}-realigned.sam \
	2> ${dir}/${sample}-realigned.bt2.log
  samtools view -bS ${dir}/${sample}-realigned.sam -o ${dir}/${sample}-realigned.bam
  # Remove intermediate files.
  rm -f ${dir}/${sample}-realigned.sam
  rm -f ${dir}/${sample}-realigned-unmapped.1
  rm -f ${dir}/${sample}-realigned-unmapped.2
fi

# Remove duplicate reads, then sort and index BAM file.
echo "Remove duplicate reads."
if [ ! -f ${dir}/${sample}-realigned.summary ] || [ ${clean} -eq "0" ];
then
  # Sort and index BAM file.
  samtools sort ${dir}/${sample}-realigned.bam -o ${dir}/${sample}-realigned-sorted.bam
  samtools index ${dir}/${sample}-realigned-sorted.bam
  rm -f ${dir}/${sample}-realigned.bam
  # Remove sequencing duplicates.
  java -Xmx2g -jar ${PICARD_DIR}/MarkDuplicates.jar \
    INPUT=${dir}/${sample}-realigned-sorted.bam \
	OUTPUT=${dir}/${sample}-realigned-rmdup.bam \
	METRICS_FILE=${dir}/${sample}.picard \
	REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=TRUE
  rm -f ${dir}/${sample}-realigned-sorted.bam
  # Sort BAM file.
  samtools sort ${dir}/${sample}-realigned-rmdup.bam -o ${dir}/${sample}-clean.bam
  rm -f ${dir}/${sample}-realigned-rmdup.bam
  # Index BAM file.
  samtools index ${dir}/${sample}-clean.bam
fi

# Tally reads at each position in the population consensus sequence.
# Use the custom script SummarizeBAM.
echo "Summarize re-aligned reads using custom scripts."
if [ ! -f ${dir}/${sample}-realigned.summary ] || [ ${clean} -eq "0" ];
then
  ${SummarizeBAM} -i <(samtools view ${dir}/${sample}-clean.bam) \
    -f ${dir}/${initialtimepoint}-clean.fasta -o ${dir}/${sample}-realigned.summary
fi

# Annotate each variant as synonymous, nonsynonymous, etc.
echo "Annotate variants."
if [ ! -f ${dir}/${sample}-realigned-annotated.summary ] || [ ${clean} -eq "0" ];
then
  ${AnnotateVariants} -i ${dir}/${sample}-realigned.summary \
    -f ${dir}/${initialtimepoint}-clean.fasta -b ${bedfile} -o ${dir}/${sample}-realigned-annotated.summary
	
  # Annotate the annotation files with the sample name.
  sed -i "s/$/\t${sampleshort}/" \
	${dir}/${sample}-realigned-annotated.summary
fi

