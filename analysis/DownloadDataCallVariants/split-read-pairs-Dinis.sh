#!/bin/bash

# Given a single .fastq.gz file with read 1 and read 2 combined
# but not interleaved,
# combine the reads onto a single line,
# sort the reads,
# split read pairs into separate files,
# and break reads back into FASTQ format.

fastq="$1"
outstem="$2"

# Use paste to combine each FASTQ read, which spans 4 lines,
# into a single line with comma-delimited parts.
zless ${fastq} | tr ' ' '|' | paste - - - - -d'\t' | sort > ${outstem}-sorted.fastq

# Iterate through the sorted reads.
# Identify paired reads and output them onto the same line.
prevRead=""
prevReadID=""
while read currRead
do
  currReadID=${currRead%%|*}
  if [ "${currReadID}" == "${prevReadID}" ]; then
	echo -e ${prevRead}"\t\t"${currRead}
  fi
  prevRead=${currRead}
  prevReadID=${currReadID}
done < ${outstem}-sorted.fastq > ${outstem}-paired.fastq

# Split read pairs 1 and 2 into separate files.
# Restore them into standard FASTQ format, with each read
# split across four lines.
awk -F"\t\t" '{print $1}' ${outstem}-paired.fastq |
  tr ' ' '\n' | tr '|' ' ' > ${outstem}_1.fastq
gzip ${outstem}_1.fastq
awk -F"\t\t" '{print $2}' ${outstem}-paired.fastq |
  tr ' ' '\n' | tr '|' ' ' > ${outstem}_2.fastq
gzip ${outstem}_2.fastq

# Remove intermediate files.
rm -f ${outstem}-sorted.fastq
rm -f ${outstem}-paired.fastq
