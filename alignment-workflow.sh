#!/bin/bash
set -e
source /opt/sentieonStartup.sh
echo "Sentieon Script Built 2-18-21"
sentieon driver --version


echo "***************Aligning " $FASTQ1 $FASTQ2
(sentieon bwa mem -M -R $READGROUP \
  -t $NUMBER_THREADS /mnt/data/reference/$REFERENCE  /mnt/data/fastqs/R1.trimmed.fastq.gz /mnt/data/fastqs/R2.trimmed.fastq.gz || echo -n 'error during alignment') \
  | sentieon util sort -r /mnt/data/reference/$REFERENCE -o /mnt/data/output/$SAMPLE"_"$FLOWCELL.sorted.bam -t $NUMBER_THREADS --sam2bam -i -
