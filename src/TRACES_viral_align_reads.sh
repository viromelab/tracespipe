#!/bin/bash
# 
# ALIGNMENT OF VIRAL READS
#
# IT ASSUMES THAT THE FOLLOWING FASTQ FILES EXIST (FOR INPUT):
# o_fw_pr.fq o_rv_pr.fq o_fw_unpr.fq o_rv_unpr.fq
#
# $1 -> REFERENCE FASTA FILE
# $2 -> ORGAN
# $3 -> VIRAL INITIALS
#
rm -f index-$2-$3-file*
#
# BUILD THE INDEX
bowtie2-build $1 index-$2-$3-file
#
# ALIGN
bowtie2 -a --threads 12 -x index-$2-$3-file -1 o_fw_pr.fq -2 o_rv_pr.fq -U o_fw_unpr.fq,o_rv_unpr.fq > aligned-$2-$3.sam
#
# SORT & BIN
samtools sort aligned-$2-$3.sam > viral_aligned_sorted-$2-$3.bam
#
# INDEX BAM
samtools index viral_aligned_sorted-$2-$3.bam viral_aligned_sorted-$2-$3.bam.bai
#
#
