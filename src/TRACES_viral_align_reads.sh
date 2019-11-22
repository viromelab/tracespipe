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
# $4 -> NUMBER OF THREADS
#
rm -f index-$2-$3-file*
#
# BUILD THE INDEX
bowtie2-build $1 index-$2-$3-file
#
# ALIGN
bowtie2 -a --threads $4 -x index-$2-$3-file -1 o_fw_pr.fq -2 o_rv_pr.fq -U o_fw_unpr.fq,o_rv_unpr.fq > aligned-$2-$3.sam
#
# SORT & BIN
samtools sort --threads $4 aligned-$2-$3.sam > viral_aligned_sorted-$2-$3.bam
rm -f aligned-$2-$3.sam
#
if [[ "$5" -eq "1" ]];
  then
  # SORT BY NAME
  samtools sort --threads $4 -n -o viral_aligned_sorted-$2-$3.bam viral_aligned_sorted_sorted-$2-$3.bam
  #
  # ADD ms AND MC FOR MARKDUP
  samtools fixmate --threads $4 -m viral_aligned_sorted_sorted-$2-$3.bam viral_aligned_sorted_sorted-$2-$3-fixmate.bam
  rm -f viral_aligned_sorted_sorted-$2-$3.bam
  #
  # SORTING POSITION ORDER
  samtools sort --threads $4 -o viral_aligned_sorted_sorted-$2-$3-fixmate-sort.bam viral_aligned_sorted_sorted-$2-$3-fixmate.bam
  rm -f viral_aligned_sorted_sorted-$2-$3-fixmate.bam
  #
  # REMOVE DUPLICATES
  samtools markdup --threads $4 -r viral_aligned_sorted_sorted-$2-$3-fixmate-sort.bam viral_aligned_sorted-$2-$3.bam
  rm -f viral_aligned_sorted_sorted-$2-$3-fixmate-sort.bam
  #
  fi
#
# INDEX BAM
samtools index --threads $4 viral_aligned_sorted-$2-$3.bam viral_aligned_sorted-$2-$3.bam.bai
#
rm -f *.bt2
#
