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
# $5 -> DUPLICATIONS
# $6 -> HIGH_SENSITIVITY
#
# BUILD THE INDEX
rm -f index-$2-$3-file* viral_aligned_sorted-$2-$3.bam.bai
bowtie2-build $1 index-$2-$3-file 
#
# ALIGN
if [[ "$6" == "1" ]];
  then
  bowtie2 -a --threads $4 --very-sensitive -x index-$2-$3-file -1 o_fw_pr.fq -2 o_rv_pr.fq -U o_fw_unpr.fq,o_rv_unpr.fq > aligned-$2-$3.sam
  else
  bowtie2 -a --threads $4 -x index-$2-$3-file -1 o_fw_pr.fq -2 o_rv_pr.fq -U o_fw_unpr.fq,o_rv_unpr.fq > aligned-$2-$3.sam
  fi
#
# SORT & BIN
samtools sort --threads $4 aligned-$2-$3.sam > viral_aligned_sorted-$2-$3.bam
rm -f aligned-$2-$3.sam
#
if [[ "$5" -eq "1" ]];
  then
  echo "Removing Duplications ...";
  # SORT BY NAME
  samtools sort -n viral_aligned_sorted-$2-$3.bam > viral_aligned_sorted_sorted-$2-$3.bam
  #
  # ADD ms AND MC FOR MARKDUP
  samtools fixmate -m viral_aligned_sorted_sorted-$2-$3.bam viral_aligned_sorted_sorted-$2-$3-fixmate.bam
  rm -f viral_aligned_sorted_sorted-$2-$3.bam
  #
  # SORTING POSITION ORDER
  samtools sort -o viral_aligned_sorted_sorted-$2-$3-fixmate-sort.bam viral_aligned_sorted_sorted-$2-$3-fixmate.bam
  rm -f viral_aligned_sorted_sorted-$2-$3-fixmate.bam
  #
  # REMOVE DUPLICATES
  samtools markdup -r viral_aligned_sorted_sorted-$2-$3-fixmate-sort.bam viral_aligned_sorted-$2-$3.bam
  rm -f viral_aligned_sorted_sorted-$2-$3-fixmate-sort.bam
  #
  fi
#
# INDEX BAM
samtools index -b viral_aligned_sorted-$2-$3.bam viral_aligned_sorted-$2-$3.bam.bai
#
mkdir -p ../output_data/TRACES_viral_statistics/
samtools flagstat viral_aligned_sorted-$2-$3.bam > ../output_data/TRACES_viral_statistics/Alignments-viral-$2-$3.txt
#
rm -f *.bt2
#
