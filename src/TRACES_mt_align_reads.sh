#!/bin/bash
# 
# ALIGNMENT OF MITOCHONDRIAL READS
#
# IT ASSUMES THAT THE FOLLOWING FASTQ FILES EXIST (FOR INPUT):
# o_fw_pr.fq o_rv_pr.fq o_fw_unpr.fq o_rv_unpr.fq
#
# $1 -> REFERENCE FASTA FILE
# $2 -> ORGAN
# $3 -> NUMBER OF THREADS
# $4 -> REMOVE DUPLICATIONS
#
# INDEX
rm -f index_file*
bowtie2-build $1 index_file
#
# ALIGN
bowtie2 -a --threads $3 --very-sensitive -x index_file -1 o_fw_pr.fq -2 o_rv_pr.fq -U o_fw_unpr.fq,o_rv_unpr.fq > aligned-$2.sam
#
# SORT & BIN
samtools sort --threads $3 aligned-$2.sam > mt_aligned_sorted-$2.bam
rm -f aligned-$2.sam
#
if [[ "$4" -eq "1" ]];
  then
  echo "Removing Duplications ...";
  # SORT BY NAME
  samtools sort --threads $3 -n mt_aligned_sorted-$2.bam > mt_aligned_sorted_sorted-$2.bam
  #
  # ADD ms AND MC FOR MARKDUP
  samtools fixmate --threads $3 -m mt_aligned_sorted_sorted-$2.bam mt_aligned_sorted_sorted-$2-fixmate.bam
  rm -f mt_aligned_sorted_sorted-$2.bam
  #
  # SORTING POSITION ORDER
  samtools sort --threads $3 -o mt_aligned_sorted_sorted-$2-fixmate-sort.bam mt_aligned_sorted_sorted-$2-fixmate.bam
  rm -f mt_aligned_sorted_sorted-$2-fixmate.bam
  #
  # REMOVE DUPLICATES
  samtools markdup --threads $3 -r mt_aligned_sorted_sorted-$2-fixmate-sort.bam mt_aligned_sorted-$2.bam
  rm -f mt_aligned_sorted_sorted-$2-fixmate-sort.bam
  #
  fi
# INDEX BAM
samtools index -@ $3 mt_aligned_sorted-$2.bam mt_aligned_sorted-$2.bam.bai
#
rm -f *.bt2
#
