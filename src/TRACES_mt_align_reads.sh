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
# $5 -> HIGH_SENSITIVITY = 1 ?
#
# INDEX
rm -f index_file* mt_aligned_sorted-$2.bam.bai
bowtie2-build $1 index_file
#
# ALIGN
if [[ "$5" == "1" ]];
  then
  bowtie2 -a --threads $3 --very-sensitive -x index_file -1 o_fw_pr.fq -2 o_rv_pr.fq -U o_fw_unpr.fq,o_rv_unpr.fq > aligned-$2.sam
  else
  bowtie2 -a --threads $3 -x index_file -1 o_fw_pr.fq -2 o_rv_pr.fq -U o_fw_unpr.fq,o_rv_unpr.fq > aligned-$2.sam
  fi
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
samtools flagstat -@ $3 mt_aligned_sorted-$2.bam > ../output_data/TRACES_mtdna_statistics/Alignments-mt-$2.txt
#
rm -f *.bt2
#
