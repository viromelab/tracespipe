#!/bin/bash
# 
# ALIGNMENT OF Y-CHROMOSOME READS
#
# IT ASSUMES THAT THE FOLLOWING FASTQ FILES EXIST (FOR INPUT):
# o_fw_pr.fq o_rv_pr.fq o_fw_unpr.fq o_rv_unpr.fq
#
# $1 -> REFERENCE FASTA FILE
# $2 -> ORGAN
# $3 -> NUMBER OF THREADS
# $4 -> HIGH_SENSITIVITY = 1 ?
#
rm -f cy_index_file* cy_aligned_sorted-$2.bam.bai
bowtie2-build $1 cy_index_file
#
if [[ "$4" == "1" ]];
  then
  bowtie2 -a --very-sensitive --threads $3 -x cy_index_file -1 o_fw_pr.fq -2 o_rv_pr.fq -U o_fw_unpr.fq,o_rv_unpr.fq > cy_aligned-$2.sam
  else
  bowtie2 -a --threads $3 -x cy_index_file -1 o_fw_pr.fq -2 o_rv_pr.fq -U o_fw_unpr.fq,o_rv_unpr.fq > cy_aligned-$2.sam
  fi
samtools sort cy_aligned-$2.sam > cy_aligned_sorted-$2.bam
#
samtools index cy_aligned_sorted-$2.bam cy_aligned_sorted-$2.bam.bai
#
