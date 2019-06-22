#!/bin/bash
# 
# ALIGNMENT OF Y-CHROMOSOME READS
#
# IT ASSUMES THAT THE FOLLOWING FASTQ FILES EXIST (FOR INPUT):
# o_fw_pr.fq o_rv_pr.fq o_fw_unpr.fq o_rv_unpr.fq
#
# $1 -> REFERENCE FASTA FILE
# $2 -> ORGAN
#
rm -f cy_index_file*
bowtie2-build $1 cy_index_file
#
bowtie2 -a --threads 12 -x cy_index_file -1 o_fw_pr.fq -2 o_rv_pr.fq -U o_fw_unpr.fq,o_rv_unpr.fq > cy_aligned-$2.sam
samtools sort cy_aligned-$2.sam > cy_aligned_sorted-$2.bam
#
samtools index cy_aligned_sorted-$2.bam cy_aligned_sorted-$2.bam.bai
#
