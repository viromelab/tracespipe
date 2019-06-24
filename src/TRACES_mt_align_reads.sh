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
#
rm -f index_file*
bowtie2-build $1 index_file
#
bowtie2 -a --threads $3 -x index_file -1 o_fw_pr.fq -2 o_rv_pr.fq -U o_fw_unpr.fq,o_rv_unpr.fq > aligned-$2.sam
samtools sort aligned-$2.sam > mt_aligned_sorted-$2.bam
rm -f aligned-$2.sam
#
samtools index mt_aligned_sorted-$2.bam mt_aligned_sorted-$2.bam.bai
#
