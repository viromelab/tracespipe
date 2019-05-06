#!/bin/bash
#
# IT ASSUMES THAT THE FOLLOWING FASTQ FILES EXIST (FOR INPUT):
# o_fw_pr.fq o_rv_pr.fq o_fw_unpr.fq o_rv_unpr.fq
#
# $1 -> REFERENCE FASTA FILE
# $2 -> ORGAN
#
bowtie2-build $1 index_file
#
bowtie2 -a --threads 12 -x index_file -1 o_fw_pr.fq -2 o_rv_pr.fq -U o_fw_unpr.fq,o_rv_unpr.fq > aligned-$ORGAN.sam
samtools sort aligned-$ORGAN.sam > aligned_sorted-$ORGAN.bam
#
samtools index aligned_sorted-$ORGAN.bam aligned_sorted-$ORGAN.bam.bai
#
