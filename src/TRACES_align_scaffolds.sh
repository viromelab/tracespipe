#!/bin/bash
#
REF="POLY.fa";
SCAFFOLDS="scaffolds.fasta";
THREADS="8";
#
bwa index $REF
bwa mem -t $THREADS $REF $SCAFFOLDS > scaffolds_aligned_$REF.sam
samtools sort scaffolds_aligned_$REF.sam > scaffolds_aligned_sorted_$REF.bam
samtools index scaffolds_aligned_sorted_$REF.bam scaffolds_aligned_sorted_$REF.bam.bai
#
