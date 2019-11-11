#!/bin/bash
#
# $1 -> REF ;	
# $2 -> SCAFFOLDS ;	
# $3 -> THREADS ;
# $4 -> ORGAN ;
#
REF="$1";
SCAFFOLDS="$2";
THREADS="$3";
ORGAN="$4";
#
bwa index $REF.fa
bwa mem -t $THREADS $REF.fa $SCAFFOLDS > scaffolds_aligned_$REF-$ORGAN.sam
samtools sort scaffolds_aligned_$REF-$ORGAN.sam > scaffolds_aligned_sorted_$REF-$ORGAN.bam
samtools index scaffolds_aligned_sorted_$REF-$ORGAN.bam scaffolds_aligned_sorted_$REF-$ORGAN.bam.bai
mv scaffolds_aligned_sorted_$REF-$ORGAN.bam ../output_data/TRACES_hybrid_$ORGAN/
mv scaffolds_aligned_sorted_$REF-$ORGAN.bam.bai ../output_data/TRACES_hybrid_$ORGAN/
cp $REF.fa ../output_data/TRACES_hybrid_$ORGAN/
cp $REF.fa.fai ../output_data/TRACES_hybrid_$ORGAN/
#
rm -f scaffolds_aligned_$REF-$ORGAN.sam;
#
