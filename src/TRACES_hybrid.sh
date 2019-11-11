#!/bin/bash
#
REF="$1";
SCAFFOLDS="$2";
THREADS="$3";
ORGAN="$4";
#
#
if [ ! -f "$REF.fa" ];
  then
  echo -e "\e[38;5;208mWARNING: $REF reference file not found!\e[0m"
  echo "TIP: before this, run: ./TRACESPipe.sh --run-meta --run-all-v-alig"
  echo -e "\e[38;5;198mDon't panic! This warning may be cause by absence of aligned results\e[0m."
  echo "For addition information, see the instructions at the web page."
  exit 1;
  fi
#
if [ ! -f "$SCAFFOLDS" ];
  then
  echo -e "\e[38;5;208mWARNING: $SCAFFOLDS not found!\e[0m"
  echo "TIP: before this, make sure this runs: ./TRACESPipe.sh --run-de-novo"
  echo -e "\e[38;5;198mAnother reason for this warning may be de-novo assembly empty results\e[0m."
  fi
#
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
#
