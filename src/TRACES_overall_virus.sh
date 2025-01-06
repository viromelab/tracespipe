#!/bin/bash
#
VIRUS="$1";
ORGAN="$2";
MAX_DEPTH=1000;
#
#TOTAL_COVERAGE=`awk '{sum += (($3-$2)*$4)} END {print sum}' ../output_data/TRACES_viral_bed/$VIRUS-coverage-$ORGAN.bed `;
TOTAL_SIZE=`tail -n 1 ../output_data/TRACES_viral_bed/$VIRUS-coverage-$ORGAN.bed | awk '{ print $3}'`;
acc=$(head -1 ../output_data/TRACES_viral_bed/$VIRUS-coverage-$ORGAN.bed | cut -f1);
#Updated Depth calculation which ignores peaks, but more importantly ignores secondary/supplementary alignments
TOTAL_COVERAGE=$(samtools bedcov --max-depth $MAX_DEPTH <(echo -e "$acc\t0\t$TOTAL_SIZE") ../output_data/TRACES_viral_alignments/viral_aligned_sorted-$ORGAN-$VIRUS.bam | awk '{print $NF}')
NORMALIZED_COVERAGE=`echo "scale=4; $TOTAL_COVERAGE / $TOTAL_SIZE" | bc -l`;
printf "$VIRUS\t$ORGAN\t$NORMALIZED_COVERAGE\n" > ../output_data/TRACES_viral_statistics/$VIRUS-total-depth-coverage-$ORGAN.txt;
#
ZERO_COVERAGE=`awk '{sum += ($3-$2)} END {print sum}' ../output_data/TRACES_viral_bed/$VIRUS-zero-coverage-$ORGAN.bed `;
if [ -z $ZERO_COVERAGE ]
  then
 ZERO_COVERAGE=0;
  fi
NORMALIZED_ZERO_COVERAGE=`echo "scale=4; (100-(($ZERO_COVERAGE / $TOTAL_SIZE)*100))" | bc -l`;
printf "$VIRUS\t$ORGAN\t$NORMALIZED_ZERO_COVERAGE\n" > ../output_data/TRACES_viral_statistics/$VIRUS-total-horizontal-coverage-$ORGAN.txt;
#
