#!/bin/bash
#
ORGAN="$1";
#
TOTAL_COVERAGE=`awk '{sum += (($3-$2)*$4)} END {print sum}' ../output_data/TRACES_mtdna_bed/mt-coverage-$ORGAN.bed `;
TOTAL_SIZE=`tail -n 1 ../output_data/TRACES_mtdna_bed/mt-coverage-$ORGAN.bed | awk '{ print $3}'`;
NORMALIZED_COVERAGE=`echo "scale=4; $TOTAL_COVERAGE / $TOTAL_SIZE" | bc -l`;
printf "mt\t$ORGAN\t$NORMALIZED_COVERAGE\n" > ../output_data/TRACES_mtdna_statistics/mt-total-depth-coverage-$ORGAN.txt;
#
ZERO_COVERAGE=`awk '{sum += ($3-$2)} END {print sum}' ../output_data/TRACES_mtdna_bed/mt-zero-coverage-$ORGAN.bed `;
if [ -z $ZERO_COVERAGE ] 
  then
 ZERO_COVERAGE=0;
  fi
NORMALIZED_ZERO_COVERAGE=`echo "scale=4; (100-(($ZERO_COVERAGE / $TOTAL_SIZE)*100))" | bc -l`;
printf "mt\t$ORGAN\t$NORMALIZED_ZERO_COVERAGE\n" > ../output_data/TRACES_mtdna_statistics/mt-total-horizontal-coverage-$ORGAN.txt;
#
