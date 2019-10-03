#!/bin/bash
#
VIRUS="$1";
ORGAN="$2";
#
TOTAL_COVERAGE=`awk '{sum += (($3-$2)*$4)} END {print sum}' ../output_data/TRACES_viral_bed/$VIRUS-coverage-$ORGAN.bed `;
TOTAL_SIZE=`tail -n 1 ../output_data/TRACES_viral_bed/$VIRUS-coverage-$ORGAN.bed | awk '{ print $3}'`;
NORMALIZED_COVERAGE=`echo "scale=4; $TOTAL_COVERAGE / $TOTAL_SIZE" | bc -l`;
#
printf "$VIRUS\t$ORGAN\t$NORMALIZED_COVERAGE\n" > ../output_data/TRACES_viral_statistics/$VIRUS-normalized-total-coverage-$ORGAN.txt;
#
