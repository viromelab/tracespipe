#!/bin/bash
#
VIRUS="$1";
ORGAN="$2";
#
TOTAL_COVERAGE=`awk '{sum += (($3-$2)*$4)} END {print sum}' ../output_data/TRACES_viral_bed/$VIRUS-coverage-$ORGAN.bed `;
TOTAL_SIZE=`tail -n 1 ../output_data/TRACES_viral_bed/$VIRUS-coverage-$ORGAN.bed | awk '{ print $3}'`;
NORMALIZED_COVERAGE=`echo "scale=4; $TOTAL_COVERAGE / $TOTAL_SIZE" | bc -l`;
printf "$VIRUS\t$ORGAN\t$NORMALIZED_COVERAGE\n" > ../output_data/TRACES_viral_statistics/$VIRUS-total-depth-coverage-$ORGAN.txt;
#
ZERO_COVERAGE=`awk '{sum += ($3-$2)} END {print sum}' ../output_data/TRACES_viral_bed/$VIRUS-zero-coverage-$ORGAN.bed `;
NORMALIZED_ZERO_COVERAGE=`echo "scale=4; $ZERO_COVERAGE / $TOTAL_SIZE" | bc -l`;
printf "$VIRUS\t$ORGAN\t$NORMALIZED_ZERO_COVERAGE\n" > ../output_data/TRACES_viral_statistics/$VIRUS-total-horizontal-coverage-$ORGAN.txt;
#
