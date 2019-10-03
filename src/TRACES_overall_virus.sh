#!/bin/bash
#
VIRUS="$1";
#
TOTAL_COVERAGE=`awk '{sum += (($3-$2)*$4)} END {print sum}' ../output_data/TRACES_viral_bed/$VIRUS-coverage-testorgan.bed `;
TOTAL_SIZE=`tail -n 1 ../output_data/TRACES_viral_bed/$VIRUS-coverage-testorgan.bed | awk '{ print $3}'`;
NORMALIZED_COVERAGE=`echo "scale=4; $TOTAL_COVERAGE / $TOTAL_SIZE" | bc -l`;
#
printf "$VIRUS\t$NORMALIZED_COVERAGE\n" > ../output_data/TRACES_viral_bed/$VIRUS-normalized-total-coverage.txt;
#
