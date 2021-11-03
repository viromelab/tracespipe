#!/bin/bash
#
if [ -z "$1" ]
  then
  echo "ERROR: No coverage BED file provided."
  echo "Example: ./TRACES_normalized_depth.sh ../output_data/TRACES_viral_bed/B19-coverage-blood.bed 500"
  exit;
fi
#
if [ -z "$2" ]
  then
  echo "ERROR: NO maximum depth provided."
  echo "Example: ./TRACES_normalized_depth.sh ../output_data/TRACES_viral_bed/B19-coverage-blood.bed 500"
  exit;
fi
#
if [ ! -f "$1" ]; then
  echo "ERROR: Input file $1 does not exists!";
  exit;
fi
#
./TRACES_project_coordinates.sh $1 $2 | awk '{sum+=$2} END { print "Depth: ",sum/NR}'
#
