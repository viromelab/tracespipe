#!/bin/bash
mapfile -t ENTRIES < cy-calls.$1.bed  
#
for line in "${ENTRIES[@]}" #
  do
  #
  POSITIONS=`echo $line | awk '{ print $3 }'`;
  # echo "POS:$POSITIONS";
  cat ../system_files/cy_halotype_data.csv | grep "$POSITIONS";
  done
