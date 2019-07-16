#!/bin/bash
#
# accessions.txt with GIDs (line by line)
#
mapfile -t CODES < accessions.txt
#
for xcode in "${CODES[@]}" # 
  do
  echo "Searching $xcode on VDB.fa..."
  grep $xcode VDB.fa
  echo " ";
  done
#
