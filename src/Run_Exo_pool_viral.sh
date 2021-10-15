#!/bin/bash
#
ls -la ../output_data/TRACES_results/top-* | awk ' { print $9;}' > all_tops.txt
#
while read p; do
  echo "$p";
  ./TRACESQuick_extra_vir_view.sh "$p"
done <all_tops.txt
#
