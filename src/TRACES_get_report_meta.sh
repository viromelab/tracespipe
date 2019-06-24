#!/bin/bash
cat ../meta_data/meta_info.txt | tr ':' '\t' | awk '{ print $1}' > names_tmp_TRACES.tmp
mapfile -t NAMES < names_tmp_TRACES.tmp
#
rm -f names2_tmp_TRACES.tmp
for name in "${NAMES[@]}" #
  do
  printf "../output_data/TRACES_results/REPORT_META_VIRAL_$name.txt " >> names2_tmp_TRACES.tmp
  done
paste -d "\t" `cat names2_tmp_TRACES.tmp` > ../output_data/TRACES_results/REPORT_META_VIRAL_ALL_SAMPLES.txt
rm -f names_tmp_TRACES.tmp names2_tmp_TRACES.tmp
