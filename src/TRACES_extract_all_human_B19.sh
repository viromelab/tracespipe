#!/bin/bash
#
cat VDB.fa | grep -a -e "B19 " -e "parvovirus" -e "rythrovirus" -e "rythroparvovirus" | grep -a -v "olyomavirus" | grep -a -v "otavirus" | grep uman | tr ' ' '\t' | awk '{print $1;}' | sed 's/^.//' > access_codes_b19.txt
mapfile -t CODES < access_codes_b19.txt
#
rm -f B19_REFS.fa;
for xcode in "${CODES[@]}" # 
  do
  echo "Extracting $xcode on VDB.fa..."
  gto_fasta_extract_read_by_pattern -p "$xcode" < VDB.fa >> B19_REFS.fa
  echo " ";
  done
#
#
