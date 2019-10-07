#!/bin/bash
ORGAN=$1;
#
RESULT=`cat top-$ORGAN.csv \
| grep -a -e "imian virus 40" -e "imian_virus_40" -e "imian_Virus_40" -e "imian Virus 40" -e "simian virus 40" \
| awk '{ if($3 > 0 && $2 > 2000 && $2 < 6000) print $3"\t"$4; }' \
| head -n 1 \
| awk '{ print $1"\t"$2;}' \
| tr '_' '\t' \
| awk '{ print $1"\t"$2;}'`;
if [ -z "$RESULT" ]
  then
  echo -e "-\t-";
  else
  echo "$RESULT";
  fi
