#!/bin/bash
ORGAN=$1;
#
RESULT=`cat top-$ORGAN.csv \
| grep -a -e "Polyomavirus_12" -e "polyomavirus_12" -e "Polyoma_12" -e "polyoma_12" -e "HPyV12" -e "NC_020890" \
| awk '{ if($3 > 0 && $2 > 2500 && $2 < 7000) print $3"\t"$4; }' \
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
