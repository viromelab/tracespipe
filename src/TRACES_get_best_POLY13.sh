#!/bin/bash
ORGAN=$1;
#
RESULT=`cat top-$ORGAN.csv \
| grep -a -e "Polyomavirus_13" -e "polyomavirus_13" -e "Polyoma_13" -e "polyoma_13" -e "NJPyV" -e "NC_024118" \
| awk '{ if($3 > 0 && $2 > 2800 && $2 < 7000) print $3"\t"$4; }' \
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
