#!/bin/bash
ORGAN=$1;
#
RESULT=`cat top-$ORGAN.csv \
| grep -a -e "BK_polyoma" -e "BK_virus" -e "BKPyV" -e "NC_001538" \
| grep "uman" \
| awk '{ if($3 > 0 && $2 > 3000 && $2 < 7000) print $3"\t"$4; }' \
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
