#!/bin/bash
ORGAN=$1;
#
RESULT=`cat top-$ORGAN.csv \
| grep -a -e "olyomavirus_3" -e "KI_polyomavirus" -e "KI_virus" -e "KIPyV" -e "NC_009238" \
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
