#!/bin/bash
ORGAN=$1;
#
RESULT=`cat top-$ORGAN.csv \
| grep -a -e "uman bocavirus 2" -e "uman bocavirus 3" -e "uman bocavirus 4" -e "HBoV2" -e "HBoV3" -e "HBoV4" \
| awk '{ if($3 > 0 && $2 > 4000 && $2 < 6000) print $3"\t"$4; }' \
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
