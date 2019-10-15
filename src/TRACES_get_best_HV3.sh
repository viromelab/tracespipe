#!/bin/bash
ORGAN=$1;
#
RESULT=`cat top-$ORGAN.csv \
| grep -a -e "erpesvirus_3" -e "erpesvirus 3" -e "HHV-3" -e "HHV3" -e "VZV" -e "Varicella" -e "Zoster" \
| awk '{ if($3 > 0 && $2 > 100000 && $2 < 220000) print $3"\t"$4; }' \
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
