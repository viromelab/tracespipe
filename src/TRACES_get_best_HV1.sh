#!/bin/bash
ORGAN=$1;
#
RESULT=`cat top-$ORGAN.csv \
| grep -a -e " HHV-1" -e " HHV1" -e " HSV1" -e " HSV-1" -e "erpesvirus_1" -e "erpesvirus 1" \
| grep -a -e "human" -e "Human" -e "omo sapiens" -e "omo_sapiens" \
| grep -a -e "omplete genome" -e "omplete_genome" \
| awk '{ if($3 > 0 && $2 > 100000 && $2 < 300000) print $3"\t"$4; }' \
| head -n 1 \
| awk '{ print $1"\t"$2;}' \
| sed "s/NC\_/NC-/" \
| tr '_' '\t' \
| awk '{ print $1"\t"$2;}'`;
if [ -z "$RESULT" ]
  then
  echo -e "-\t-";
  else
  echo "$RESULT" | sed "s/NC-/NC\_/"
  fi

