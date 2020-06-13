#!/bin/bash
ORGAN=$1;
#
RESULT=`cat top-$ORGAN.csv \
| grep -a -e "bocavirus 2" -e "bocavirus 3" -e "bocavirus 4" -e "Bocavirus 2" -e "Bocavirus 3" -e "Bocavirus 4" -e " HBoV2" -e " HBoV3" -e " HBoV4" \
| grep -a -e "human" -e "Human" \
| grep -a -e "omplete genome" -e "omplete_genome" \
| awk '{ if($3 > 0 && $2 > 4000 && $2 < 6000) print $3"\t"$4; }' \
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

