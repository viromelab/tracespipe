#!/bin/bash
ORGAN=$1;
#
RESULT=`cat top-$ORGAN.csv \
| grep -a -e "_TT_virus" -e " TT virus" -e "anello" -e "Anello" -e " torque " -e " Torque " -e "_torque_" -e "_Torque_" \
| grep -a -e "human" -e "Human" -e "omo sapiens" -e "omo_sapiens" \
| grep -a -e "omplete genome" -e "omplete_genome" \
| awk '{ if($3 > 0 && $2 > 2000 && $2 < 5000) print $3"\t"$4; }' \
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

