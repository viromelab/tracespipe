#!/bin/bash
ORGAN=$1;
#
RESULT=`cat top-$ORGAN.csv \
| grep -a -e " B19 " -e "parvovirus" -e "Parvovirus" -e "erythrovirus" -e "Erythrovirus" \
| grep -a -e "complete genome" -e "complete_genome" -e "Complete genome" -e "Complete_genome" \
| grep -a -v "Polyomavirus" \ 
| grep -a -v "polyomavirus" \ 
| grep -a -v "bocavirus" \
| grep -a -v "Bocavirus" \
| grep -a -v "bocaparvovirus" \
| grep -a -v "Bocaparvovirus" \
| awk '{ if($3 > 0 && $2 > 5300 && $2 < 5900) print $3"\t"$4; }' \
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
