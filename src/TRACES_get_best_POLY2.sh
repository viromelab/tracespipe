#!/bin/bash
ORGAN=$1;
#
RESULT=`cat top-$ORGAN.csv \
| grep -a -e "JC polyomavirus" -e "JC_polyomavirus" -e "JC_virus" -e "JC virus" -e "olyomavirus_2 " -e "olyomavirus 2 " -e "olyomavirus_2|" -e "olyomavirus 2|" -e "JCPyV" -e "NC_001699" -e "MF662195.1" -e "MF662183.1" -e "MF662192.1" -e "MF662193.1" \
| grep -a -e "omplete genome" -e "omplete_genome" \
| awk '{ if($3 > 0 && $2 > 3500 && $2 < 7000) print $3"\t"$4; }' \
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

