#!/bin/bash
ORGAN=$1;
#
RESULT=`cat top-$ORGAN.csv \
| grep -a -e "olyomavirus_5" -e "olyomavirus 5" -e "MC_polyomavirus" -e "MC polyomavirus" -e "MC_virus" -e "MC virus" -e "Merkel_cell_polyomavirus" -e "Merkel cell polyomavirus" -e "MCPyV" -e "NC_010277" \
| grep -a -e "complete genome" -e "complete_genome" \
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

