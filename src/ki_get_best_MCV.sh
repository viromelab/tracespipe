#!/bin/bash
ORGAN=$1;
#
cat top-$ORGAN.csv \
| grep -e "MC_polyomavirus" -e "MC_virus" -e "Merkel_cell_polyomavirus" \
| awk '{ if($3 > 1 && $2 > 4000 && $2 < 7000) print $3"\t"$4; }' \
| head -n 1 \
| awk '{ print $1"\t"$2;}' \
| tr '_' '\t' \
| awk '{ print $1"\t"$2;}'