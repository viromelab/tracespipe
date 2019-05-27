#!/bin/bash
ORGAN=$1;
#
cat top-$ORGAN.csv \
| grep -e "nello" -e "orque" \
| awk '{ if($3 > 1 && $2 > 1500 && $2 < 5000) print $3"\t"$4; }' \
| head -n 1 \
| awk '{ print $2;}' \
| tr '_' '\t' \
| awk '{ print $1;}'
