#!/bin/bash
ORGAN=$1;
#
cat top-$ORGAN.csv \
| grep -e "B19" -e "parvovirus" -e "Erythrovirus"  \
| awk '{ if($3 > 1 && $2 > 4000 && $2 < 6000) print $3"\t"$4; }' \
| head -n 1 \
| awk '{ print $1"\t"$2;}' \
| tr '_' '\t' \
| awk '{ print $1"\t"$2;}'
