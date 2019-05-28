#!/bin/bash
ORGAN=$1;
#
cat top-$ORGAN.csv \
| grep -e "Epstein" -e "epstein" -e "herpesvirus_4"" \
| awk '{ if($3 > 1 && $2 > 110000 && $2 < 300000) print $3"\t"$4; }' \
| head -n 1 \
| awk '{ print $2;}' \
| tr '_' '\t' \
| awk '{ print $1;}'
