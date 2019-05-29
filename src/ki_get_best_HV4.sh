#!/bin/bash
ORGAN=$1;
#
cat top-$ORGAN.csv \
| grep -a -e "Epstein" -e "epstein" -e "herpesvirus_4" -e "HHV-4" -e "HHV4" -e "EBV" \
| awk '{ if($3 > 0 && $2 > 100000 && $2 < 220000) print $3"\t"$4; }' \
| head -n 1 \
| awk '{ print $1"\t"$2;}' \
| tr '_' '\t' \
| awk '{ print $1"\t"$2;}'
