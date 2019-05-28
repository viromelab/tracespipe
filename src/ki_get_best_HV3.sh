#!/bin/bash
ORGAN=$1;
#
cat top-$ORGAN.csv \
| grep -e "herpesvirus_3" -e "HHV-3"  \
| awk '{ if($3 > 1 && $2 > 100000 && $2 < 300000) print $3"\t"$4; }' \
| head -n 1 \
| awk '{ print $1"\t"$2;}' \
| tr '_' '\t' \
| awk '{ print $1"\t"$2;}'
