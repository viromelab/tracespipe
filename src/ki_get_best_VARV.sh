#!/bin/bash
ORGAN=$1;
#
cat top-$ORGAN.csv \
| grep -a -e "VARV" -e "variola" -e "Variola"  \
| awk '{ if($3 > 0 && $2 > 100000 && $2 < 300000) print $3"\t"$4; }' \
| head -n 1 \
| awk '{ print $1"\t"$2;}' \
| tr '_' '\t' \
| awk '{ print $1"\t"$2;}'
