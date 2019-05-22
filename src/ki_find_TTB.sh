#!/bin/bash
ORGAN=$1;
#
cat top-$ORGAN.csv \
| grep -e "nello" -e "orque" \
| awk '{ if($3 > 1 && $2 > 1500 && $2 < 5000) print $1"\t"$2"\t"$3"\t"$4; }'
#
