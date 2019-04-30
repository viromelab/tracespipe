#!/bin/bash
# 
# IT ASSUMES THAT THE FOLLOWING OUTPUT FILES EXIST:
# NP-o_fw_pr.fq NP-o_fw_unpr.fq NP-o_rv_pr.fq NP-o_rv_unpr.fq
#

ORGAN="$1";
RNAMES="NP-o_fw_pr.fq:NP-o_fw_unpr.fq:NP-o_rv_pr.fq:NP-o_rv_unpr.fq";

## RUN METAGENOMIC COMPOSITION
echo "RUNNING FALCON ...";
FALCON -v -n 12 -t 38 -F -Z -l 47 -c 20 -x top-$ORGAN.csv -y $ORGAN.com $RNAMES VDB.fa
FALCON-filter -v -F -t 1.0 -o $ORGAN.pos $ORGAN.com
FALCON-filter-visual -v -e 1 -F -o $ORGAN.svg $ORGAN.pos

## CONVERT SVG OUTPUT TO PDF
rsvg-convert -f pdf -o $ORGAN.pdf $ORGAN.svg

## RUN GULL FOR INTER-GENOMIC SIMILARITY ANALYSIS
cat top-$ORGAN.csv | awk '{ if($3 > 1) print $1"\t"$2"\t"$3"\t"$4; }' \
| awk '{ print $4;}' | tr '|' '\t' | tr '_' '\t' | awk '{ print $1;}' > GIS-$ORGAN;
idx=0;
cat GIS-$ORGAN | while read line
  do
  namex=`echo $line | tr ' ' '_'`;
  if [[ "$idx" -eq "0" ]];
    then
    printf "%s" "$namex" > $ORGAN-FNAMES.fil;
    else
    printf ":%s" "$namex" >> $ORGAN-FNAMES.fil;
    fi
  gto_fasta_extract_read_by_pattern -p "$line" < VDB.fa > $namex;
  ((idx++));
  done
FALCON-inter -v -m 6:1:1:0/0 -m 13:20:1:3/10 -m 20:100:1:5/10 -c 30 -n 8 -x $ORGAN-MATRIX.csv `cat $ORGAN-FNAMES.fil`
FALCON-inter-visual -v -w 25 -a 8 -x $ORGAN-HEAT.svg $ORGAN-MATRIX.csv

## CONVERT SVG OUTPUT TO PDF
rsvg-convert -f pdf -o $ORGAN-HEAT.pdf $ORGAN-HEAT.svg

