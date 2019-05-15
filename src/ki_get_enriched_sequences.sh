#!/bin/bash
rm -f en_seqs.fa;
while read line; do
efetch -db nucleotide -format fasta -id $line >> en_seqs.fa
done < ids_enrichment.txt 
cat $1 en_seqs.fa > TMP_VDB.fa
mv TMP_VDB.fa $1
