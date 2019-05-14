#!/bin/bash
rm -f en_seqs.fa;
while read line; do
efetch -db nucleotide -format fasta -id $line >> en_seqs.fa
done < ids_enrichment.txt 
cat VDB.fa en_seqs.fa > TMP_VDB.fa
mv TMP_VDB.fa VDB.fa
