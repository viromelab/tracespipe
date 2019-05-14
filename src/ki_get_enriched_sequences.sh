#!/bin/bash
rm -f en_seqs.fa;
while read line; do
efetch -db nucleotide -format fasta -id $line >> en_seqs.fa
done < ids_enrichment.txt 
