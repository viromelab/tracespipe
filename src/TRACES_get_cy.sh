#!/bin/bash
#
echo -e "\e[34m[TRACESPipe]\e[32m Download Reference Human Y chromosome (cy.fa) from NCBI ...\e[0m";
#
efetch -db nucleotide -format fasta -id "NC_000024.10" > cy.fa
#
#rm -f hs_ref_GRCh38.p12_chrY.fa.gz
#wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/H_sapiens/Assembled_chromosomes/seq/hs_ref_GRCh38.p12_chrY.fa.gz
#gunzip -f hs_ref_GRCh38.p12_chrY.fa.gz
#mv hs_ref_GRCh38.p12_chrY.fa cy.fa
#
echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
