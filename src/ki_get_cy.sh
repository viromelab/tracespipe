#!/bin/bash
#
rm -f hs_ref_GRCh38.p12_chrY.fa.gz
#
echo "Downloading Reference Human CY chromosome (cy.fa) from NCBI ..."
#
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/H_sapiens/Assembled_chromosomes/seq/hs_ref_GRCh38.p12_chrY.fa.gz
gunzip hs_ref_GRCh38.p12_chrY.fa.gz
mv hs_ref_GRCh38.p12_chrY.fa.gz cy.fa
#
echo "Done!";
