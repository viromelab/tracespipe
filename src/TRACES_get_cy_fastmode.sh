#!/bin/bash
#
echo "Downloading Reference Human Y chromosome (cy.fa) from NCBI ..."
#
rm -f hs_ref_GRCh38.p12_chrY.fa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/H_sapiens/Assembled_chromosomes/seq/hs_ref_GRCh38.p12_chrY.fa.gz
gunzip -f hs_ref_GRCh38.p12_chrY.fa.gz
mv hs_ref_GRCh38.p12_chrY.fa cy.fa
#
echo "Done!";
