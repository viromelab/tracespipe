#!/bin/bash
#
rm -f hs_ref_GRCh38.p12_chrMT.fa.gz
#
echo "Downloading Reference Human mithocondrial genome (mtDNA.fa) from NCBI ..."
#
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/H_sapiens/Assembled_chromosomes/seq/hs_ref_GRCh38.p12_chrMT.fa.gz
gunzip hs_ref_GRCh38.p12_chrMT.fa.gz
mv hs_ref_GRCh38.p12_chrMT.fa mtDNA.fa
#
echo "Done!";
