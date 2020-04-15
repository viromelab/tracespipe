#!/bin/bash
#
echo "Downloading Reference Human mithocondrial genome (mtDNA.fa) from NCBI ..."
#
rm -f mtDNA.fa;
efetch -db nucleotide -format fasta -id "NC_012920.1" > mtDNA.fa
#
echo "Done!";
#
