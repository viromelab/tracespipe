#!/bin/bash
#
echo "Changing Reference Human mithocondrial genome (mtDNA.fa) to $1 ..."
#
rm -f mtDNA.fa;
efetch -db nucleotide -format fasta -id "$1" > mtDNA.fa
#
echo "Done!";
#
