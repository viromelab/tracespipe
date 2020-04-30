#!/bin/bash
#
echo -e "\e[34m[TRACESPipe]\e[32m Downloading Reference Human mithocondrial genome (mtDNA.fa) from NCBI ...\e[0m";
#
rm -f mtDNA.fa;
efetch -db nucleotide -format fasta -id "NC_012920.1" > mtDNA.fa
#
echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
#
