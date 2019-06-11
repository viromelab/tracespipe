#!/bin/bash
#
echo "Downloading Phix reference genomes (F_PHIX.fa) from NCBI ...";
#
#gto_fasta_extract_read_by_pattern -p "phiX174" < VDB.fa > F_PHIX.fa
efetch -db nucleotide -format fasta -id "J02482.1,NC_001422.1" > F_PHIX.fa
#
echo "Done!";
