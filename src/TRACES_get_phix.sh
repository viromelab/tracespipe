#!/bin/bash
#
echo "Extracting PhiX genomes from VDB.fa into F_PHIX.fa ...";
#
#gto_fasta_extract_read_by_pattern -p "phiX174" < VDB.fa > F_PHIX.fa
efetch -db nucleotide -format fasta -id "J02482.1,NC_001422.1" > F_PHIX.fa
#
echo "Done!";
