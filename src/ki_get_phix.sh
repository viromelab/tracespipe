#!/bin/bash
#
echo "Extracting PhiX genomes from VDB.fa into F_PHIX.fa ...";
#
gto_fasta_extract_read_by_pattern -p "phiX174" < VDB.fa > F_PHIX.fa
#
echo "Done!";
