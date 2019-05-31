#!/bin/bash
#
echo "Extracting HV3 genome from VDB.fa into HV3.fa ...";
#
gto_fasta_extract_read_by_pattern -p "X04370.1" < VDB.fa > HV3.fa
#
echo "Done!";
