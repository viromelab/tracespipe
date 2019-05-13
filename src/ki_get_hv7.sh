#!/bin/bash
#
echo "Extracting HV7 genome from VDB.fa into HV7.fa ...";
#
gto_fasta_extract_read_by_pattern -p "AF037218" < VDB.fa > HV7.fa
#
echo "Done!";
