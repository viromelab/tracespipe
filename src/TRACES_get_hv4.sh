#!/bin/bash
#
echo "Extracting HV4 genome from VDB.fa into HV4.fa ...";
#
gto_fasta_extract_read_by_pattern -p "DQ279927.1" < VDB.fa > HV4.fa
#
echo "Done!";
