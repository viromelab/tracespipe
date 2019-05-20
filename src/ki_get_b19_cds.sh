#!/bin/bash
#
echo "Extracting B19 coding sequence from VDB.fa into B19_cds.fa ...";
#
gto_fasta_extract_read_by_pattern -p "AY044266.2" < VDB.fa > B19_cds.fa
#
echo "Done!";
