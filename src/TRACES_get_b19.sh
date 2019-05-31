#!/bin/bash
#
echo "Extracting B19 genome from VDB.fa into B19.fa ...";
#
gto_fasta_extract_read_by_pattern -p "AY386330.1" < VDB.fa > B19.fa
#
echo "Done!";
