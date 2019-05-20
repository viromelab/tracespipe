#!/bin/bash
#
echo "Extracting JCV genome from VDB.fa into JCV.fa ...";
#
gto_fasta_extract_read_by_pattern -p "AY044266.2" < VDB.fa > JCV.fa
#
echo "Done!";
