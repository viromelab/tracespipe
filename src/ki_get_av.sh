#!/bin/bash
#
echo "Extracting AV genome from VDB.fa into AV.fa ...";
#
gto_fasta_extract_read_by_pattern -p "MH648911.1" < VDB.fa > AV.fa
#
echo "Done!";
