#!/bin/bash
#  
zcat VDB.fa.gz | ./gto_fasta_extract_read_by_pattern -p "phiX174" > F_PHIX.fa
#
