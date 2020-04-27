#!/bin/bash
# 
# SEE ALL DBS: update_blastdb.pl --showall
#
# INSTALL ======================================================================
#
mkdir -p ../system_files/blast_dbs
mkdir -p ../system_files/blast_dbs/nt
cd ../system_files/blast_dbs/nt/
update_blastdb.pl --decompress --passive --timeout 500 --force --verbose nt
cd ../../../src
#
# ==============================================================================
