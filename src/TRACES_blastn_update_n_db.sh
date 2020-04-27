#!/bin/bash
# 
# UPDATE =======================================================================
#
mkdir -p ../system_files/blast_dbs
mkdir -p ../system_files/blast_dbs/nt
cd ../system_files/blast_dbs/nt/
if [ ! -f nt..ndb ];
  then
  echo -e "\e[31mERROR: Unable to find nt blast database!\e[0m"
  echo "TIP: before this, run: ./TRACESPipe.sh --install-blast"
  echo "For addition information, see the instructions at the web page."
  cd ../../../src
  exit 1;
else
  update_blastdb.pl --decompress --passive --timeout 500 --force --verbose nt
  cd ../../../src
fi
#
# ==============================================================================
