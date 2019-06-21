#!/bin/bash
#
# DOWNLOAD: https://www.ncbi.nlm.nih.gov/genome/viruses/
#
# https://www.biostars.org/p/259792/
#
# CHECK ALSO THIS: https://www.viprbrc.org/brc/vipr_allSpecies_search.spg?method=SubmitForm&decorator=vipr
#
tail -n +3 ../system_files/taxid10239.nbr > DATA.nbr
#
rm -f VDB.fa
#
mapfile -t ENTRIES < DATA.nbr
#
for line in "${ENTRIES[@]}" #
  do
  #
  GID=`echo $line | awk '{ print $2 }'`;
  efetch -db nucleotide -format fasta -id "$GID" >> VDB.fa
  done
#
