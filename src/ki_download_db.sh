#!/bin/bash
#
# DOWNLOAD: https://www.ncbi.nlm.nih.gov/genome/viruses/
# CHECK ALSO THIS: https://www.viprbrc.org/brc/vipr_allSpecies_search.spg?method=SubmitForm&decorator=vipr
#
tail -n +3 taxid10239.nbr > DATA.nbr
#
rm -f DATABASE.fa
#
mapfile -t ENTRIES < DATA.nbr
#
for line in "${ENTRIES[@]}" #
  do
  #
  GID=`echo $line | awk '{ print $2 }'`;
  efetch -db nucleotide -format fasta -id "$GID" >> DATABASE.fa
  done
#
