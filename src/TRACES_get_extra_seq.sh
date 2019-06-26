#!/bin/bash
#
gto_fasta_extract_read_by_pattern -p "$2" < $1 > tmp_exists_or_not.fa
FIRST_LINE=`head -n 1 tmp_exists_or_not.fa > tmp_first_line.fa`;
if [ `wc -l tmp_exists_or_not.fa` -ge "2" ]
then
  echo "ERROR: This sequence seems to already exist in the VDB.fa.";
  echo "Header: $FIRST_LINE";
  exit 1;
fi
rm -f en_seqs.fa;
efetch -db nucleotide -format fasta -id $2 > en_seqs.fa
cat $1 en_seqs.fa > TMP_VDB.fa
mv TMP_VDB.fa $1
#
