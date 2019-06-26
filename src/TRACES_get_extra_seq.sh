#!/bin/bash
#
gto_fasta_extract_read_by_pattern -p "$2" < $1 > tmp_exists_or_not.fa
FIRST_LINE=`head -n 1 tmp_exists_or_not.fa`;
NUMBER_OF_LINES=`wc -l tmp_exists_or_not.fa |awk '{print $1;}'`;
if [ $NUMBER_OF_LINES -ge "2" ]
then
  echo -e "\e[31mERROR: This sequence seems to already exist in the VDB.fa.\e[0m"
  echo "Header: $FIRST_LINE";
  exit 1;
else
  rm -f en_seqs.fa;
  efetch -db nucleotide -format fasta -id $2 > en_seqs.fa
  NUMBER_OF_DW_LINES=`wc -l en_seqs.fa |awk '{print $1;}'`;
  if [ $NUMBER_OF_DW_LINES -ge "2" ]
    then
    echo "Adding sequence ...";
    cat $1 en_seqs.fa > TMP_VDB.fa.tmp
    mv TMP_VDB.fa.tmp $1
    echo "Done!";
    fi
fi
#
