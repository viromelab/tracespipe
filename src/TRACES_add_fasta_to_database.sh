#!/bin/bash
#
# check if $2
FIRST_SYMBOL=`head -c 1 $2`;
if [ $FIRST_SYMBOL == '>' ]
  then 
  echo "It is a fasta file";
  else
  echo "ERROR: the file $2 is not a FASTA file!";
  fi
#
#cat $1 $2 > tmp-add-fasta-$1
#mv tmp-add-fasta-$1 $1
#
