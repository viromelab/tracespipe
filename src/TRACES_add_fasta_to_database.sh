#!/bin/bash
#
# check if $2
FIRST_SYMBOL=`head -c 1 $2`;
#echo "First symbol: $FIRST_SYMBOL";
#
if [ "$FIRST_SYMBOL" = ">" ]
  then 
  echo "Adding the FASTA file to the $1 database";
  cat $1 $2 > tmp-add-fasta-$1
  mv tmp-add-fasta-$1 $1
  echo "Done!"
  else
  echo "ERROR: the file $2 is not a FASTA file!";
  exit 1;
  fi
#
