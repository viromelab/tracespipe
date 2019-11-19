#!/bin/bash
#
mapfile -t POSITIONS < $1
#
for read in "${POSITIONS[@]}" #
  do
  #
  IPOS=`echo $read | awk '{ print $2 }'`;
  EPOS=`echo $read | awk '{ print $3 }'`;
  COVE=`echo $read | awk '{ print $4 }'`;
  #
  for(( x = $IPOS ; x < $EPOS ; ++x ));
    do
    printf "%u\t%u\n" "$x" "$COVE";
    done
  done
