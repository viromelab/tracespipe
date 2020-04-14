#!/bin/bash
#
# TRACES_project_coordinates.sh positions_files max_coverage
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
    if [[ "$2" -ne "0" ]];
      then
      if [[ "$COVE" -gt "$2" ]];
        then
        printf "%u\t%u\n" "$x" "$2";
        else
        printf "%u\t%u\n" "$x" "$COVE";
        fi
      else
      printf "%u\t%u\n" "$x" "$COVE";
      fi
    done
  done
#
