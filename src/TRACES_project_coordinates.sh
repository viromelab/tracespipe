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
        printf "%u\t%.0f\n" "$x" `echo $2 | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' | bc -l`;
        else
        printf "%u\t%.0f\n" "$x" `echo $COVE | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' | bc -l`;
        fi
      else
      printf "%u\t%.0f\n" "$x" `echo $COVE | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' | bc -l`;
      fi
    done
  done
#
