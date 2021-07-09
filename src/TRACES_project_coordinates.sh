#!/bin/bash
#
# TRACES_project_coordinates.sh positions_files max_coverage
#
if [ ! -f "$1" ]; then
  echo "ERROR: Input file $1 does not exists!";
  exit;
fi
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
  for((x = $IPOS ; x < $EPOS ; ++x));
    do
    if [[ "$2" -ne "0" ]];
      then
      COVE2=`echo $COVE | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' | bc -l`;
      COVE_PARSED=`printf "%.0f" "$COVE2"`;
      if [[ "$COVE_PARSED" -gt "$2" ]];
        then
        printf "%u\t%u\n" "$x" `echo $2`;
        else
        printf "%u\t%u\n" "$x" "$COVE_PARSED"     
        fi
      else
      printf "%u\t%u\n" "$x" "$COVE" 
      fi
    done
  done
#
