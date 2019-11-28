#!/bin/bash
#
declare -a INDIVIDUALS=('1' '3' '4' '6' '7' '8' '12' '14' '15' '16' '18' '19' '29' '30');
#
for x in "${INDIVIDUALS[@]}"
  do
  echo "Running I$x ...";
  cd novaseq_1/I$x/tracespipe_template/src/ 
  ./TRACESPipe.sh --run-meta --run-mito --run-all-v-alig --run-hybrid
  cd ../../../../
  echo "Done!";
  done
#
