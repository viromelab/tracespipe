#!/bin/bash
#
declare -a INDIVIDUALS=('1' '3' '4' '6' '7' '8' '12' '14' '15' '16' '18' '19' '29' '30');
#
for x in "${INDIVIDUALS[@]}"
  do
  echo "Running I$x ...";
  cd novaseq_1_no_duplications/I$x/tracespipe_template/src/ 
  ./TRACESPipe.sh --run-meta --run-all-v-alig --remove-dup
  cd ../../../../
  echo "Done!";
  done
#
