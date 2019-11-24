#!/bin/bash
#
declare -a INDIVIDUALS=('1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' '23' '24' '25' '26' '27' '28' '29' '30' '31' '32' '33' '34' '35' '36' '37' '38' '40' '41' '42' '43' '44' '45' '46' '47' '48' '49' '50' '51' '52' '53' '54' 'POOL');
#
for x in "${INDIVIDUALS[@]}"
  do
  echo "Building I$x ...";
  cp tracespipe_template/src/TRACES_assemble_all.sh novaseq_2/I$x/tracespipe_template/src 
  echo "Done!";
  done
