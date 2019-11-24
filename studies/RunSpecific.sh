#!/bin/bash
#
# RUNNING SPECIFIC REFERENCE FOR MULTIPLE INDIVIDUALS AND SEQUENCED SAMPLES
#
declare -a INDIVIDUALS1=('I1' 'I3' 'I6' 'I7' 'I8' 'I14' 'I15' 'I16' 'I18')
declare -a INDIVIDUALS2=('I23' 'I24' 'I25' 'I26' 'I27' 'I28' 'I36' 'I40' 'I46' 'I47' 'I54')
#
# HISEQ
echo "Running hiSeq ..."
cd hiseq/tracespipe/src/
./TRACESPipe.sh --run-specific AB550331.1;
cd ../../../
echo "Done!";
#
# NOVASEQ 1
for x in "${INDIVIDUALS1[@]}"
  do
  echo "Running $x ..."
  cd novaseq_1/$x/tracespipe_template/src/
  ./TRACESPipe.sh --run-specific AB550331.1
  cd ../../../../
  echo "Done!";
  done
#
# NOVASEQ 2
for x in "${INDIVIDUALS2[@]}"
  do
  echo "Running $x ..."
  cd novaseq_2/$x/tracespipe_template/src/
  ./TRACESPipe.sh --run-specific AB550331.1
  cd ../../../../
  echo "Done!";
  done
#  
