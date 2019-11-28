#!/bin/bash
#
declare -a INDIVIDUALS1=('I1' 'I3' 'I4' 'I6' 'I7' 'I8' 'I12' 'I14' 'I15' 'I16' 'I18' 'I19' 'I29' 'I30')
declare -a INDIVIDUALS2=('I1' 'I2' 'I3' 'I4' 'I5' 'I6' 'I7' 'I8' 'I9' 'I10' 'I11' 'I12' 'I13' 'I14' 'I15' 'I16' 'I17' 'I18' 'I19' 'I20' 'I21' 'I22' 'I23' 'I24' 'I25' 'I26' 'I27' 'I28' 'I29' 'I30' 'I31' 'I32' 'I33' 'I34' 'I35' 'I36' 'I37' 'I38' 'I40' 'I41' 'I42' 'I43' 'I44' 'I45' 'I46' 'I47' 'I48' 'I49' 'I50' 'I51' 'I52' 'I53' 'I54' 'IPOOL')
#
# HISEQ
echo "Updating hiseq ...";
cp tracespipe/src/*.sh hiseq_no_duplications/tracespipe/src/
echo "Done!";
# 
# NOVASE 1
for name in "${INDIVIDUALS1[@]}" #
  do
  echo "Updating: $name ...";
  cp tracespipe/src/*.sh novaseq_1_no_duplications/$name/tracespipe_template/src/
  echo "Done!";
  done
#
# NOVASE 2
for name in "${INDIVIDUALS2[@]}" #
  do
  echo "Updating: $name ...";
  cp tracespipe/src/*.sh novaseq_2_no_duplications/$name/tracespipe_template/src/
  echo "Done!";
  done
#
