#!/bin/bash
#
EMAIL="$1";
TOP_SIZE=500;
#
declare -a INDIVIDUALS1=('I1' 'I3' 'I4' 'I6' 'I7' 'I8' 'I12' 'I14' 'I15' 'I16' 'I18' 'I19' 'I29' 'I30')
declare -a INDIVIDUALS2=('I1' 'I2' 'I3' 'I4' 'I5' 'I6' 'I7' 'I8' 'I9' 'I10' 'I11' 'I12' 'I13' 'I14' 'I15' 'I16' 'I17' 'I18' 'I19' 'I20' 'I21' 'I22' 'I23' 'I24' 'I25' 'I26' 'I27' 'I28' 'I29' 'I30' 'I31' 'I32' 'I33' 'I34' 'I35' 'I36' 'I37' 'I38' 'I40' 'I41' 'I42' 'I43' 'I44' 'I45' 'I46' 'I47' 'I48' 'I49' 'I50' 'I51' 'I52' 'I53' 'I54' 'IPOOL')
#
# HISEQ
echo "" > hiseq_coverage.csv;
echo "" > hiseq_tops.csv;
echo "" > hiseq_best.csv;
echo "Running hiSeq ..."
echo "HiSeq:" >> hiseq_coverage.csv
echo "HiSeq:" >> hiseq_tops.csv
echo "HiSeq:" >> hiseq_best.csv
cd hiseq/tracespipe/src/
./TRACESPipe.sh --coverage-csv >> ../../../hiseq_coverage.csv;
cat ../meta_data/meta_info.txt | tr ':' '\t' | awk '{ print $1}' > names_tmp_TRACES_fetch.tmp
mapfile -t NAMES < names_tmp_TRACES_fetch.tmp
for name in "${NAMES[@]}" #
  do
  echo "Organ: $name" >> ../../../hiseq_tops.csv;
  head -n $TOP_SIZE ../output_data/TRACES_results/top-$name.csv >> ../../../hiseq_tops.csv;
  done
rm -f names_tmp_TRACES_fetch.tmp
#  
cat ../output_data/TRACES_results/REPORT_META_VIRAL_ALL_SAMPLES.txt >> ../../../hiseq_best.csv;
echo "" >> ../../../hiseq_coverage.csv
echo "" >> ../../../hiseq_tops.csv
echo "" >> ../../../hiseq_best.csv
cd ../../../
echo "Done!";
#
# NOVASEQ 1
echo "" > novaseq1_coverage.csv;
echo "" > novaseq1_tops.csv;
echo "" > novaseq1_best.csv;
for x in "${INDIVIDUALS1[@]}"
  do
  echo "Running $x ..."
  echo "$x:" >> novaseq1_coverage.csv
  echo "$x:" >> novaseq1_tops.csv
  echo "$x:" >> novaseq1_best.csv
  cd novaseq_1/$x/tracespipe_template/src/
  ./TRACESPipe.sh --coverage-csv >> ../../../../novaseq1_coverage.csv; 
  #
  cat ../meta_data/meta_info.txt | tr ':' '\t' | awk '{ print $1}' > names_tmp_TRACES_fetch.tmp
  mapfile -t NAMES < names_tmp_TRACES_fetch.tmp
  for name in "${NAMES[@]}" #
    do
    echo "Organ: $name" >> ../../../../novaseq1_tops.csv;
    head -n $TOP_SIZE ../output_data/TRACES_results/top-$name.csv >> ../../../../novaseq1_tops.csv;
    done
  rm -f names_tmp_TRACES_fetch.tmp
  #  
  cat ../output_data/TRACES_results/REPORT_META_VIRAL_ALL_SAMPLES.txt >> ../../../../novaseq1_best.csv;
  echo "" >> ../../../../novaseq1_coverage.csv
  echo "" >> ../../../../novaseq1_tops.csv
  echo "" >> ../../../../novaseq1_best.csv
  cd ../../../../
  echo "Done!";
  done
#
# NOVASEQ 2
echo "" > novaseq2_coverage.csv;
echo "" > novaseq2_tops.csv;
echo "" > novaseq2_best.csv;
for x in "${INDIVIDUALS2[@]}"
  do
  echo "Running $x ..."
  echo "$x:" >> novaseq2_coverage.csv
  echo "$x:" >> novaseq2_tops.csv
  echo "$x:" >> novaseq2_best.csv
  cd novaseq_2/$x/tracespipe_template/src/
  ./TRACESPipe.sh --coverage-csv >> ../../../../novaseq2_coverage.csv; 
  #
  cat ../meta_data/meta_info.txt | tr ':' '\t' | awk '{ print $1}' > names_tmp_TRACES_fetch.tmp
  mapfile -t NAMES < names_tmp_TRACES_fetch.tmp
  for name in "${NAMES[@]}" #
    do
    echo "Organ: $name" >> ../../../../novaseq2_tops.csv;
    head -n $TOP_SIZE ../output_data/TRACES_results/top-$name.csv >> ../../../../novaseq2_tops.csv;
    done
  rm -f names_tmp_TRACES_fetch.tmp
  #  
  cat ../output_data/TRACES_results/REPORT_META_VIRAL_ALL_SAMPLES.txt >> ../../../../novaseq2_best.csv;
  echo "" >> ../../../../novaseq2_coverage.csv
  echo "" >> ../../../../novaseq2_tops.csv
  echo "" >> ../../../../novaseq2_best.csv
  cd ../../../../
  echo "Done!";
  done
#  
# EMAIL RESULTS:
mail -s "hiseq tops" $EMAIL < hiseq_tops.csv; 
mail -s "novaseq1 tops" $EMAIL < novaseq1_tops.csv; 
mail -s "novaseq2 tops" $EMAIL < novaseq2_tops.csv; 
#
mail -s "hiseq best" $EMAIL < hiseq_best.csv; 
mail -s "novaseq1 best" $EMAIL < novaseq1_best.csv; 
mail -s "novaseq2 best" $EMAIL < novaseq2_best.csv; 
#
mail -s "hiseq coverage" $EMAIL < hiseq_coverage.csv; 
mail -s "novaseq1 coverage" $EMAIL < novaseq1_coverage.csv; 
mail -s "novaseq2 coverage" $EMAIL < novaseq2_coverage.csv; 
#
