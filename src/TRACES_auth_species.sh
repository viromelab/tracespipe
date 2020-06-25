#!/bin/bash
#
THREADS="$1";
OUTPUT="../output_data/TRACES_mtdna_authentication";
#
# ==============================================================================
#
CHECK_FILE () {
  if [ ! -f $1 ];
    then
    echo -e "\e[31mERROR: file $1 not found!\e[0m"
    echo "For addition information, see the instructions at the web page."
    exit 1;
    fi
  }
#
# ==============================================================================
#
CHECK_FILE "../meta_data/meta_info.txt";
CHECK_FILE "SPECIES-MTs.fa";
#
mkdir -p $OUTPUT;
#
mapfile -t READS < ../meta_data/meta_info.txt
#
for read in "${READS[@]}" #
  do
  #
  ORGAN=`echo $read | tr ':' '\t' | awk '{ print $1 }'`;
  SPL_Forward=`echo $read | tr ':' '\t' | awk '{ print $2 }'`;
  SPL_Reverse=`echo $read | tr ':' '\t' | awk '{ print $3 }'`;
  #
  echo -e "\e[34m[TRACESPipe]\e[32m Running $ORGAN sample ...\e[0m";
  zcat ../input_data/$SPL_Forward ../input_data/$SPL_Reverse > reads_human_auth.fq
  FALCON -F -c 30 -l 47 -t 200 -x $OUTPUT/MT-SPECIES-$ORGAN.txt -n $THREADS reads_human_auth.fq SPECIES-MTs.fa 1>> ../logs/Log-stdout-$ORGAN.txt 2>> ../logs/Log-stderr-$ORGAN.txt;
  done
#
# ==============================================================================
#
