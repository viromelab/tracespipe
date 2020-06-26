#!/bin/bash
#
THREADS="$3";
OUTPUT="../output_data/TRACES_specific_diff";
#
# ==============================================================================
#
CHECK_FILE () {
  if [ ! -f $1 ];
    then
    echo -e "\e[34m[TRACESPipe]\e[31m ERROR: file $1 not found!\e[0m"
    echo -e "\e[34m[TRACESPipe]\e[33m For addition information, see the instructions at the web page.\e[0m"
    exit 1;
    fi
  }
#
#
CHECK_FILE_REC () {
  if [ ! -f $1 ];
    then
    echo -e "\e[34m[TRACESPipe]\e[31m ERROR: file $1 not found!\e[0m"
    echo -e "\e[34m[TRACESPipe]\e[32m TIP: this is possibly given by absence of reconstructed data.\e[0m"
    echo -e "\e[34m[TRACESPipe]\e[33m For addition information, see the instructions at the web page.\e[0m"
    exit 1;
    fi
  }
#
# ==============================================================================
#
gto_fasta_extract_read_by_pattern -p "$2" < VDB.fa | awk "/^>/ {n++} n>1 {exit} 1" > SPECIFIC-$2.fa 2>> ../logs/Log-system.txt;
#
CHECK_FILE "../meta_data/meta_info.txt";
CHECK_FILE "SPECIFIC-$2.fa";
#
NLINES=`grep ">" SPECIFIC-$2.fa | wc -l`;
if [[ "$NLINES" != "1" ]];
  then
	  echo -e "\e[34m[TRACESPipe]\e[31m ERROR: pattern or GID ($2) not found in VDB.fa!\e[0m";
  exit 1;
  fi
#
mkdir -p $OUTPUT;
#
mapfile -t READS < ../meta_data/meta_info.txt
#
echo -e "Organ\tVirus\tGID\tAligned bases\tidentity\tSNPs" > ../output_data/TRACES_specific_diff/$1-$2_Diff.txt;
#
for read in "${READS[@]}" #
  do
  #
  ORGAN=`echo $read | tr ':' '\t' | awk '{ print $1 }'`;
  SPL_Forward=`echo $read | tr ':' '\t' | awk '{ print $2 }'`;
  SPL_Reverse=`echo $read | tr ':' '\t' | awk '{ print $3 }'`;
  #
  echo -e "\e[34m[TRACESPipe]\e[32m Running $ORGAN sample ...\e[0m";
  #
  CHECK_FILE_REC "../output_data/TRACES_hybrid_R5_consensus/$1-consensus-$ORGAN.fa";
  #
  cp ../output_data/TRACES_hybrid_R5_consensus/$1-consensus-$ORGAN.fa G_A.fa;
  cp SPECIFIC-$2.fa G_B.fa;
  dnadiff G_A.fa G_B.fa 1>> ../logs/Log-stdout-$ORGAN.txt 2>> ../logs/Log-stderr-$ORGAN.txt;
  IDEN=`cat out.report | grep "AvgIdentity " | head -n 1 | awk '{ print $2;}'`;
  ALBA=`cat out.report | grep "AlignedBases " | head -n 1 | awk '{ print $2;}'`;
  SNPS=`cat out.report | grep TotalSNPs | awk '{ print $2;}'`;
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$ORGAN" "$1" "$2" "$ALBA" "$IDEN" "$SNPS" 1>> ../output_data/TRACES_specific_diff/$1-$2_Diff.txt
  rm -f G_A.fa G_B.fa ;
  #
  done
  #
  cat ../output_data/TRACES_specific_diff/$1-$2_Diff.txt;
  #
#
# ==============================================================================
#
