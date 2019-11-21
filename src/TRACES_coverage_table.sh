#!/bin/bash
#
declare -a VIRUSES=('B19' 'HV1' 'HV2' 'HV3' 'HV4' 'HV5' 'HV6' 'HV6A' 'HV6B' 'HV7' 'HV8' 'POLY1' 'POLY2' 'POLY3' 'POLY4' 'POLY5' 'POLY6' 'POLY7' 'POLY8' 'POLY9' 'POLY10' 'POLY11' 'POLY12' 'POLY13' 'POLY14' 'TTV' 'HBOV1' 'HBOVNOT1' 'HBV' 'HPV' 'VARV' 'SV40' 'CUTA' 'HERV');
declare -a ORGANS=( $(cat ../meta_data/meta_info.txt | tr ':' '\t' | awk '{ print $1;}') );
#
printf "; ";
for organ in "${ORGANS[@]}" #
  do
  printf "%s;;" "$organ";
  done
printf "\n";
#
printf ";";
for organ in "${ORGANS[@]}" #
  do
  printf "D;H;";
  done
printf "\n";
#
x=0;
for virus in "${VIRUSES[@]}" #
  do
  y=0;
  printf "\\multirow{2}{*}{%s}  & \\cellcolor[gray]{.9} - &\\cellcolor[gray]{.9} - & \\cellcolor[gray]{.9} - &\\cellcolor[gray]{.9} - & \\cellcolor[gray]{.9} - &\\cellcolor[gray]{.9} - &\\cellcolor[gray]{.9} - &\\cellcolor[gray]{.9} - &\\cellcolor[gray]{.9} - &\\cellcolor[gray]{.9} - &\\cellcolor[gray]{.9} - &\\cellcolor[gray]{.9} - &\\cellcolor[gray]{.9} - &\\cellcolor[gray]{.9} - &\\cellcolor[gray]{.9} - &\\cellcolor[gray]{.9} - &\\cellcolor[gray]{.9} - &\\cellcolor[gray]{.9} - &\\cellcolor[gray]{.9} - & \\cellcolor[gray]{.9} - &\\cellcolor[gray]{.9} - &\\cellcolor[gray]{.9} -  \\\\\ " "$virus";
  for organ in "${ORGANS[@]}" #
    do
    #
    if [ ! -f ../output_data/TRACES_viral_statistics/$virus-total-depth-coverage-$organ.txt ];
      then
      D_coverage="0";
      else
      D_coverage=`cat ../output_data/TRACES_viral_statistics/$virus-total-depth-coverage-$organ.txt | awk '{ print $3;}'`;
      fi
    #
    if [ ! -f ../output_data/TRACES_viral_statistics/$virus-total-horizontal-coverage-$organ.txt ];
      then
      H_coverage="0";
      else
      H_coverage=`cat ../output_data/TRACES_viral_statistics/$virus-total-horizontal-coverage-$organ.txt | awk '{ print $3;}'`;
      fi
    #
    if [[ "$H_coverage" == "" ]];
      then
      H_coverage="0";
      fi
    #
    if [[ "$D_coverage" == "" ]];
      then
      D_coverage="0";
      fi
    #
    LC_NUMERIC="en_US.UTF-8" printf "%2.1f & %2.0f\t" "$D_coverage" "$H_coverage";
    ((y++));
    done
  printf "\\\\\ \\hline\n";
  ((x++));
  done
#

