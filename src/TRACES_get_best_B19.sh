#!/bin/bash
ORGAN=$1;
#
RESULT=`cat top-$ORGAN.csv \
| grep -a -e "irus B19 " -e "irus_B19_" -e " B19V " -e "B19 virus" -e "B19_virus" -e "parvovirus B19" -e "Parvovirus B19" -e "erythrovirus V9" -e "Erythrovirus V9" -e "AY386330.1" -e "AY661661" -e "AY028237.1" -e "AY504945.1" -e "KM393163.1" -e "KM393164.1" -e "KM393165.1" -e "KM393166.1" -e "KM393167.1" -e "KM393168.1" -e "KM393169.1" -e "DQ357065.1" -e "DQ408301.1" -e "DQ225150.1" -e "DQ225149.1" -e "DQ225148.1" -e "AJ781038.1" -e "AJ781036.1" -e "AJ781037.1" -e "AJ781035.1" -e "AJ781034.1" -e "AJ781033.1" -e "AJ781032.1" -e "AJ781031.1" -e "FN669502.1" -e "FN598217.1" -e "FN669503.1" -e "FN669507.1" -e "FN669504.1" -e "FN669505.1" -e "FN598218.1" -e "Z68146.1" -e "Z70528.1" -e "Z70560.1" -e "Z70599.1" -e "M24682.1" -e "M13178.1" -e "AB030673.1" -e "AB030694.1" -e "AB126262.1" -e "AB126263.1" -e "AB126264.1" -e "AB126266.1" -e "AB126267.1" -e "AB126268.1" -e "AB126269.1" -e "AB126271.1" -e "AB030693.1" -e "GM703961.1" -e "GM703832.1" -e "KY940273.1" -e "KT268312" -e "KT310174" -e "MK989716.1" -e "MH201455" -e "MH201456" -e "GP700454.1" -e "AF113323.1" -e "AF161226.1" -e "AF161225.1" -e "AF161224.1" -e "AF162273.1" -e "FJ591158.1." -e "FB715682.1" -e "000883.1" -e "DQ293995.2" -e "MT410185" -e "MT410187" -e "MT410189" -e "AB550331.1" -e "HQ340601" -e "HQ340602" -e "AJ717293.1" -e "LN680968.2" -e "AY044266.2" -e "AY064476.1" -e "AY903437.1" -e "KF724386" -e "KF724387" -e "EF216869.1" -e "DQ333426.1" -e "DQ333427.1" -e "DQ333428.1" -e "MT410184" -e "MT410186" -e "MT410188" -e "MT410190" -e "AJ249437.1" -e "AX003421.1" -e "AY582125.2" -e "DQ234775.2" -e "DQ234771.2" -e "DQ234769.2" -e "FJ265736.1" -e "AY582124.2" -e "AY647977.1" -e "DQ234778.2" -e "DQ234779.2" -e "DQ408305.1" -e "DQ408302.1" -e "DQ408304.1" -e "DQ408303.1" -e "AY083234.1" \
| grep -a -e "omplete genome" -e "omplete_genome" \
| grep -v "Polyomavirus" \
| grep -v "polyomavirus" \
| grep -v "bocavirus" \
| grep -v "Bocavirus" \
| grep -v "bocaparvovirus" \
| grep -v "Bocaparvovirus" \
| awk '{ if($3 > 0 && $2 > 5300 && $2 < 5900) print $3"\t"$4; }' \
| head -n 1 \
| awk '{ print $1"\t"$2;}' \
| sed "s/NC\_/NC-/" \
| tr '_' '\t' \
| awk '{ print $1"\t"$2;}'`;
if [ -z "$RESULT" ]
  then
  echo -e "-\t-";
  else
  echo "$RESULT" | sed "s/NC-/NC\_/"
  fi
