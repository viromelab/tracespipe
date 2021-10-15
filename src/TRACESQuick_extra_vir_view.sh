#!/bin/bash
#
cat $1 | grep -v "erpes" | grep -v "papilloma" | grep -v "erythroparvovirus" \
       | grep -v "parvovirus" | grep -v "erythrovirus" \
       | grep -v "endogenous_retrovirus" | grep -v "_B19_" | grep -v "_B19V_" \
       | grep -v "JC_virus" | grep -v "polyomavirus" | grep -v "Hepatitis_B" \
       | grep -v "_HBV_" | grep -v "WWII_" | grep -v "GT3A_" | grep -v "GT3B_" \
       | grep -v "GT2_" | grep -v "GT1_" | grep -v "Torque_teno_virus" \
       | awk '{ if($3 > 5.0) print $1"\t"$2"\t"$3"\t"$4; }'
#
