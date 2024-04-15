#!/bin/bash
#
# $1 -> ORGAN NAME
# $2 -> NUMBER OF THREADS
#
# IT ASSUMES THAT THE FOLLOWING INPUT FILES EXIST:
# o_fw_pr.fq o_rv_pr.fq o_fw_unpr.fq o_rv_unpr.fq
#
# ASSEMBLE
if [ -s o_fw_unpr.fq ] && [ -s o_rv_unpr.fq ]; then
  spades.py --meta --threads $2 -o ../output_data/TRACES_denovo_$1 -1 o_fw_pr.fq -2 o_rv_pr.fq -s o_fw_unpr.fq -s o_rv_unpr.fq
elif [ -s o_fw_unpr.fq ]; then
  spades.py --meta --threads $2 -o ../output_data/TRACES_denovo_$1 -1 o_fw_pr.fq -2 o_rv_pr.fq -s o_fw_unpr.fq
elif [ -s o_rv_unpr.fq ]; then
  spades.py --meta --threads $2 -o ../output_data/TRACES_denovo_$1 -1 o_fw_pr.fq -2 o_rv_pr.fq -s o_rv_unpr.fq
else
  spades.py --meta --threads $2 -o ../output_data/TRACES_denovo_$1 -1 o_fw_pr.fq -2 o_rv_pr.fq
fi
#
