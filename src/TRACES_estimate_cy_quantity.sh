#!/bin/bash
#
ORGAN="$1";
THREADS="$2";
#
# IT ASSUMES THAT THE FOLLOWING OUTPUT FILES EXIST:
# o_fw_pr.fq o_fw_unpr.fq o_rv_pr.fq o_rv_unpr.fq
#
mkdir -p ../output_data/TRACES_results;
#
cat o_fw_pr.fq o_rv_pr.fq o_fw_unpr.fq o_rv_unpr.fq > READS_TO_CY.fq
#
MAGNET -v -F -n $THREADS -m 6:1:1:0/0 -m 13:20:1:5/1 -t 0.3 -o READS_CY_FILTERED_cy.fq cy.fa READS_TO_CY.fq
#
gto_geco -v -rm 8:1:1:0:0.9/0:0:0 -rm 13:20:1:0:0.94/0:0:0 -rm 16:200:1:15:0.95/5:10:0.95 -r READS_CY_FILTERED_cy.fq cy.fa > ../output_data/TRACES_results/REP_CY_$ORGAN.txt
#
rm -f READS_TO_CY.fq READS_CY_FILTERED_cy.fq;
#
