#!/bin/bash
#
# IT ASSUMES THAT THE FOLLOWING OUTPUT FILES EXIST:
# o_fw_pr.fq o_fw_unpr.fq o_rv_pr.fq o_rv_unpr.fq
#
cat o_fw_pr.fq o_rv_pr.fq o_fw_unpr.fq o_rv_unpr.fq > REFERENCE_CY.fq
#
gto_geco -v -rm 8:1:1:0:0.9/0:0:0 -rm 14:20:1:0:0.94/0:0:0 -rm 16:200:1:30:0.95/5:10:0.95 -r REFERENCE_CY.fq cy.fa > REP_CY_$1.txt
#
rm -f REFERENCE_CY.fq;
#
