#!/bin/bash
#
# IT ASSUMES THAT THE FOLLOWING OUTPUT FILES EXIST:
# NP-o_fw_pr.fq NP-o_fw_unpr.fq NP-o_rv_pr.fq NP-o_rv_unpr.fq
#
## EXTRACT MTDNA FROM THE SAMPLES
MAGNET -v -F -n 12 -m 6:1:1:0/0 -m 13:20:1:3/10 -t 0.3 -o CY-o_fw_pr.fq cy.fa NP-o_fw_pr.fq
MAGNET -v -F -n 12 -m 6:1:1:0/0 -m 13:20:1:3/10 -t 0.3 -o CY-o_rv_pr.fq cy.fa NP-o_rv_pr.fq
MAGNET -v -F -n 12 -m 6:1:1:0/0 -m 13:20:1:3/10 -t 0.3 -o CY-o_fw_unpr.fq cy.fa NP-o_fw_unpr.fq
MAGNET -v -F -n 12 -m 6:1:1:0/0 -m 13:20:1:3/10 -t 0.3 -o CY-o_rv_unpr.fq cy.fa NP-o_rv_unpr.fq
#
cat CY-o_fw_pr.fq CY-o_rv_pr.fq CY-o_fw_unpr.fq CY-o_rv_unpr.fq > REFERENCE_CY.fq
#
gto_geco -v -rm 13:20:1:0:0.94/0:0:0 -rm 16:200:1:20:0.95/5:10:0.95 -r REFERENCE_CY.fq cy.fa &> REP_$1_cy.txt
#
