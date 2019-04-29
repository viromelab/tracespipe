#!/bin/bash
#
gunzip o_fw_pr.fq.gz 
gunzip o_rv_pr.fq.gz 
gunzip o_fw_unpr.fq.gz 
gunzip o_rv_unpr.fq.gz
#
## REMOVE PHIX FROM THE SAMPLES
echo "RUNNING MAGNET ...";
MAGNET -v -F -n 12 -m 6:1:1:0/0 -m 13:20:1:3/10 -t 0.9 -i -o NP-o_fw_pr.fq F_PHIX.fa o_fw_pr.fq
MAGNET -v -F -n 12 -m 6:1:1:0/0 -m 13:20:1:3/10 -t 0.9 -i -o NP-o_rv_pr.fq F_PHIX.fa o_rv_pr.fq
MAGNET -v -F -n 12 -m 6:1:1:0/0 -m 13:20:1:3/10 -t 0.9 -i -o NP-o_fw_unpr.fq F_PHIX.fa o_fw_unpr.fq
MAGNET -v -F -n 12 -m 6:1:1:0/0 -m 13:20:1:3/10 -t 0.9 -i -o NP-o_rv_unpr.fq F_PHIX.fa o_rv_unpr.fq
#
