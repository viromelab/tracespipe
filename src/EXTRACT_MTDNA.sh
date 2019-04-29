#!/bin/bash
#
## EXTRACT MTDNA FROM THE SAMPLES
echo "RUNNING MAGNET ...";
MAGNET -v -F -n 12 -m 6:1:1:0/0 -m 13:20:1:3/10 -t 0.2 -o MT-o_fw_pr.fq mtDNA.fa o_fw_pr.fq
MAGNET -v -F -n 12 -m 6:1:1:0/0 -m 13:20:1:3/10 -t 0.2 -o MT-o_rv_pr.fq mtDNA.fa o_rv_pr.fq
MAGNET -v -F -n 12 -m 6:1:1:0/0 -m 13:20:1:3/10 -t 0.2 -o MT-o_fw_unpr.fq mtDNA.fa o_fw_unpr.fq
MAGNET -v -F -n 12 -m 6:1:1:0/0 -m 13:20:1:3/10 -t 0.2 -o MT-o_rv_unpr.fq mtDNA.fa o_rv_unpr.fq
#
