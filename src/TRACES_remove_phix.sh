#!/bin/bash
#
# $1 -> NUMBER OF THREADS
#
# IT ASSUMES THAT THE FOLLOWING OUTPUT FILES EXIST:
# o_fw_pr.fq o_fw_unpr.fq o_rv_pr.fq o_rv_unpr.fq
#
## REMOVE PHIX FROM THE SAMPLES
MAGNET -v -F -n $1 -m 6:1:1:0/0 -m 13:20:1:3/10 -t 0.9 -i -o NP-o_fw_pr.fq F_PHIX.fa o_fw_pr.fq
MAGNET -v -F -n $1 -m 6:1:1:0/0 -m 13:20:1:3/10 -t 0.9 -i -o NP-o_rv_pr.fq F_PHIX.fa o_rv_pr.fq
MAGNET -v -F -n $1 -m 6:1:1:0/0 -m 13:20:1:3/10 -t 0.9 -i -o NP-o_fw_unpr.fq F_PHIX.fa o_fw_unpr.fq
MAGNET -v -F -n $1 -m 6:1:1:0/0 -m 13:20:1:3/10 -t 0.9 -i -o NP-o_rv_unpr.fq F_PHIX.fa o_rv_unpr.fq
#
