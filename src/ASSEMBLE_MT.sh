#!/bin/bash
#
spades.py -t 12 -m 20 --careful -o out_spades_$1 -1 MT-o_fw_pr.fq -2 MT-o_rv_pr.fq -s MT-o_fw_unpr.fq -s MT-o_rv_unpr.fq
#
