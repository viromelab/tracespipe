#!/bin/bash
#
rm -f o_fw_pr.fq.gz o_fw_unpr.fq.gz o_rv_pr.fq.gz o_rv_unpr.fq.gz;
#
trimmomatic PE -threads 12 -phred33 FW_READS.fq.gz RV_READS.fq.gz o_fw_pr.fq o_fw_unpr.fq o_rv_pr.fq o_rv_unpr.fq ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35
#
