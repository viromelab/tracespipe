#!/bin/bash
USER="dpratas";
SERVER="atlas.fimm.fi";
PATH="/fas/NGS/pipes/fastq/fimm_ngs_mtoppinen/autopsy_NGS_novaseq1";
#
scp $USER@$SERVER:$PATH/\{V1_*,V2_*,V3_*,V4_*,V5_*,V6_*,V7_*,V8_*,V47_*,V48_*\} .
#
