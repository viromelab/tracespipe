#!/bin/bash
#
# ==============================================================================
#
rm -f alg_tot.fasta alg_pa_tot.fasta alg_tot.fasta alg_pa_tot.fasta
#
wget https://www.hmtdb.uniba.it/download_file?dataset=alg_tot.zip -O alg_tot.zip
wget https://www.hmtdb.uniba.it/download_file?dataset=alg_pa_tot.zip -O alg_pa_tot.zip
#
unzip alg_tot.zip
unzip alg_pa_tot.zip
#
# ==============================================================================
#
