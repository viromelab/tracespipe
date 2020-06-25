#!/bin/bash
#
# ==============================================================================
#
rm -f mitochondrion.1.1.genomic.fna.gz mitochondrion.2.1.genomic.fna.gz MTs.fa
#
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.1.genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.2.1.genomic.fna.gz
zcat mitochondrion.1.1.genomic.fna.gz mitochondrion.2.1.genomic.fna.gz > SPECIES-MTs.fa
#
# ==============================================================================
#
