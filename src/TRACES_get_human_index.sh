#!/bin/bash
#
# ==============================================================================
# GET HUMAN HOST INDEXES
#
rm -f hg19.zip
wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip
unzip hg19.zip
./make_hg19.sh
#
# ==============================================================================
