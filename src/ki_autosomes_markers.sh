#!/bin/bash
# 
# ==============================================================================
# |  POPULATION INFERENCE FROM THE ANALYSIS OF SPECIFIC PATTERNS IN AUTOSOMES  |
# ==============================================================================
#
# REQUIREMENTS:
#   [x] autosome_snps_data.csv
#
# INSTALLATION:
#   [x] conda install -c bioconda bwa2
#   [x] conda install -c cobilab gto
#
# ==============================================================================
#
echo "Searching for specific SNPs in Autosomes ..." ;
#
# IT ADMITS A FILE WITH: 
# GID \t POSITIONS \t SNP_PATTERN \t POPULATION
mapfile -t SNIPS < autosome_snps_data.csv
#
FLANK=50;
INIT_POS=0;
#
for entry in "${SNIPS[@]}"
  do
  #
  GID=`echo $entry | awk '{ print $1; }'`;
  SNP_POSITION=`echo $entry | awk '{ print $2; }'`;
  SNP_PATTERN=`echo $entry | awk '{ print $3; }'`;
  POPULATION=`echo $entry | awk '{ print $4; }'`;
  #
  if [[ $SNP_POSITION -gt $FLANK ]];
    then
    INIT_POS=`echo "$SNP_POSITION-$FLANK" | bc -l`;
    else
    INIT_POS=0;
    echo "WARNING: flanking position below the limits! ($SNP_POSITION)";
    fi
  #
  END_POS=`echo "$SNP_POSITION+$FLANK" | bc -l`;
  #
  efetch -db nucleotide -id "$GID" -seq_start "$INIT_POS" -seq_stop "$END_POS" -format fasta
  #
  done
#
echo "Done";
##
# ==============================================================================
