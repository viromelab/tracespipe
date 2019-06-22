#!/bin/bash
# 
# ==============================================================================
#  HAPLOGROUP INFERENCE THROW THE ANALYSIS OF SPECIFIC PATTERNS IN Y CHROMOSOME 
# ==============================================================================
#
# REQUIREMENTS:
#   [x] cy.fa [GRC38] needed : run ./ki_get_cy.sh
#   [x] cy_halotype_data.csv
#
# INSTALLATION:
#   [x] conda install -c bioconda bwa2
#   [x] conda install -c cobilab gto
#
# ==============================================================================
#
# ALIGN READS OF CY
echo "Searching for specific SNPs in Y chromosome ..." ;
./TRACES_cy_align_reads.sh cy.fa $1     #$1->ORGAN
./TRACES_cy_consensus.sh cy.fa cy_aligned_sorted-$1.bam $1
./TRACES_search_halotypes.sh $1
echo "Done";
##
# mapfile -t ENTRIES < cy_halotype_data.csv # http://www.phylotree.org/Y/marker_list.htm
##
#for read in "${ENTRIES[@]}" # 
#  do
#  #
#  POSITIONS=`echo $read | awk '{ print $5 }'`;
#  MUTATION=`echo $read | awk '{ print $6 }'`;
#  CLADE=`echo $read | awk '{ print $7 }'`;
#  #
#  echo "Searching for positions $POSITIONS with mutation $MUTATION belonging to major clade $CLADE ..." ;
#  gto_fasta_extract -i 333 -e 333 < in > out
#  
#
#  echo "Done";
#  #
#  done
#
# ==============================================================================
