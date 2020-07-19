#!/bin/bash
# 
# ./TRACES_blastn_n_db.sh query.fa
#
# ==============================================================================
#
SPATH="../system_files/blast_dbs/nt";
QUERY="$1";
DBASE="nt";
#
#STYLE="6 qacc sacc pident mismatch btop";
#
# qseqid    Query Seq-id
# sseqid    Subject Seq-id
# qstart    Start of alignment in query
# qend      End of alignment in query
# sstart    Start of alignment in subject
# send      End of alignment in subject
# evalue    Expect value
# bitscore  Bit score
# length    Alignment length
# pident    Percentage of identical matches
# mismatch  Number of mismatches
# gapopen   Number of gap openings
# btop      Blast traceback operations (BTOP)
#
STYLE="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore btop";
#
# RUN ==========================================================================
#
blastn -db $SPATH/$DBASE -query $QUERY -outfmt "$STYLE";
#
# ==============================================================================
