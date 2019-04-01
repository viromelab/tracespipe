#!/bin/bash
#
GET_SMASHPP=0;
GET_GTO=0;
GET_GOOSE=0;
BUILD_DATABASE=0;
#
FILTER_READS=0;
RUN_COMPARISON=1;
#
# ==============================================================================
# ================================= GET SMASHPP ================================
if [[ "$GET_SMASHPP" -eq "1" ]];
  then
  rm -fr smashpp/ Smashpp
  git clone https://github.com/smortezah/smashpp.git
  cd smashpp
  cmake .
  make
  cp smashpp ../Smashpp
  cd ..
  fi
#
# ==============================================================================
# ================================== GET GTO ===================================
if [[ "$GET_GTO" -eq "1" ]];
  then
  rm -fr gto/ gto_fasta_extract_read_by_pattern
  git clone https://github.com/bioinformatics-ua/gto.git
  cd gto/src/
  make
  cp ../bin/gto_fasta_extract_read_by_pattern ../../
  cp ../bin/gto_fasta_rand_extra_chars ../../
  cp ../bin/gto_fasta_to_seq ../../
  cp ../bin/gto_filter ../../
  cp ../bin/gto_reverse ../../
  cd ../../
  fi
#
# ==============================================================================
# ================================= GET GOOSE ==================================
if [[ "$GET_GOOSE" -eq "1" ]];
  then
  rm -fr goose/
  git clone https://github.com/pratas/goose.git
  cd goose/src/
  make
  cp goose-min ../../
  cd ../../
  fi
#
# ==============================================================================
# =============================== BUILD DATABASE ===============================
if [[ "$BUILD_DATABASE" -eq "1" ]];
  then
  rm -f assembly_summary.txt;
  #
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/assembly_summary.txt
  awk -F '\t' '{if($12=="Complete Genome") print $20}' assembly_summary.txt > ASCG.txt
  #
  # '{if($2=="Complete Genome"||$2=="Complete sequence") # IF ADDITIONAL FIELDS
  # (awk -F '\t' '{print $12}' assembly_summary.txt | sort -d | uniq)
  # KNOWN FIELDS: "assembly_level", "Chromosome", "Complete Genome", "Scaffold"
  #
  mkdir -p GB_DB_VIRAL/
  rm -f GB_DB_VIRAL/*.fa.gz
  #
  cat ASCG.txt | xargs -I{} -n1 -P8 wget -P GB_DB_VIRAL {}/*_genomic.fna.gz
  #
  # MOVE CDS AND RNA SEQUENCES TO RESPECTIVE DIRS
  mkdir -p GB_DB_VIRAL_CDS/
  mkdir -p GB_DB_VIRAL_RNA/
  rm -f GB_DB_VIRAL_CDS/*.fa.gz
  rm -f GB_DB_VIRAL_RNA/*.fa.gz
  mv GB_DB_VIRAL/*_cds_from_genomic.fna.gz GB_DB_VIRAL_CDS/
  mv GB_DB_VIRAL/*_rna_from_genomic.fna.gz GB_DB_VIRAL_RNA/
  #
  rm -f VDB.fa.gz;
  #zcat GB_DB_VIRAL/*.fna.gz | gzip -9 > VDB.fa.gz
  zcat GB_DB_VIRAL/*.fna.gz > VDB.fa.gz
  fi
#  
# ==============================================================================
# ================================ FILTER READS ================================
if [[ "$FILTER_READS" -eq "1" ]];
  then
  ./gto_fasta_extract_read_by_pattern -p "AY386330" < VDB.fa > B19.fa
  ./gto_fasta_extract_read_by_pattern -p "X04370.1" < VDB.fa > HHV3.fa
  ./gto_fasta_extract_read_by_pattern -p "DQ279927.1" < VDB.fa > HHV4.fa
  ./gto_fasta_extract_read_by_pattern -p "AF037218" < VDB.fa > HHV7.fa
  ./gto_fasta_extract_read_by_pattern -p "MH649129.1" < VDB.fa > AV.fa
  fi
#
# ==============================================================================
# =============================== RUN COMPARISON ===============================
if [[ "$RUN_COMPARISON" -eq "1" ]];
  then
  ./gto_fasta_rand_extra_chars < HHV3.fa | ./gto_fasta_to_seq > SEQ
  ./gto_reverse < SEQ > SEQ_R;
  ./GeCo -v -tm 16:200:1:5/10 -e -g 0.95 SEQ
  ./GeCo -v -tm 16:200:1:5/10 -e -g 0.95 SEQ_R
  ./gto_filter -w 51 < SEQ.iae   | awk '{ print $2;} ' > FIL.x
  ./gto_filter -w 51 < SEQ_R.iae | awk '{ print $2;} ' > FIL_R.x
  tac FIL_R.x > FIL_N.x
  ./goose-min FIL.x FIL_N.x > Merged.x
gnuplot << EOF
    reset
    set terminal pdfcairo enhanced color font 'Verdana,12'
    set output "profile1.pdf"
    set style line 101 lc rgb '#000000' lt 1 lw 4
    set border 3 front ls 101
    set tics nomirror out scale 0.75
    set format '%g'
    set size ratio 0.1
    unset key
    set yrange [:2] 
    set xrange [:]
    set xtics auto
    set ytics 0.5
    set grid 
    set ylabel "Bps"
    set xlabel "Length"
    set border linewidth 1.5
    set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 5 ps 0.4 # --- blue
    set style line 2 lc rgb '#0060ad' lt 1 lw 4 pt 6 ps 0.4 # --- green
    plot "Merged.x" using 1 with lines ls 1
EOF

  fi
#
# ==============================================================================
# ==============================================================================
