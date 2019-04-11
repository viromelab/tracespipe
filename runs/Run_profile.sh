#!/bin/bash
#
# ==============================================================================
# |                                                                            |
# |     THIS PROGRAM COMPUTES COMPLEXITY AND RELATIVE COMPLEXITY PROFILES      |                                                            # |     =================================================================      |
# |                                                                            |
# |        FILES NEEDED TO THE WHOLE COMPUTATION:                              |
# |                                                                            |
# |        ids_of_genomes.txt > FILE WITH THE GENOMES ID                       |
# |        VDB.fa             > REFERENCE DATABASE WITH THE GENOMES            |
# |                             WHEN BUILD_DATABASE=1 IT WILL BE CREATED       |
# |        reads.fq           > READS FROM METAGENOMICS [RELATIVE MODE]        |
# |                                                                            |
# ==============================================================================
#
# ==============================================================================
# ================================ DEFINITIONS =================================
GET_GTO=1;
GET_GECO=1;
#
BUILD_DATABASE=0; # USE ONLY IF THERE IS NO VDB.fa (FASTA REFERENCE DATABASE)
#
FILTER_READS=1;
RUN_COMPARISON=1;
#
RUN_RELATIVE=1; # reads.fq ARE NEEDED FOR THIS COMPUTATION
#
# ==============================================================================
# ================================== GET GTO ===================================
if [[ "$GET_GTO" -eq "1" ]];
  then
  rm -fr gto/ gto_fasta_extract_read_by_pattern
  git clone https://github.com/bioinformatics-ua/gto.git
  cd gto/src/
  make
  cp ../bin/* ../../
  cd ../../
  fi
#
# ==============================================================================
# ================================= GET GOOSE ==================================
if [[ "$GET_GECO" -eq "1" ]];
  then
  rm -fr geco/
  git clone https://github.com/pratas/geco.git
  cd geco/src/
  cmake .
  make
  cp GeCo ../../
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
  IFS=$'\r\n' GLOBIGNORE='*' command eval 'IDS_ARRAY=($(cat ids_of_genomes.txt))'
  for ids_gen in ${IDS_ARRAY[@]};
    do
    ./gto_fasta_extract_read_by_pattern -p "$ids_gen" < VDB.fa > S$ids_gen.fa
    done
  fi
#
# ==============================================================================
# =============================== RUN COMPARISON ===============================
if [[ "$RUN_COMPARISON" -eq "1" ]];
  then
  for ids_gen in ${IDS_ARRAY[@]};
    do
    seqx="S$ids_gen.fa";
    #
    ./gto_fasta_rand_extra_chars < $seqx | ./gto_fasta_to_seq > SEQ;
    ./gto_reverse < SEQ > SEQ_R;
    #
    # GET WINDOW SIZE BY SEQUENCE SIZE
    LENGTH=`./gto_info < SEQ | grep "Number of sym" | awk '{ print $5}'`;
    WINDOW_SIZE=`echo "$LENGTH / 285" | bc`;
    echo "WINDOW SIZE: $WINDOW_SIZE";
    #
    ./GeCo -v -tm 16:200:1:5/10 -e -g 0.95 SEQ
    ./GeCo -v -tm 16:200:1:5/10 -e -g 0.95 SEQ_R
    #
    ./gto_upper_bound -u 2 < SEQ.iae   > SEQ_UB
    ./gto_upper_bound -u 2 < SEQ_R.iae > SEQ_R_UB
    #
    ./gto_filter -d 2 -w $WINDOW_SIZE -c < SEQ_UB   > FIL_UB.x
    ./gto_filter -d 2 -w $WINDOW_SIZE -c < SEQ_R_UB > FIL_UB_R.x
    #
    tac FIL_UB_R.x > FIL_UB_N.x
    awk '{print $1;}' FIL_UB.x   > IDXES
    awk '{print $2;}' FIL_UB.x   > A_D
    awk '{print $2;}' FIL_UB_R.x > A_R
    #
    ./gto_min -f A_D -s A_R > A_min
    #
    paste -d '\t' IDXES A_min > PROFILE_N
    #
    gnuplot << EOF
      reset
      set terminal pdfcairo enhanced color font 'Verdana,12'
      set output "profile_self-$ids_gen.pdf"
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
      plot "PROFILE_N" using 1:2 with lines ls 1
EOF
    if [[ "$RUN_RELATIVE" -eq "1" ]];  # HERE reads.fq WILL BE NEEDED!
      then
      #
      ./gto_fasta_rand_extra_chars < $seqx | ./gto_fasta_to_seq > SEQ;
      ./gto_reverse < SEQ > SEQ_R;
      #
      ./GeCo -v -rm 16:200:1:5/10 -e -g 0.95 -c 20 -r reads.fq SEQ
      ./GeCo -v -rm 16:200:1:5/10 -e -g 0.95 -c 20 -r reads.fq SEQ_R
      #
      ./gto_upper_bound -u 2 < SEQ.iae   > SEQ_UB
      ./gto_upper_bound -u 2 < SEQ_R.iae > SEQ_R_UB
      #
      ./gto_filter -d 2 -w $WINDOW_SIZE -c < SEQ_UB   > FIL_UB.x
      ./gto_filter -d 2 -w $WINDOW_SIZE -c < SEQ_R_UB > FIL_UB_R.x
      #
      tac FIL_UB_R.x > FIL_UB_N.x
      awk '{print $1;}' FIL_UB.x   > IDXES
      awk '{print $2;}' FIL_UB.x   > A_D
      awk '{print $2;}' FIL_UB_R.x > A_R
      #
      ./gto_min -f A_D -s A_R > A_min
      #
      paste -d '\t' IDXES A_min > PROFILE_N_REL
      #
      gnuplot << EOF
        reset
        set terminal pdfcairo enhanced color font 'Verdana,12'
        set output "profile_relative-$ids_gen.pdf"
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
        set style line 2 lc rgb '#228B22' lt 1 lw 2 pt 5 ps 0.4 # --- green
        plot "PROFILE_N_REL" using 1:2 with lines ls 2
EOF
      # JOINT PLOT
      paste -d '\t' PROFILE_N A_min > PROFILES_N_AND_N_REL
      #
      gnuplot << EOF
        reset
        set terminal pdfcairo enhanced color font 'Verdana,12'
        set output "profiles-$ids_gen.pdf"
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
        set style line 2 lc rgb '#228B22' lt 1 lw 2 pt 5 ps 0.4 # --- green
        plot "PROFILES_N_AND_N_REL" using 1:2 with lines ls 1, "PROFILES_N_AND_N_REL" using 1:3 with lines ls 2 
EOF
      fi
    done
    #
  fi
#
# ==============================================================================
# ==============================================================================
