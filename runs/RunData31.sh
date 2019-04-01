#!/bin/bash
#
# ##############################################################################
# # ========================================================================== #
# # =                                                                        = #
# # =     AUTOMATIC DETECTION OF SAMPLE COMPOSITION IN THE READS OF #31      = #
# # =                                                                        = #
# # ========================================================================== #
# ##############################################################################
#
# 1 -> TURNS THE RESPECTIVE COMPUTATION ON
# 0 -> IGNORES THE RESPECTIVE ACTION
#
GET_UTILS=0;
GET_TRIMMOMATIC=0;
GET_FALCON=0;
GET_MAGNET=0;
GET_GULL=0;
GET_GTO=0;
#
BUILD_DATABASE=0;
GET_PHIX=0;
GET_READS=0;
#
RUN_TRIM_METAGENOMICS=1;
RUN_GULL=1;
#
EMAIL_TO="diogo.pratas@helsinki.fi";
EMAIL_FROM="kone@helsinki.fi";
REP_NAME="I31_REPORT.tex";
WRITE_LATEX_HEADER=1;
WRITE_LAC=1;
WRITE_LATEX_TAIL=1;
COMPILE_PDF=1;
#
# ==============================================================================
# ##############################################################################
# ==============================================================================
#
#
# ==============================================================================
# ================================= GET UTILS ==================================
if [[ "$GET_UTILS" -eq "1" ]];
  then
  sudo apt-get install librsvg2-bin cmake git python3-pip
  pip3 install conda
  fi
#
# ==============================================================================
# ================================= GET FALCON =================================
if [[ "$GET_FALCON" -eq "1" ]];
  then
  rm -fr falcon FALCON FALCON-FILTER FALCON-EYE
  git clone https://github.com/pratas/falcon.git
  cd falcon/src/
  cmake .
  make
  cp FALCON ../../
  cp FALCON-FILTER ../../
  cp FALCON-EYE ../../
  cd ../../
  fi
#
# ==============================================================================
# ================================= GET MAGNET =================================
if [[ "$GET_MAGNET" -eq "1" ]];
  then
  rm -fr magnet Magnet
  git clone https://github.com/pratas/Magnet.git
  cd Magnet/src
  cmake .
  make
  cp Magnet ../../MAGNET
  cd ../../
  fi
#
# ==============================================================================
# ================================== GET GULL ==================================
if [[ "$GET_GULL" -eq "1" ]]; 
  then
  rm -fr GULL GULL-map GULL-visual
  git clone https://github.com/pratas/GULL.git
  cd GULL/src/
  cmake .
  make
  cp GULL-map ../../
  cp GULL-visual ../../
  cd ../../
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
  cd ../../
  fi
#
# ==============================================================================
# ============================== GET TRIMMOMATIC ===============================
if [[ "$GET_TRIMMOMATIC" -eq "1" ]];
  then
  conda install -c bioconda trimmomatic
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
  zcat GB_DB_VIRAL/*.fna.gz | gzip -9 > VDB.fa.gz
  #
  # CHECK FOR ">" SYMBOLS INSIDE THE HEADER:
  # zcat VDB.fa.gz | grep ">" | cut -c2- | grep ">"
  fi
#  
# ==============================================================================
# ================================ GET PHIX SEQ ================================
if [[ "$GET_PHIX" -eq "1" ]];
  then
  zcat VDB.fa.gz | ./gto_fasta_extract_read_by_pattern -p "phiX174" > F_PHIX.fa
  fi
#
# ==============================================================================
# ================================= GET READS ==================================
if [[ "$GET_READS" -eq "1" ]]; 
  then
  #
  # FOR COMPLETE AUTOMATIC LOGIN, USE EXPECT SHELL [[ IT WILL ADD PASSWORD ]]
  scp -r dpratas@ssh.fimm.fi:/fas/NGS/pipes/fastq/fimm_ngs_mtoppinen/autopsy_NGS_pilot_31/*_001.fastq.gz .
  #
  fi
#
# ==============================================================================
# ============================= WRITE LATEX HEADER =============================
if [[ "$WRITE_LATEX_HEADER" -eq "1" ]];
  then
  echo "\\documentclass{article}" > $REP_NAME
  echo "\\usepackage{graphicx}" >> $REP_NAME
  echo "\\begin{document}" >> $REP_NAME
  echo "%" >> $REP_NAME
  echo "\\title{SAMPLE COMPOSITION OF 31}" >> $REP_NAME
  echo "\\author{T. Ietokone}" >> $REP_NAME
  echo "%" >> $REP_NAME
  echo "\\maketitle" >> $REP_NAME
  echo "%" >> $REP_NAME
  echo "\\begin{abstract} Report on the next-generation sequencing reads generated by a HiSeq V2 Rapid PE101 run using different organs from an individual named 31. The reads were enriched according to several viral genomes from GeneBank. The reads are trimmed with TRIMMOMATIC. The PhiX (sequencing control) is removed from the reads using MAGNET. The database includes the latest complete genomes from the NCBI virual database. The metagenomic analysis is perfom with FALCON.\\end{abstract}" >>  $REP_NAME 
  echo "%" >> $REP_NAME
  fi
#
# ==============================================================================
# ========================== RUN TRIM & METAGENOMICS ===========================
if [[ "$RUN_TRIM_METAGENOMICS" -eq "1" ]];
  then
  #
  gunzip VDB.fa.gz
  declare -a READS=("healthy-skin:MT18021B_S1_L001_R1_001.fastq.gz:MT18021B_S1_L001_R2_001.fastq.gz:MT18021B_S1_L002_R1_001.fastq.gz:MT18021B_S1_L002_R2_001.fastq.gz" "skin:MT18021C_S2_L001_R1_001.fastq.gz:MT18021C_S2_L001_R2_001.fastq.gz:MT18021C_S2_L002_R1_001.fastq.gz:MT18021C_S2_L002_R2_001.fastq.gz" "bone:MT18021D_S3_L001_R1_001.fastq.gz:MT18021D_S3_L001_R2_001.fastq.gz:MT18021D_S3_L002_R1_001.fastq.gz:MT18021D_S3_L002_R2_001.fastq.gz" "colon:MT18021E_S4_L001_R1_001.fastq.gz:MT18021E_S4_L001_R2_001.fastq.gz:MT18021E_S4_L002_R1_001.fastq.gz:MT18021E_S4_L002_R2_001.fastq.gz" "heart:MT18021F_S5_L001_R1_001.fastq.gz:MT18021F_S5_L001_R2_001.fastq.gz:MT18021F_S5_L002_R1_001.fastq.gz:MT18021F_S5_L002_R2_001.fastq.gz" "liver:MT18021G_S6_L001_R1_001.fastq.gz:MT18021G_S6_L001_R2_001.fastq.gz:MT18021G_S6_L002_R1_001.fastq.gz:MT18021G_S6_L002_R2_001.fastq.gz" "spleen:MT18022B_S7_L001_R1_001.fastq.gz:MT18022B_S7_L001_R2_001.fastq.gz:MT18022B_S7_L002_R1_001.fastq.gz:MT18022B_S7_L002_R2_001.fastq.gz" "kidney:MT18022C_S8_L001_R1_001.fastq.gz:MT18022C_S8_L001_R2_001.fastq.gz:MT18022C_S8_L002_R1_001.fastq.gz:MT18022C_S8_L002_R2_001.fastq.gz" "lung:MT18022D_S9_L001_R1_001.fastq.gz:MT18022D_S9_L001_R2_001.fastq.gz:MT18022D_S9_L002_R1_001.fastq.gz:MT18022D_S9_L002_R2_001.fastq.gz" "plasma:MT18022E_S10_L001_R1_001.fastq.gz:MT18022E_S10_L001_R2_001.fastq.gz:MT18022E_S10_L002_R1_001.fastq.gz:MT18022E_S10_L002_R2_001.fastq.gz" "blood:MT18022F_S11_L001_R1_001.fastq.gz:MT18022F_S11_L001_R2_001.fastq.gz:MT18022F_S11_L002_R1_001.fastq.gz:MT18022F_S11_L002_R2_001.fastq.gz" "bone-marrow:MT18022G_S12_L001_R1_001.fastq.gz:MT18022G_S12_L001_R2_001.fastq.gz:MT18022G_S12_L002_R1_001.fastq.gz:MT18022G_S12_L002_R2_001.fastq.gz" "teeth:MT18023B_S13_L001_R1_001.fastq.gz:MT18023B_S13_L001_R2_001.fastq.gz:MT18023B_S13_L002_R1_001.fastq.gz:MT18023B_S13_L002_R2_001.fastq.gz" "brain:MT18023C_S14_L001_R1_001.fastq.gz:MT18023C_S14_L001_R2_001.fastq.gz:MT18023C_S14_L002_R1_001.fastq.gz:MT18023C_S14_L002_R2_001.fastq.gz")

  ## LOOP TRHOUGH THE ARRAY READS AND SPLITS THE DATA INTO 5 FIELDS
  for read in "${READS[@]}"
    do
    ORGAN_T=`echo $read | tr ':' '\t' | awk '{ print $1 }'`;
    SPL_R1A=`echo $read | tr ':' '\t' | awk '{ print $2 }'`;
    SPL_R2A=`echo $read | tr ':' '\t' | awk '{ print $3 }'`;
    SPL_R1B=`echo $read | tr ':' '\t' | awk '{ print $4 }'`;
    SPL_R2B=`echo $read | tr ':' '\t' | awk '{ print $5 }'`;
    echo "RUN ORGAN=$ORGAN_T F1=$SPL_R1A R1=$SPL_R2A F2=$SPL_R1B R2=$SPL_R2B";

    ## MERGE FILES
    echo "MERGGING FILES ...";
    rm -f FW_READS.fq.gz RV_READS.fq.gz
    zcat $SPL_R1A $SPL_R1B | gzip > FW_READS.fq.gz
    zcat $SPL_R2A $SPL_R2B | gzip > RV_READS.fq.gz

    ## LOCATE TRIMMOMATIC PRIMERS:
    ADAPTERS_PATH=`locate TruSeq3-PE-2.fa | grep conda | head -n 1`;

    ## RUN TRIMMING AND FILTERING
    echo "RUNNING TRIMMOMATIC ...";
    rm -f o_fw_pr.fq.gz o_fw_unpr.fq.gz o_rv_pr.fq.gz o_rv_unpr.fq.gz;
    trimmomatic PE -threads 8 -phred33 FW_READS.fq.gz RV_READS.fq.gz o_fw_pr.fq.gz o_fw_unpr.fq.gz o_rv_pr.fq.gz o_rv_unpr.fq.gz ILLUMINACLIP:$ADAPTERS_PATH:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35 CROP:96
    #XXX: YOU MAY REMOVE CROP:96 ALTHOUGH SOME ACCEPTABLE(?) NOISE WAS FOUND THERE

    ## CONCATENATE ALL THE FILES
    echo "CONCATENATING FILES ...";
    zcat o_fw_pr.fq.gz o_fw_unpr.fq.gz o_rv_pr.fq.gz o_rv_unpr.fq.gz > $ORGAN_T.fq

    ## REMOVE PHIX FROM THE SAMPLES
    echo "RUNNING MAGNET ...";
    ./MAGNET -v -F -n 8 -m 6:1:1:0/0 -m 13:20:1:3/10 -t 0.9 -i -o $ORGAN_T-NO_PHIX.fq F_PHIX.fa $ORGAN_T.fq

    ## RUN METAGENOMIC COMPOSITION
    echo "RUNNING FALCON ...";
    ./FALCON -v -n 8 -t 38 -F -Z -l 47 -c 20 -x top-$ORGAN_T.csv -y $ORGAN_T.com $ORGAN_T-NO_PHIX.fq VDB.fa
    ./FALCON-FILTER -v -F -t 1.0 -o $ORGAN_T.pos $ORGAN_T.com
    ./FALCON-EYE -v -e 1 -F -o $ORGAN_T.svg $ORGAN_T.pos

    ## CONVERT SVG OUTPUT TO PDF
    rsvg-convert -f pdf -o $ORGAN_T.pdf $ORGAN_T.svg

    ## RUN GULL FOR INTER-GENOMIC SIMILARITY ANALYSIS
    if [[ "$RUN_GULL" -eq "1" ]]; 
      then
      cat top-$ORGAN_T.csv | awk '{ if($3 > 0.01) print $1"\t"$2"\t"$3"\t"$4; }' \
      | awk '{ print $4;}' | tr '|' '\t' | tr '_' '\t' | awk '{ print $1;}' > GIS-$ORGAN_T;
      idx=0;
      cat GIS-$ORGAN_T | while read line
        do
        namex=`echo $line | tr ' ' '_'`;
        if [[ "$idx" -eq "0" ]]; 
	  then
          printf "%s" "$namex" > $ORGAN_T-FNAMES.fil;
          else
          printf ":%s" "$namex" >> $ORGAN_T-FNAMES.fil;
          fi
        ./gto_fasta_extract_read_by_pattern -p "$line" < VDB.fa > $namex;
        ((idx++));
        done
      ./GULL-map -v -m 6:1:1:0/0 -m 13:20:1:3/10 -m 20:100:1:5/10 -c 30 -n 8 -x $ORGAN_T-MATRIX.csv `cat $ORGAN_T-FNAMES.fil`
      ./GULL-visual -v -w 25 -a 8 -x $ORGAN_T-HEAT.svg $ORGAN_T-MATRIX.csv
    fi

    ## CONVERT SVG OUTPUT TO PDF
    if [[ "$RUN_GULL" -eq "1" ]];
      then
      rsvg-convert -f pdf -o $ORGAN_T-HEAT.pdf $ORGAN_T-HEAT.svg
      fi

    ## WRITE TO LATEX FILE
    if [[ "$WRITE_LAC" -eq "1" ]];
      then
      echo "\\section{$ORGAN_T}" >> $REP_NAME
      echo "%" >> $REP_NAME
      #
      echo "\\begin{figure}[h]" >> $REP_NAME
      echo "\\caption{Metagenomics analysis map for I31 using $ORGAN_T.} " >> $REP_NAME
      echo "\\centering" >> $REP_NAME
      echo "\\includegraphics[width=0.5\\textwidth]{$ORGAN_T.pdf}" >> $REP_NAME
      echo "\\end{figure}" >> $REP_NAME
      #
      echo "%" >> $REP_NAME
      #
      if [[ "$RUN_GULL" -eq "1" ]];
        then
        echo "%" >> $REP_NAME
	echo "\\begin{figure}[h]" >> $REP_NAME
        echo "\\caption{Intra-genomic similarity analysis map for organisms with some level of similarity in the metagenomic analsys of I31 (using $ORGAN_T).} " >> $REP_NAME
        echo "\\centering" >> $REP_NAME
        echo "\\includegraphics[width=0.5\\textwidth]{$ORGAN_T-HEAT.pdf}" >> $REP_NAME
        echo "\\end{figure}" >> $REP_NAME
        echo "%" >> $REP_NAME
        fi
      #
      fi
    done
  gzip VDB.fa
  fi
#
# ==============================================================================
# =============================== WRITE LATEX TAIL =============================
if [[ "$WRITE_LATEX_TAIL" -eq "1" ]];
  then
  echo "%" >> $REP_NAME
  echo "\\end{document}" >> $REP_NAME
  fi
#
# ==============================================================================
# ================================= COMPILE PDF ================================
if [[ "$COMPLIE_PDF" -eq "1" ]];
  then
  pdflatex --jobname=I31.pdf $REP_NAME
  fi
#
# ==============================================================================
# ================================== EMAIL PDF =================================
if [[ "$EMAIL_PDF" -eq "1" ]];
  then
  #sudo apt-get install mailutils
  #dpkg-reconfigure -plow postfix ## Internet connection > ... > ... > ...
  echo "Results attach" | mail -s "Report #31" -A I31.pdf $EMAIL_FROM $EMAIL_TO
  fi
#
# ==============================================================================
# ==============================================================================
