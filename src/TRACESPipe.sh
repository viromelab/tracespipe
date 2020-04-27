#!/bin/bash
#
##################################################################################
# ============================================================================== #
# =                                                                            = #
# =                              T R A C E S P i p e                           = #
# =                                                                            = #
# =     A Next-Generation Sequencing pipeline for identification, assembly,    = #
# =     and analysis of viral and human-host genomes at multi-organ level.     = #
# =                                                                            = #
# ============================================================================== #
##################################################################################
#
SHOW_HELP=0;
SHOW_VERSION=0;
FORCE=0;
FLUSH_LOGS=0;
GET_THREADS=0;
THREADS=0;
#
INSTALL=0;
UPDATE=0;
BUILD_VDB_ALL=0;
BUILD_VDB_REF=0;
BUILD_UDB=0;
#
GEN_ADAPTERS=0;
GET_PHIX=0;
GET_MITO=0;
GET_CY=0;
GET_EXTRA=0;
ADD_EXTRA_SPECIFIC=0;
ADD_FASTA=0;
NEW_FASTA="";
#
CREATE_BLAST_DB=0;
UPDATE_BLAST_DB=0;
SEARCH_BLAST_DB=0;
BLAST_QUERY="";
#
RUN_ANALYSIS=0;
#
RUN_META_ON=0;
RUN_PROFILES_ON=0;
RUN_META_NON_VIRAL_ON=0;
#
RUN_MITO_ON=0;
RUN_MITO_DAMAGE_ON=0;
#
VIEW_TOP=0;
#
REMOVE_DUPLICATIONS=0;
#
RUN_B19_ON=0;
RUN_HV1_ON=0;
RUN_HV2_ON=0;
RUN_HV3_ON=0;
RUN_HV4_ON=0;
RUN_HV5_ON=0;
RUN_HV6_ON=0;
RUN_HV6A_ON=0;
RUN_HV6B_ON=0;
RUN_HV7_ON=0;
RUN_HV8_ON=0;
#
RUN_POLY1_ON=0;
RUN_POLY2_ON=0;
RUN_POLY3_ON=0;
RUN_POLY4_ON=0;
RUN_POLY5_ON=0;
RUN_POLY6_ON=0;
RUN_POLY7_ON=0;
RUN_POLY8_ON=0;
RUN_POLY9_ON=0;
RUN_POLY10_ON=0;
RUN_POLY11_ON=0;
RUN_POLY12_ON=0;
RUN_POLY13_ON=0;
RUN_POLY14_ON=0;
#
RUN_TTV_ON=0;
RUN_HBOV1_ON=0;
RUN_HBOVNOT1_ON=0;
RUN_HBV_ON=0;
RUN_HPV_ON=0;
RUN_VARV_ON=0;
RUN_SV40_ON=0;
RUN_CUTA_ON=0;
RUN_HERV_ON=0;
#
RUN_VISUAL_ALIGN=0;
#
RUN_COVERAGE_TABLE=0;
RUN_COVERAGE_TABLE_CSV=0;
RUN_COVERAGE_PROFILE=0;
COVERAGE_NAME="";
MAX_COVERAGE_PROFILE=0;
#
RUN_CHANGE_MT=0;
#
RUN_DECRYPT=0;
RUN_ENCRYPT=0;
#
RUN_SPECIFIC=0;
RUN_SPECIFIC_SENSITIVE=0;
#
RUN_CY_ON=0;
RUN_CY_QUANT_ON=0;
#
RUN_DE_NOVO_ASSEMBLY=0;
#
RUN_HYBRID=0;
#
# ==============================================================================
#
# DEFAULT VALUES:
#
MINIMAL_SIMILARITY_VALUE=1.0;
TSIZE=10;
#
# ==============================================================================
#
declare -a VIRUSES=("B19" "HV1" "HV2" "HV3" "HV4" "HV5" "HV6" "HV6A" "HV6B" 
                    "HV7" "HV8" "POLY1" "POLY2" "POLY2" "POLY3" "POLY4" "POLY5" 
		    "POLY6" "POLY7" "POLY8" "POLY9" "POLY10" "POLY11" "POLY12" 
		    "POLY13" "POLY14" "HBV" "HPV" "TTV" "HBOV1" "HBOVNOT1" 
		    "VARV" "SV40" "CUTA" "HERV");
#
# ==============================================================================
# CHECK INTERNAL SYSTEM FILES
#
CHECK_FILTERING_SYSTEM_FILES () {
  for VIRUS in "${VIRUSES[@]}"
    do
    if [ ! -f TRACES_get_best_$VIRUS.sh ];
      then
      echo -e "\e[31mERROR: TRACES_get_best_$VIRUS.sh file not found!\e[0m"
      echo "This file may have been deleted acidentally or VIRAL array changed."
      echo "System is corrupted!"
      exit 1;
    fi
    done
  }
#
# ==============================================================================
# CHECK IF FILES EXIST
#
CHECK_META_INFO () {
  if [ ! -f ../meta_data/meta_info.txt ];
    then
    echo -e "\e[31mERROR: meta_info.txt file not found!\e[0m"
    echo "Please create a meta information file before the run."
    echo "For addition information, see the instructions at the web page."
    exit 1;
    fi
  }
#
#
CHECK_GZIP_FILES () {
  FILE_GZ1=`file ../input_data/$1 | grep gzip | wc -l`;
  FILE_GZ2=`file ../input_data/$2 | grep gzip | wc -l`;
  if [[ "$FILE_GZ1" != "1" ]] || [[ "$FILE_GZ2" != "1" ]];
    then
    echo -e "\e[31mERROR: the input reads are not gzipped!\e[0m"
    echo "Please compress the reads (input_data) with gzip."
    echo "For addition information, see the instructions at the web page."
    exit 1;
    fi
  }
#
#
CHECK_VDB () {
  if [ ! -f VDB.fa ];
    then
    echo -e "\e[31mERROR: viral database (VDB.fa) not found!\e[0m"
    echo "TIP: before this, run: ./TRACESPipe.sh --build-viral"
    echo "For addition information, see the instructions at the web page."
    exit 1;
    fi
  }
#
#
CHECK_DB () {
  if [ ! -f DB.fa ];
    then
    echo -e "\e[31mERROR: Non-viral database FASTA file (DB.fa) not found!\e[0m"
    echo "TIP: first run ./TRACESPipe --build-unviral"
    exit 1;
    fi
  }
#
#
CHECK_PHIX () {
  if [ ! -f F_PHIX.fa ];
    then
    echo -e "\e[31mERROR: viral PhiX (F_PHIX.fa) not found!\e[0m"
    echo "TIP: before this, run: ./TRACESPipe.sh --get-phix"
    echo "For addition information, see the instructions at the web page."
    exit 1;
    fi
  }
#
#
CHECK_ADAPTERS () {
  if [ ! -f adapters.fa ];
    then
    echo -e "\e[31mERROR: adapter sequences (adapters.fa) not found!\e[0m"
    echo "TIP: before this, run: ./TRACESPipe.sh --gen-adapters"
    echo "For addition information, see the instructions at the web page."
    exit 1;
    fi
  }
#
#
CHECK_ADAPTERS_AR () {
  if [ ! -f adapters_ar.fa ];
    then
    echo -e "\e[31mERROR: adapter sequences (adapters_ar.fa) not found!\e[0m"
    echo "TIP: before this, include valid adapters for AdapterRemoval"
    echo "For addition information, see the instructions at the web page."
    exit 1;
    fi
  }
#
#
CHECK_MT_DNA () {
  if [ ! -f mtDNA.fa ];
    then
    echo -e "\e[31mERROR: reference mitochondrial DNA (mtDNA.fa) not found!\e[0m"
    echo "TIP: before this, run: ./TRACESPipe.sh --get-mito"
    echo "For addition information, see the instructions at the web page."
    exit 1;
    fi
  }
#
#
CHECK_CY_DNA () {
  if [ ! -f cy.fa ];
    then
    echo -e "\e[31mERROR: reference y-chromosome DNA (cy.fa) not found!\e[0m"
    echo "TIP: before this, run: ./TRACESPipe.sh --get-y-chromo"
    echo "For addition information, see the instructions at the web page."
    exit 1;
    fi
  }
#
#
CHECK_TOP () {
  if [ ! -f top-$1.csv ];
    then
    echo -e "\e[31mERROR: top-$1.csv not found!\e[0m"
    echo "Viral alignments are only possible after metagenomic analysis".
    echo "(Unless is a specific viral alignment by ID/PATTERN)."
    echo "TIP: before this, run: ./TRACESPipe.sh --run-meta"
    echo "For addition information, see the instructions at the web page."
    exit 1;
    fi
  }
#
#
CHECK_ENRICH () {
  if [ ! -f ../system_files/ids_enrichment.txt ];
    then
    echo -e "\e[31mERROR: ids_enrichment.txt file not found!\e[0m"
    echo "Please create a file ids_enrichment.txt with the GIS (each one line by line)"
    echo "and add the file at ../system_files/ folder before the run."
    echo "For addition information, see the instructions at the web page."
    exit 1;
    fi
  }
#
#
CHECK_E_FILE () {
  if [ ! -f $1 ];
    then
    echo -e "\e[31mERROR: $1 file not found!\e[0m"
    echo "For addition information, see the instructions at the web page."
    exit 1;
    fi
  }
#
# ==============================================================================
#
ALIGN_AND_CONSENSUS () {
  #
  V_TAG="$1";
  echo -e "\e[34m[TRACESPipe]\e[32m Assessing $V_TAG best reference ...\e[0m";
  CHECK_TOP "$ORGAN_T";
  V_INFO=`./TRACES_get_best_$V_TAG.sh $ORGAN_T`;
  V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
  V_VAL=`echo "$V_INFO" | awk '{ print $1; }'`;
  echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
  if [[ "$V_GID" != "-" && "$V_VAL" > "$2" ]];
    then
    echo -e "\e[34m[TRACESPipe]\e[96m Similarity best match: $V_INFO\e[0m";
    echo -e "\e[34m[TRACESPipe]\e[32m Extracting sequence from VDB.fa\e[0m";
    CHECK_VDB;
    gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-$V_TAG.fa 2>> ../logs/Log-$ORGAN_T.txt;
    echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
    echo -e "\e[34m[TRACESPipe]\e[32m Aligning reads to $V_TAG best reference with bowtie2 ...\e[0m";
    ./TRACES_viral_align_reads.sh $ORGAN_T-$V_TAG.fa $ORGAN_T $V_TAG $THREADS $REMOVE_DUPLICATIONS 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
    echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
    #
    echo -e "\e[34m[TRACESPipe]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
    ./TRACES_viral_consensus.sh $ORGAN_T-$V_TAG.fa viral_aligned_sorted-$ORGAN_T-$V_TAG.bam $ORGAN_T $V_TAG 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
    #
    mkdir -p ../output_data/TRACES_viral_alignments;
    cp $ORGAN_T-$V_TAG.fa ../output_data/TRACES_viral_alignments/
    cp $ORGAN_T-$V_TAG.fa.fai ../output_data/TRACES_viral_alignments/
    mv viral_aligned_sorted-$ORGAN_T-$V_TAG.bam ../output_data/TRACES_viral_alignments/
    mv viral_aligned_sorted-$ORGAN_T-$V_TAG.bam.bai ../output_data/TRACES_viral_alignments/
    mkdir -p ../output_data/TRACES_viral_consensus;
    #rm -f ../output_data/TRACES_viral_consensus/*
    mv $V_TAG-consensus-$ORGAN_T.fa ../output_data/TRACES_viral_consensus/
    mkdir -p ../output_data/TRACES_viral_bed;
    #rm -f ../output_data/TRACES_viral_bed/*
    mv $V_TAG-calls-$ORGAN_T.bed ../output_data/TRACES_viral_bed/
    mv $V_TAG-coverage-$ORGAN_T.bed ../output_data/TRACES_viral_bed/
    mv $V_TAG-zero-coverage-$ORGAN_T.bed ../output_data/TRACES_viral_bed/
    mkdir -p ../output_data/TRACES_viral_statistics;
    echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
    echo -e "\e[34m[TRACESPipe]\e[32m Calculating coverage ...\e[0m";
    ./TRACES_overall_virus.sh $V_TAG $ORGAN_T
    C_BREADTH=`cat ../output_data/TRACES_viral_statistics/$V_TAG-total-horizontal-coverage-$ORGAN_T.txt`;
    C_DEPTH=`cat ../output_data/TRACES_viral_statistics/$V_TAG-total-depth-coverage-$ORGAN_T.txt`;
    echo -e "\e[34m[TRACESPipe]\e[1m Breadth (H) coverage: $C_BREADTH \e[0m";
    echo -e "\e[34m[TRACESPipe]\e[1m Depth-x (V) coverage: $C_DEPTH \e[0m";
    echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
    fi
  }
#
# ==============================================================================
#
if [ "$#" -eq 0 ];
  then
  SHOW_HELP=1;
  fi
#
POSITIONAL=();
#
while [[ $# -gt 0 ]]
  do
  i="$1";
  case $i in
    -h|--help|?)
      SHOW_HELP=1;
      shift
    ;;
    -v|-V|--version)
      SHOW_VERSION=1;
      shift
    ;;
    -f|-F|--force)
      FORCE=1;
      shift
    ;;
    -flog|--flush-logs)
      FLUSH_LOGS=1;
      shift
    ;;
    -gmt|--get-max-threads)
      GET_THREADS=1;
      shift
    ;;
    -t|--threads)
      THREADS=$2;
      shift 2;
    ;;
    -i|--install)
      INSTALL=1;
      SHOW_HELP=0;
      shift
    ;;
    -up|--update)
      UPDATE=1;
      SHOW_HELP=0;
      shift
    ;;
    -vdb|--build-viral)
      BUILD_VDB_ALL=1;
      SHOW_HELP=0;
      shift
    ;;
    -vdbr|--build-viral-r)
      BUILD_VDB_REF=1;
      SHOW_HELP=0;
      shift
    ;;
    -udb|--build-unviral)
      BUILD_UDB=1;
      SHOW_HELP=0;
      shift
    ;;
    -gad|--gen-adapters)
      GEN_ADAPTERS=1;
      SHOW_HELP=0;
      shift
    ;;
    -top|--view-top)
      VIEW_TOP=1;
      VTOP_SIZE="$2";
      SHOW_HELP=0;
      shift 2
    ;;
    -gp|--get-phix)
      GET_PHIX=1;
      SHOW_HELP=0;
      shift
    ;;
    -cbn|--create-blast-db)
      CREATE_BLAST_DB=1;
      SHOW_HELP=0;
      shift
    ;;
    -ubn|--update-blast-db)
      UPDATE_BLAST_DB=1;
      SHOW_HELP=0;
      shift
    ;;
    -sfs|--search-blast-db)
      SEARCH_BLAST_DB=1;
      BLAST_QUERY="$2";
      SHOW_HELP=0;
      shift 2
    ;;
    -rdup|--remove-dup)
      REMOVE_DUPLICATIONS=1;
      SHOW_HELP=0;
      shift
    ;;    
    -gm|--get-mito)
      GET_MITO=1;
      SHOW_HELP=0;
      shift
    ;;
    -mis|--min-similarity)
      MINIMAL_SIMILARITY_VALUE="$2";
      SHOW_HELP=0;
      shift 2
    ;;
    -cmt|--change-mito)
      RUN_CHANGE_MT=1;
      NEW_MT="$2";
      SHOW_HELP=0;
      shift 2
    ;;
    -gy|--get-y-chromo)
      GET_CY=1;
      SHOW_HELP=0;
      shift
    ;;
    -gax|--get-all-aux)
      GEN_ADAPTERS=1;
      GET_PHIX=1;
      GET_MITO=1;
      GET_CY=1;
      SHOW_HELP=0;
      shift
    ;;
    -gx|--get-extra-vir)
      GET_EXTRA=1;
      SHOW_HELP=0;
      shift
    ;;
    -aes|--add-extra-seq)
      ADD_EXTRA_SEQ=1;
      NEW_SEQ_ID="$2";
      SHOW_HELP=0;
      shift 2
    ;;
    -afs|--add-fasta)
      ADD_FASTA=1;
      NEW_FASTA="$2";
      SHOW_HELP=0;
      shift 2
    ;;
    -ra|--run-analysis)
      RUN_ANALYSIS=1;
      RUN_META_ON=1;
      RUN_PROFILES_ON=1;
      RUN_META_NON_VIRAL_ON=1;
      RUN_B19_ON=1;
      RUN_HV1_ON=1;
      RUN_HV2_ON=1;
      RUN_HV3_ON=1;
      RUN_HV4_ON=1;
      RUN_HV5_ON=1;
      RUN_HV6_ON=1;
      RUN_HV6A_ON=1;
      RUN_HV6B_ON=1;
      RUN_HV7_ON=1;
      RUN_HV8_ON=1;
      RUN_POLY1_ON=1;
      RUN_POLY2_ON=1;
      RUN_POLY3_ON=1;
      RUN_POLY4_ON=1;
      RUN_POLY5_ON=1;
      RUN_POLY6_ON=1;
      RUN_POLY7_ON=1;
      RUN_POLY8_ON=1;
      RUN_POLY9_ON=1;
      RUN_POLY10_ON=1;
      RUN_POLY11_ON=1;
      RUN_POLY12_ON=1;
      RUN_POLY13_ON=1;
      RUN_POLY14_ON=1;
      RUN_TTV_ON_ON=1;
      RUN_HBOV1_ON=1;
      RUN_HBOVNOT1_ON=1;
      RUN_HBV_ON=1;
      RUN_HPV_ON=1;
      RUN_VARV_ON=1;
      RUN_SV40_ON=1;
      RUN_CUTA_ON=1;
      RUN_HERV_ON=1;
      RUN_CY_ON=1;
      RUN_CY_QUANT_ON=1;
      RUN_MITO_ON=1;
      RUN_DE_NOVO_ASSEMBLY=1;
      RUN_HYBRID=1;
      SHOW_HELP=0;
      shift
    ;;
    -iss|--inter-sim-size)
      TSIZE="$2";
      SHOW_HELP=0;
      shift 2;
    ;;
    -rsr|--run-specific)
      RUN_ANALYSIS=1;
      RUN_SPECIFIC=1;
      SPECIFIC_ID="$2";
      SHOW_HELP=0;
      shift 2;
    ;;
    -cmax|--max-coverage)
      MAX_COVERAGE_PROFILE=1;
      COVERAGE_MAX="$2";
      SHOW_HELP=0;
      shift 2;
    ;;
    -covp|--coverage-profile)
      RUN_COVERAGE_PROFILE=1;
      COVERAGE_NAME="$2";
      SHOW_HELP=0;
      shift 2;
    ;;
    -rsx|--run-extreme)
      RUN_ANALYSIS=1;
      RUN_SPECIFIC_SENSITIVE=1;
      SPECIFIC_ID="$2";
      SHOW_HELP=0;
      shift 2;
    ;;
    -rm|--run-meta)
      RUN_ANALYSIS=1;
      RUN_META_ON=1;
      SHOW_HELP=0;
      shift
    ;;      
    -ro|--run-meta-nv)
      RUN_ANALYSIS=1;
      RUN_META_NON_VIRAL_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rmt|--run-mito)
      RUN_ANALYSIS=1;
      RUN_MITO_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rmtd|--run-mito-dam)
      RUN_ANALYSIS=1;
      RUN_MITO_DAMAGE_ON=1;
      SHOW_HELP=0;
      shift
    ;;      
    -rava|--run-all-v-alig)
      RUN_ANALYSIS=1;
      RUN_B19_ON=1;
      RUN_HV1_ON=1;
      RUN_HV2_ON=1;
      RUN_HV3_ON=1;
      RUN_HV4_ON=1;
      RUN_HV5_ON=1;
      RUN_HV6_ON=1;
      RUN_HV6A_ON=1;
      RUN_HV6B_ON=1;
      RUN_HV7_ON=1;
      RUN_HV8_ON=1;
      RUN_POLY1_ON=1;
      RUN_POLY2_ON=1;
      RUN_POLY3_ON=1;
      RUN_POLY4_ON=1;
      RUN_POLY5_ON=1;
      RUN_POLY6_ON=1;
      RUN_POLY7_ON=1;
      RUN_POLY8_ON=1;
      RUN_POLY9_ON=1;
      RUN_POLY10_ON=1;
      RUN_POLY11_ON=1;
      RUN_POLY12_ON=1;
      RUN_POLY13_ON=1;
      RUN_POLY14_ON=1;
      RUN_TTV_ON=1;
      RUN_HBOV1_ON=1;
      RUN_HBOVNOT1_ON=1;
      RUN_HBV_ON=1;
      RUN_HPV_ON=1;
      RUN_VARV_ON=1;
      RUN_SV40_ON=1;
      RUN_CUTA_ON=1;
      RUN_HERV_ON=1;
      shift
    ;;
    -rb19|--run-b19)
      RUN_ANALYSIS=1;
      RUN_B19_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rh1|--run-hv1)
      RUN_ANALYSIS=1;
      RUN_HV1_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rh2|--run-hv2)
      RUN_ANALYSIS=1;
      RUN_HV2_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rh3|--run-hv3)
      RUN_ANALYSIS=1;
      RUN_HV3_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rh4|--run-hv4)
      RUN_ANALYSIS=1;
      RUN_HV4_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rh5|--run-hv5)
      RUN_ANALYSIS=1;
      RUN_HV5_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rh6|--run-hv6)
      RUN_ANALYSIS=1;
      RUN_HV6_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rh6a|--run-hv6a)
      RUN_ANALYSIS=1;
      RUN_HV6A_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rh6b|--run-hv6b)
      RUN_ANALYSIS=1;
      RUN_HV6B_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rh7|--run-hv7)
      RUN_ANALYSIS=1;
      RUN_HV7_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rh8|--run-hv8)
      RUN_ANALYSIS=1;
      RUN_HV8_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rp1|--run-poly1)
      RUN_ANALYSIS=1;
      RUN_POLY1_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rp2|--run-poly2)
      RUN_ANALYSIS=1;
      RUN_POLY2_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rp3|--run-poly3)
      RUN_ANALYSIS=1;
      RUN_POLY3_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rp4|--run-poly4)
      RUN_ANALYSIS=1;
      RUN_POLY4_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rp5|--run-poly5)
      RUN_ANALYSIS=1;
      RUN_POLY5_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rp6|--run-poly6)
      RUN_ANALYSIS=1;
      RUN_POLY6_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rp7|--run-poly7)
      RUN_ANALYSIS=1;
      RUN_POLY7_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rp8|--run-poly8)
      RUN_ANALYSIS=1;
      RUN_POLY8_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rp9|--run-poly9)
      RUN_ANALYSIS=1;
      RUN_POLY9_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rp10|--run-poly10)
      RUN_ANALYSIS=1;
      RUN_POLY10_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rp11|--run-poly11)
      RUN_ANALYSIS=1;
      RUN_POLY11_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rp12|--run-poly12)
      RUN_ANALYSIS=1;
      RUN_POLY12_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rp13|--run-poly13)
      RUN_ANALYSIS=1;
      RUN_POLY13_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rp14|--run-poly14)
      RUN_ANALYSIS=1;
      RUN_POLY14_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rtt|--run-ttv)
      RUN_ANALYSIS=1;
      RUN_TTV_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rbv1|--run-hbov1)
      RUN_ANALYSIS=1;
      RUN_HBOV1_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rbv0|--run-hbovnot1)
      RUN_ANALYSIS=1;
      RUN_HBOVNOT1_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rhbv|--run-hbv)
      RUN_ANALYSIS=1;
      RUN_HBV_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rhpv|--run-hpv)
      RUN_ANALYSIS=1;
      RUN_HPV_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rvar|--run-varv)
      RUN_ANALYSIS=1;
      RUN_VARV_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rsv40|--run-sv40)
      RUN_ANALYSIS=1;
      RUN_SV40_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rcuta|--run-cuta)
      RUN_ANALYSIS=1;
      RUN_CUTA_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rherv|--run-herv)
      RUN_ANALYSIS=1;
      RUN_HERV_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rya|--run-cy-align)
      RUN_ANALYSIS=1;
      RUN_CY_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -ryq|--run-cy-quant)
      RUN_ANALYSIS=1;
      RUN_CY_QUANT_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -vis|--visual-align)
      RUN_VISUAL_ALIGN=1;
      SHOW_HELP=0;
      shift
    ;;
    -covl|--coverage-latex)
      RUN_COVERAGE_TABLE=1;
      SHOW_HELP=0;
      shift
    ;;
    -covc|--coverage-csv)
      RUN_COVERAGE_TABLE_CSV=1;
      SHOW_HELP=0;
      shift
    ;;
    -rpro|--run-profiles)
      RUN_PROFILES_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rda|--run-de-novo)
      RUN_ANALYSIS=1;
      RUN_DE_NOVO_ASSEMBLY=1;
      SHOW_HELP=0;
      shift
    ;;
    -rhyb|--run-hybrid)
      RUN_ANALYSIS=1;
      RUN_HYBRID=1;
      SHOW_HELP=0;
      shift
    ;;
    -enc|--encrypt)
      RUN_ENCRYPT=1;
      SHOW_HELP=0;
      shift
    ;;
    -dec|--decrypt)
      RUN_DECRYPT=1;
      SHOW_HELP=0;
      shift
    ;;
    -all|--run-all)
      INSTALL=1;
      BUILD_VDB_ALL=1;
      BUILD_UDB=1;
      GEN_ADAPTERS=1;
      GET_PHIX=1;
      GET_MITO=1;
      GET_CY=1;
      RUN_ANALYSIS=1;
      #
      RUN_META_ON=1;
      RUN_PROFILES_ON=1;
      RUN_META_NON_VIRAL_ON=1;
      RUN_B19_ON=1;
      RUN_HV1_ON=1;
      RUN_HV2_ON=1;
      RUN_HV3_ON=1;
      RUN_HV4_ON=1;
      RUN_HV5_ON=1;
      RUN_HV6_ON=1;
      RUN_HV6A_ON=1;
      RUN_HV6B_ON=1;
      RUN_HV7_ON=1;
      RUN_HV8_ON=1;
      RUN_POLY1_ON=1;
      RUN_POLY2_ON=1;
      RUN_POLY3_ON=1;
      RUN_POLY4_ON=1;
      RUN_POLY5_ON=1;
      RUN_POLY6_ON=1;
      RUN_POLY7_ON=1;
      RUN_POLY8_ON=1;
      RUN_POLY9_ON=1;
      RUN_POLY10_ON=1;
      RUN_POLY11_ON=1;
      RUN_POLY12_ON=1;
      RUN_POLY13_ON=1;
      RUN_POLY14_ON=1;
      RUN_TTV_ON=1;
      RUN_HBOV1_ON=1;
      RUN_HBOVNOT1_ON=1;
      RUN_HBV_ON=1;
      RUN_HPV_ON=1;
      RUN_VARV_ON=1;
      RUN_SV40_ON=1;
      RUN_CUTA_ON=1;
      RUN_HERV_ON=1;
      RUN_MITO_ON=1;
      RUN_MITO_DAMAGE_ON=1;
      RUN_CY_ON=1;
      RUN_CY_QUANT_ON=1;
      RUN_DE_NOVO_ASSEMBLY=1;
      RUN_HYBRID=1;
      SHOW_HELP=0;
      shift
    ;;
    -*) # unknown option with small
    echo "Invalid arg ($1)!";
    echo "For help, try: ./TRACESPipe.sh -h"
    exit 1;
    ;;
  esac
  done
#
set -- "${POSITIONAL[@]}" # restore positional parameters
#
# ==============================================================================
# HELP
#
if [ "$SHOW_HELP" -eq "1" ];
  then
  echo "                                                                "
  echo -e "\e[34m                                                         "
  echo "         ████████╗ ██████╗   █████╗   ██████╗ ███████╗ ███████╗   "
  echo "         ╚══██╔══╝ ██╔══██╗ ██╔══██╗ ██╔════╝ ██╔════╝ ██╔════╝   "
  echo "            ██║    ██████╔╝ ███████║ ██║      █████╗   ███████╗   "
  echo "            ██║    ██╔══██╗ ██╔══██║ ██║      ██╔══╝   ╚════██║   "
  echo "            ██║    ██║  ██║ ██║  ██║ ╚██████╗ ███████╗ ███████║   "
  echo "            ╚═╝    ╚═╝  ╚═╝ ╚═╝  ╚═╝  ╚═════╝ ╚══════╝ ╚══════╝   "
  echo "                                                                  "
  echo "                             P I P E L I N E                            "
  echo "                                                                "
  echo -e "    \e[32mA Next-generation sequencing pipeline for identification, assembly,\e[0m" 
  echo -e "    \e[32mand analysis of viral and human-host genomes at multi-organ level\e[0m."
  echo "                                                                "
  echo -e "\e[93m    Usage: ./TRACESPipe.sh [options]                             \e[0m"
  echo "                                                                   "
  echo "    -h,     --help            Show this help message and exit,     "
  echo "    -v,     --version         Show the version and some information,  "
  echo "    -f,     --force           Force running and overwrite of files,  "
  echo "    -flog,  --flush-logs      Flush logs (delete logs),              "
  echo "    -i,     --install         Installation of all the tools,       "
  echo "    -up,    --update          Update all the tools in TRACESPipe,  "
  echo "                                                                   "
  echo "    -gmt,   --get-max-threads Get the number of maximum machine threads,"
  echo "    -t <THREADS>, --threads <THREADS>                              "
  echo "                              Number of threads to use,            "
  echo "                                                                   "
  echo "    -dec,   --decrypt         Decrypt (all files in ../encrypted_data), "
  echo "    -enc,   --encrypt         Encrypt (all files in ../to_encrypt_data),"
  echo "                                                                   "
  echo "    -vdb,   --build-viral     Build viral database (all) [Recommended], "
  echo "    -vdbr,  --build-viral-r   Build viral database (references only),  "
  echo "    -udb,   --build-unviral   Build non viral database (control),  "
  echo "                                                                   "
  echo "    -afs <FASTA>, --add-fasta <FASTA>                               "
  echo "                              Add a FASTA sequence to the VDB.fa,  "
  echo "    -aes <ID>, --add-extra-seq <ID>                                "
  echo "                              Add extra sequence to the VDB.fa,    "
  echo "    -gx,    --get-extra-vir   Downloads/appends (VDB) extra viral seq, "
  echo "                                                                   "
  echo "    -gad,   --gen-adapters    Generate FASTA file with adapters,   "
  echo "    -gp,    --get-phix        Extracts PhiX genomes (Needs viral DB),  "
  echo "    -gm,    --get-mito        Downloads human Mitochondrial genome,"
  echo "                                                                   "
  echo "    -cmt <ID>, --change-mito <ID>                                  "
  echo "                              Set any Mitochondrial genome by ID,  "
  echo "                                                                   "
  echo "    -gy,    --get-y-chromo    Downloads human Y-chromosome,        "
  echo "    -gax,   --get-all-aux     Runs -gad -gp -gm -gy,               "
  echo "                                                                   "
  echo "    -cbn,   --create-blast-db It creates a nucleotide blast database, "
  echo "    -ubn,   --update-blast-db It updates a nucleotide blast database, "
  echo "    -sfs <FASTA>, --search-blast-db <FASTA>                           "
  echo "                              It blasts the nucleotide (nt) blast DB, "
  echo "                                                                   "
  echo "    -rdup,  --remove-dup      Remove duplications (e.g. PCR dup),  "
  echo "                                                                   "
  echo "    -iss <SIZE>, --inter-sim-size <SIZE>                                  "
  echo "                              Inter-genome similarity top size (control), "
  echo "                                                                   "
  echo "    -rpro,  --run-profiles    Run complexity and relative profiles (control), "
  echo "                                                                   "
  echo "    -rm,    --run-meta        Run viral metagenomic identification,    "
  echo "    -ro,    --run-meta-nv     Run NON-viral metagenomic identification,"
  echo "                                                                  "
  echo "    -mis <VALUE>, --min-similarity <VALUE>                         "
  echo "                              Minimum similarity value to consider the "
  echo "                              sequence for alignment-consensus (filter), "
  echo "                                                                       "
  echo "    -top <VALUE>, --view-top <VALUE>                                   "
  echo "                              Display the top <VALUE> with the highest "
  echo "                              similarity (by descending order),        "
  echo "                                                                       "
  echo "    -rava,  --run-all-v-alig  Run all viral align/sort/consensus seqs, "
  echo "                                                                       "
  echo "    -rb19,  --run-b19         Run B19  align and consensus seq,        "
  echo "    -rh1,   --run-hv1         Run HHV1   align and consensus seq,    "
  echo "    -rh2,   --run-hv2         Run HHV2   align and consensus seq,    "
  echo "    -rh3,   --run-hv3         Run HHV3   align and consensus seq,    "
  echo "    -rh4,   --run-hv4         Run HHV4   align and consensus seq,    "
  echo "    -rh5,   --run-hv5         Run HHV5   align and consensus seq,    "
  echo "    -rh6,   --run-hv6         Run HHV6   align and consensus seq,    "
  echo "    -rh6a,  --run-hv6a        Run HHV6A  align and consensus seq,    "
  echo "    -rh6b,  --run-hv6b        Run HHV6B  align and consensus seq,    "
  echo "    -rh7,   --run-hv7         Run HHV7   align and consensus seq,    "
  echo "    -rh8,   --run-hv8         Run HHV8   align and consensus seq,    "
  echo "    -rh8,   --run-hv8         Run HHV8   align and consensus seq,    "
 
  echo "    -rp1,   --run-poly1       Run Polyoma 1  align and consensus seq, "
  echo "    -rp2,   --run-poly2       Run Polyoma 2  align and consensus seq, "
  echo "    -rp3,   --run-poly3       Run Polyoma 3  align and consensus seq, "
  echo "    -rp4,   --run-poly4       Run Polyoma 4  align and consensus seq, "
  echo "    -rp5,   --run-poly5       Run Polyoma 5  align and consensus seq, "
  echo "    -rp6,   --run-poly6       Run Polyoma 6  align and consensus seq, "
  echo "    -rp7,   --run-poly7       Run Polyoma 7  align and consensus seq, "
  echo "    -rp8,   --run-poly8       Run Polyoma 8  align and consensus seq, "
  echo "    -rp9,   --run-poly9       Run Polyoma 9  align and consensus seq, "
  echo "    -rp10,  --run-poly10      Run Polyoma 10 align and consensus seq, "
  echo "    -rp11,  --run-poly11      Run Polyoma 11 align and consensus seq, "
  echo "    -rp12,  --run-poly12      Run Polyoma 12 align and consensus seq, "
  echo "    -rp13,  --run-poly13      Run Polyoma 13 align and consensus seq, "
  echo "    -rp14,  --run-poly14      Run Polyoma 14 align and consensus seq, "
 
  echo "    -rtt,   --run-ttv         Run TTV   align and consensus seq,    "
  echo "    -rbv1,  --run-hbov1       Run HBoV1 align and consensus seq,    "
  echo "    -rbv0,  --run-hbovnot1    Run HBoV (2,3,...) align/consensus seq, "
  echo "    -rhbv,  --run-hbv         Run HBV   align and consensus seq,    "
  echo "    -rhpv,  --run-hpv         Run HPV   align and consensus seq,    "
  echo "    -rvar,  --run-varv        Run VARV  align and consensus seq,    "
  echo "    -rsv40, --run-sv40        Run Simian 40 align and consensus seq,  "
  echo "    -rcuta, --run-cuta        Run Cutavirys align and consensus seq,  "
  echo "    -rherv, --run-herv        Run H Endo Retro align and consensus seq,  "
  echo "                                                                  "
  echo "    -rsr <ID>, --run-specific <ID/PATTERN>                        "
  echo "                              Run specific reference align/consensus, "
  echo "    -rsx <ID>, --run-extreme <ID/PATTERN>                            "
  echo "                              Run specific reference align/consensys"
  echo "                              using extreme sensitivity,            "
  echo "                                                                 "
  echo "    -rmt,   --run-mito        Run Mito align and consensus seq,   "
  echo "    -rmtd,  --run-mito-dam    Run Mito damage only,               "
  echo "                                                                 "
  echo "    -rya,   --run-cy-align    Run CY align and consensus seq,    "
  echo "    -ryq,   --run-cy-quant    Estimate the quantity of CY DNA,    "
  echo "                                                                  "
  echo "    -rda,   --run-de-novo     Run de-novo assembly,               "
  echo "                                                                  "
  echo "    -rhyb,  --run-hybrid      Run hybrid assembly (align/de-novo), "
  echo "                                                                  "
  echo "    -vis,   --visual-align    Run Visualization tool for alignments, "
  echo "    -covl,  --coverage-latex  Run coverage table in Latex format,   "
  echo "    -covc,  --coverage-csv    Run coverage table in CSV format,    "
  echo "    -covp <NAME>, --coverage-profile <BED_NAME_FILE>                      "
  echo "                              Run coverage profile for specific BED file, "
  echo "    -cmax <MAX>,  --max-coverage <MAX_COVERAGE>                           "
  echo "                              Maximum depth coverage (depth normalization), "
  echo "                                                                  "
  echo "    -ra,    --run-analysis    Run data analysis,                   "
  echo "    -all,   --run-all         Run all the options.                 "
  echo "                                                                "
  echo -e "\e[93m    Example: ./TRACESPipe.sh --run-meta --run-b19 --run-mito \e[0m"
  echo "                                                                "
  echo "    Add the file meta_info.txt at ../meta_data/ folder. Example:      "
  echo "    meta_info.txt -> 'organ:reads_forward.fa.gz:reads_reverse.fa.gz'  "
  echo "    The reads must be GZIPed in the ../input_data/ folder.            "
  echo "    The output results are at ../output_data/ folder.                 "
  echo "                                                                "
  echo -e "\e[32m    Contact: projectraces@gmail.com                  \e[0m"
  echo "                                                                "
  exit 1;
  fi
#
# ==============================================================================
# VERSION
#
if [ "$SHOW_VERSION" -eq "1" ];
  then
  echo "                                                                      ";
  echo "                              TRACESPipe                              ";
  echo "                                                                      ";
  echo "                            Version: 1.0.1                            ";
  echo "                                                                      ";
  echo "                      Department of Virology and                      ";
  echo "                   Department of Forensic Medicine,                   ";
  echo "                   University of Helsinki, Finland.                   ";
  echo "                                  &                                   ";
  echo "                             IEETA/DETI,                              ";
  echo "                    University of Aveiro, Portugal.                   ";
  echo "                                                                      ";
  echo "                        projectraces@gmail.com                        ";
  echo "                                                                      ";
  exit 0;
  fi
#
# ==============================================================================
# MAKE SURE FOLDERS EXIST
#
mkdir -p ../logs/
mkdir -p ../output_data/
#
# ==============================================================================
#
if [[ "$FLUSH_LOGS" -eq "1" ]];
  then
  echo -e "\e[34m[TRACESPipe]\e[32m Flushing logs ...\e[0m";
  rm -f ../logs/Log-stdout-*.txt
  rm -f ../logs/Log-stderr-*.txt
  echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
  fi
#  
# ==============================================================================
#
if [[ "$RUN_VISUAL_ALIGN" -eq "1" ]];
  then
  ./TRACES_run_visual_alignment.sh
  exit 0;
  fi
#
# ==============================================================================
#
if [[ "$RUN_COVERAGE_TABLE" -eq "1" ]];
  then
  ./TRACES_coverage_table.sh
  exit 0;
  fi
#
# ==============================================================================
#
if [[ "$RUN_COVERAGE_TABLE_CSV" -eq "1" ]];
  then
  ./TRACES_coverage_table_csv.sh
  exit 0;
  fi
#
# ==============================================================================
#
if [[ "$CREATE_BLAST_DB" -eq "1" ]];
  then
  ./TRACES_blastn_create_n_db.sh
  exit 0;
  fi
#
# ==============================================================================
#
if [[ "$UPDATE_BLAST_DB" -eq "1" ]];
  then
  ./TRACES_blastn_update_n_db.sh
  exit 0;
  fi
#
# ==============================================================================
#
if [[ "$SEARCH_BLAST_DB" -eq "1" ]];
  then
  ./TRACES_blastn_n_db.sh $BLAST_QUERY
  exit 0;
  fi
#
# ==============================================================================
#
if [[ "$RUN_COVERAGE_PROFILE" -eq "1" ]];
  then
  CHECK_E_FILE $COVERAGE_NAME
  ./TRACES_project_coordinates.sh $COVERAGE_NAME $COVERAGE_MAX > x.projectd.profile;
  #XXX: OPTIONALLY, LOW-PASS FILTER CAN BE APPLIED HERE (for larger sequences)!
gnuplot << EOF
    reset
    set terminal pdfcairo enhanced color font 'Verdana,12'
    set output "$COVERAGE_NAME.pdf"
    set style line 101 lc rgb '#000000' lt 1 lw 4
    set border 3 front ls 101
    set tics nomirror out scale 0.75
    set format '%g'
    set size ratio 0.2
    set key outside horiz center top
    set yrange [:]
    set xrange [:]
    set xtics auto
    set grid
    set ylabel "Depth"
    set xlabel "Position"
    set border linewidth 1.5
    set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 5 ps 0.4 # --- blue
    set style line 2 lc rgb '#228B22' lt 1 lw 2 pt 6 ps 0.4 # --- green
    set style line 3 lc rgb '#dd181f' lt 1 lw 4 pt 7 ps 0.4 # --- ?
    set style line 4 lc rgb '#4d1811' lt 1 lw 4 pt 8 ps 0.4 # --- ?
    set style line 5 lc rgb '#1d121f' lt 1 lw 4 pt 9 ps 0.4 # --- ?
    plot "x.projectd.profile" using 2 with lines ls 2
EOF
  exit 0;
  fi
#
# ==============================================================================
#
if [[ "$GET_THREADS" -eq "1" ]];
  then
  ./TRACES_get_max_threads.sh
  exit 1;
  fi
#
# ==============================================================================
#
if [[ "$THREADS" -eq "0" ]];
  then
  THREADS=`./TRACES_get_max_threads.sh |awk '{print $8;}'`;
  echo -e "\e[34m[TRACESPipe]\e[32m Running with $THREADS threads.\e[0m";
  else
  echo -e "\e[34m[TRACESPipe]\e[32m Running with $THREADS threads.\e[0m";
  fi
#
# ==============================================================================
#
if [[ "$INSTALL" -eq "1" ]];
  then
  ./TRACES_install.sh
  fi
#
# ==============================================================================
#
if [[ "$UPDATE" -eq "1" ]];
  then
  ./TRACES_update.sh
  fi 
#
# ==============================================================================
#
if [[ "$RUN_ENCRYPT" -eq "1" ]];
  then
  ./TRACES_encrypt_data.sh
  fi
#
# ==============================================================================
#
if [[ "$RUN_DECRYPT" -eq "1" ]];
  then
  ./TRACES_decrypt_data.sh
  fi
#
# ==============================================================================
#
if [[ "$BUILD_VDB_REF" -eq "1" ]];
  then
  gto_build_dbs.sh -vi
  gunzip VDB.fa.gz
  fi
#
# ==============================================================================
#
if [[ "$BUILD_VDB_ALL" -eq "1" ]];
  then
  ./TRACES_download_db.sh
  fi
#
# ==============================================================================
#
if [[ "$BUILD_UDB" -eq "1" ]];
  then
  gto_build_dbs.sh -ba -ar -pr -fu -pl -in -mi -ps 
  zcat BDB.fa.gz ADB.fa.gz PDB.fa.gz FDB.fa.gz TDB.fa.gz TDB.fa.gz IDB.fa.gz MTDB.fa.gz PLDB.fa.gz > DB.fa
  # ADDITIONAL: -vm -vo
  fi
#
# ==============================================================================
#
if [[ "$GEN_ADAPTERS" -eq "1" ]];
  then
  ./TRACES_generate_adapters.sh
  fi
#
# ==============================================================================
#
if [[ "$GET_PHIX" -eq "1" ]];
  then
  ./TRACES_get_phix.sh
  fi
#
# ==============================================================================
#
if [[ "$GET_MITO" -eq "1" ]];
  then
  ./TRACES_get_mito.sh
  fi
#
# ==============================================================================
#
if [[ "$RUN_CHANGE_MT" -eq "1" ]];
  then
  ./TRACES_change_mito.sh $NEW_MT
  fi
#
# ==============================================================================
#
if [[ "$GET_CY" -eq "1" ]];
  then
  ./TRACES_get_cy.sh
  fi
#
# ==============================================================================
#
if [[ "$ADD_FASTA" -eq "1" ]];
  then
  CHECK_VDB;
  ./TRACES_add_fasta_to_database.sh VDB.fa $NEW_FASTA
  fi
#
# ==============================================================================
#
if [[ "$GET_EXTRA" -eq "1" ]];
  then
  CHECK_ENRICH;
  ./TRACES_get_enriched_sequences.sh VDB.fa
  fi
#
# ==============================================================================
#
if [[ "$ADD_EXTRA_SEQ" -eq "1" ]];
  then
  CHECK_VDB;
  ./TRACES_get_extra_seq.sh VDB.fa $NEW_SEQ_ID
  fi
#
# ==============================================================================
#
if [[ "$RUN_ANALYSIS" -eq "1" ]];
  then
  #
  CHECK_FILTERING_SYSTEM_FILES;
  CHECK_META_INFO;
  #
  mapfile -t READS < ../meta_data/meta_info.txt
  #
  for read in "${READS[@]}" # 
    do
    #
    ORGAN_T=`echo $read | tr ':' '\t' | awk '{ print $1 }'`;
    SPL_Forward=`echo $read | tr ':' '\t' | awk '{ print $2 }'`;
    SPL_Reverse=`echo $read | tr ':' '\t' | awk '{ print $3 }'`;
    #
    CHECK_GZIP_FILES $SPL_Forward $SPL_Reverse;
    #
    echo -e "\e[34m[TRACESPipe]\e[93m Organ=$ORGAN_T Forward=$SPL_Forward Reverse=$SPL_Reverse\e[0m";
    #
    rm -f FW_READS.fq.gz RV_READS.fq.gz
    echo -e "\e[34m[TRACESPipe]\e[32m Copping an instance of the files ...\e[0m";
    cp ../input_data/$SPL_Forward FW_READS.fq.gz;
    cp ../input_data/$SPL_Reverse RV_READS.fq.gz;
    echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
    #
    # ========================================================================
    # TRIM AND FILTER READS TRIMMOMATIC
    #
    CHECK_ADAPTERS;
    touch o_fw_pr.fq o_fw_unpr.fq o_rv_pr.fq o_rv_unpr.fq;
    #
    echo -e "\e[34m[TRACESPipe]\e[32m Trimming and filtering with Trimmomatic ...\e[0m";
    ./TRACES_trim_filter_reads.sh $THREADS 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
    echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
    #
    # THE OUTPUT OF TRIMMING IS:
    # o_fw_pr.fq  o_fw_unpr.fq  o_rv_pr.fq  o_rv_unpr.fq
    #
    # ========================================================================== 
    # METAGENOMICS USING ONLY A VIRAL DATABASE
    #
    if [[ "$RUN_META_ON" -eq "1" ]];
      then
      CHECK_VDB;
      CHECK_PHIX;
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Removing PhiX from the samples with MAGNET ...\e[0m";
      #./TRACES_remove_phix.sh $THREADS 
      # XXX: REMOVE FROM THE DB PHIX OR CORRECT MAGNET AT READ LEVEL
      cp o_fw_pr.fq NP-o_fw_pr.fq;
      cp o_fw_unpr.fq NP-o_fw_unpr.fq;
      cp o_rv_pr.fq NP-o_rv_pr.fq;
      cp o_rv_unpr.fq NP-o_rv_unpr.fq;      
      # IT IS USED ONLY FOR FALCON
      #
      # fastq_pair test_R1.fastq test_R2.fastq: [needs adaptation]
      # IF you want to remove Phix also before assembly
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Running viral metagenomic analysis with FALCON-meta ...\e[0m";
      mkdir -p ../output_data/TRACES_results
      ./TRACES_metagenomics_viral.sh $ORGAN_T VDB.fa 10000 $THREADS $TSIZE 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
      mv $ORGAN_T.svg ../output_data/TRACES_results/
      mv $ORGAN_T-HEAT.svg ../output_data/TRACES_results/
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Finding the best references ...\e[0m";
      #
      cp top-$ORGAN_T.csv ../output_data/TRACES_results/
      #
      rm -f ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt;
      for VIRUS in "${VIRUSES[@]}"
        do
        ./TRACES_get_best_$VIRUS.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
        done
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      #
      fi
    #
    # ==========================================================================
    # COMPLEXITY PROFILES FOR TOP IDENTIFIED VIRAL SEQUENCES
    #
    if [[ "$RUN_PROFILES_ON" -eq "1" ]];
      then
      #	      
      CHECK_VDB;
      #
      mkdir -p ../output_data/TRACES_results
      mkdir -p ../output_data/TRACES_results/profiles/
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Building complexity profiles with gto ...\e[0m";
      cat NP-o_fw_pr.fq NP-o_fw_unpr.fq NP-o_rv_pr.fq NP-o_rv_unpr.fq > P_TRACES_sample_reads.fq
      ./TRACES_profiles.sh GIS-$ORGAN_T VDB.fa P_TRACES_sample_reads.fq $ORGAN_T 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      fi
    #
    # ========================================================================== 
    # METAGENOMICS WITH ALL NON VIRAL
    #
    if [[ "$RUN_META_NON_VIRAL_ON" -eq "1" ]];
      then
      #
      CHECK_DB;
      #	
      echo -e "\e[34m[TRACESPipe]\e[32m Running NON viral metagenomic analysis with FALCON ...\e[0m";
      ./TRACES_metagenomics.sh $ORGAN_T DB.fa 12000 $THREADS $TSIZE 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
      mkdir -p ../output_data/TRACES_results;
      #rm -f ../output_data/TRACES_results/*
      mv NV-$ORGAN_T.svg ../output_data/TRACES_results/
      mv NV-$ORGAN_T-HEAT.svg ../output_data/TRACES_results/
      mv REPORT_META_NON_VIRAL_$ORGAN_T.txt ../output_data/TRACES_results/
      cp top-non-viral-$ORGAN_T.csv  ../output_data/TRACES_results/
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      fi
    #
    # ==============================================================================
    #
    if [[ "$VIEW_TOP" -eq "1" ]];
      then
      CHECK_E_FILE "../output_data/TRACES_results/top-$ORGAN_T.csv";
      F_TOP_SIZE=`wc -l ../output_data/TRACES_results/top-$ORGAN_T.csv | awk '{ print $1; }'`;
      if [[ "$F_TOP_SIZE" -eq "0" ]];
        then
	echo -e "\e[34m[TRACESPipe]\e[31m Empty top-$ORGAN_T.csv results!\e[0m";
	else
        head -n $VTOP_SIZE ../output_data/TRACES_results/top-$ORGAN_T.csv
        fi
      fi
    #
    # ==========================================================================
    # RUN SPECIFIC SPECIFIC ALIGN/CONSENSUS
    #
    if [[ "$RUN_SPECIFIC" -eq "1" ]];
      then
      echo -e "\e[34m[TRACESPipe]\e[32m Aligning reads to specific viral ref(s) with pattern \"$SPECIFIC_ID\" using bowtie2 ...\e[0m";
      #
      CHECK_VDB;
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Extracting sequence with pattern \"$SPECIFIC_ID\" from VDB.fa ...\e[0m";
      gto_fasta_extract_read_by_pattern -p "$SPECIFIC_ID" < VDB.fa > SPECIFIC-$SPECIFIC_ID.fa
      echo -e "\e[34m[TRACESPipe]\e[32m Aligning ... \e[0m";
      ./TRACES_viral_align_reads.sh SPECIFIC-$SPECIFIC_ID.fa $ORGAN_T $SPECIFIC_ID $THREADS $REMOVE_DUPLICATIONS 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
      ./TRACES_viral_consensus.sh SPECIFIC-$SPECIFIC_ID.fa viral_aligned_sorted-$ORGAN_T-$SPECIFIC_ID.bam $ORGAN_T $SPECIFIC_ID 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
      mkdir -p ../output_data/TRACES_specific_alignments;
      #rm -f ../output_data/TRACES_specific_alignments/*
      cp SPECIFIC-$SPECIFIC_ID.fa ../output_data/TRACES_specific_alignments/
      cp SPECIFIC-$SPECIFIC_ID.fa.fai ../output_data/TRACES_specific_alignments/
      mv viral_aligned_sorted-$ORGAN_T-$SPECIFIC_ID.bam ../output_data/TRACES_specific_alignments/
      mv viral_aligned_sorted-$ORGAN_T-$SPECIFIC_ID.bam.bai ../output_data/TRACES_specific_alignments/
      mkdir -p ../output_data/TRACES_specific_consensus;
      #rm -f ../output_data/TRACES_specific_consensus/*
      mv $SPECIFIC_ID-consensus-$ORGAN_T.fa ../output_data/TRACES_specific_consensus/
      mkdir -p ../output_data/TRACES_specific_bed;
      #rm -f ../output_data/TRACES_specific_bed/*
      mv $SPECIFIC_ID-calls-$ORGAN_T.bed ../output_data/TRACES_specific_bed/
      mv $SPECIFIC_ID-coverage-$ORGAN_T.bed ../output_data/TRACES_specific_bed/
      mv $SPECIFIC_ID-zero-coverage-$ORGAN_T.bed ../output_data/TRACES_specific_bed/
      mkdir -p ../output_data/TRACES_specific_statistics;
      echo -e "\e[34m[TRACESPipe]\e[32m Calculating coverage ...\e[0m";
      ./TRACES_overall_specific.sh $SPECIFIC_ID $ORGAN_T
      C_BREADTH=`cat ../output_data/TRACES_specific_statistics/$SPECIFIC_ID-total-horizontal-coverage-$ORGAN_T.txt`;
      C_DEPTH=`cat ../output_data/TRACES_specific_statistics/$SPECIFIC_ID-total-depth-coverage-$ORGAN_T.txt`;
      echo -e "\e[34m[TRACESPipe]\e[1m Breadth (H) coverage: $C_BREADTH \e[0m";
      echo -e "\e[34m[TRACESPipe]\e[1m Depth-x (V) coverage: $C_DEPTH \e[0m";
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      fi
    #
    # ==========================================================================
    # RUN SPECIFIC SPECIFIC ALIGN/CONSENSUS WITH EXTREME HIGH SENSITIVITY
    #
    if [[ "$RUN_SPECIFIC_SENSITIVE" -eq "1" ]];
      then
      echo -e "\e[34m[TRACESPipe]\e[32m Aligning reads to specific viral ref(s) with pattern \"$SPECIFIC_ID\" using bowtie2 with EXTREME sensitivity ...\e[0m";
      #
      CHECK_VDB;
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Extracting sequence with pattern \"$SPECIFIC_ID\" from VDB.fa ...\e[0m";
      gto_fasta_extract_read_by_pattern -p "$SPECIFIC_ID" < VDB.fa > SPECIFIC-$SPECIFIC_ID.fa
      echo -e "\e[34m[TRACESPipe]\e[32m Aligning ...\e[0m";
      ./TRACES_viral_sensitive_align_reads.sh SPECIFIC-$SPECIFIC_ID.fa $ORGAN_T $SPECIFIC_ID $THREADS $REMOVE_DUPLICATIONS 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
      ./TRACES_viral_consensus.sh SPECIFIC-$SPECIFIC_ID.fa viral_aligned_sorted-$ORGAN_T-$SPECIFIC_ID.bam $ORGAN_T $SPECIFIC_ID 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
      mkdir -p ../output_data/TRACES_specific_alignments;
      #rm -f ../output_data/TRACES_specific_alignments/*
      cp SPECIFIC-$SPECIFIC_ID.fa ../output_data/TRACES_specific_alignments/
      cp SPECIFIC-$SPECIFIC_ID.fa.fai ../output_data/TRACES_specific_alignments/
      mv viral_aligned_sorted-$ORGAN_T-$SPECIFIC_ID.bam ../output_data/TRACES_specific_alignments/
      mv viral_aligned_sorted-$ORGAN_T-$SPECIFIC_ID.bam.bai ../output_data/TRACES_specific_alignments/
      mkdir -p ../output_data/TRACES_specific_consensus;
      #rm -f ../output_data/TRACES_specific_consensus/*
      mv $SPECIFIC_ID-consensus-$ORGAN_T.fa ../output_data/TRACES_specific_consensus/
      mkdir -p ../output_data/TRACES_specific_bed;
      #rm -f ../output_data/TRACES_specific_bed/*
      mv $SPECIFIC_ID-calls-$ORGAN_T.bed ../output_data/TRACES_specific_bed/
      mv $SPECIFIC_ID-coverage-$ORGAN_T.bed ../output_data/TRACES_specific_bed/
      mv $SPECIFIC_ID-zero-coverage-$ORGAN_T.bed ../output_data/TRACES_specific_bed/
      mkdir -p ../output_data/TRACES_specific_statistics;
      echo -e "\e[34m[TRACESPipe]\e[32m Calculating coverage ...\e[0m";
      ./TRACES_overall_specific.sh $SPECIFIC_ID $ORGAN_T
      C_BREADTH=`cat ../output_data/TRACES_specific_statistics/$SPECIFIC_ID-total-horizontal-coverage-$ORGAN_T.txt`;
      C_DEPTH=`cat ../output_data/TRACES_specific_statistics/$SPECIFIC_ID-total-depth-coverage-$ORGAN_T.txt`;
      echo -e "\e[34m[TRACESPipe]\e[1m Breadth (H) coverage: $C_BREADTH \e[0m";
      echo -e "\e[34m[TRACESPipe]\e[1m Depth-x (V) coverage: $C_DEPTH \e[0m";
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      fi
    #
    # ========================================================================== 
    # DETAILED VIRAL ALIGN/CONSENSUS
    #
    if [[ "$RUN_B19_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "B19" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_HV1_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HV1" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_HV2_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HV2" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_HV3_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HV3" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_HV4_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HV4" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_HV5_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HV5" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_HV6_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HV6" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_HV6A_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HV6A" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_HV6B_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HV6B" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_HV7_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HV7" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_HV8_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HV8" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_POLY1_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY1" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_POLY2_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY2" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_POLY2_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY2" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_POLY3_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY3" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_POLY4_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY4" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_POLY5_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY5" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_POLY6_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY6" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_POLY7_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY7" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_POLY8_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY8" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_POLY9_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY9" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_POLY10_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY10" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_POLY11_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY11" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_POLY12_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY12" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_POLY13_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY13" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_POLY14_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY14" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_TTV_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "TTV" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_HBOV1_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HBOV1" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_HBOVNOT1_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HBOVNOT1" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_HBV_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HBV" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_HPV_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HPV" "$MINIMAL_SIMILARITY_VALUE";
      fi
    # 
    if [[ "$RUN_VARV_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "VARV" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_SV40_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "SV40" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_CUTA_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "CUTA" "$MINIMAL_SIMILARITY_VALUE";
      fi
    #
    if [[ "$RUN_HERV_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HERV" "$MINIMAL_SIMILARITY_VALUE";
      fi
    # 
    # ==========================================================================
    # MITOCHONDRIAL GENOME ALIGN AND CONSENSUS
    #
    if [[ "$RUN_MITO_ON" -eq "1" ]];
      then
      #
      CHECK_MT_DNA;
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Aligning reads to mitochondrial ref with bowtie2 ...\e[0m";
      ./TRACES_mt_align_reads.sh mtDNA.fa $ORGAN_T $THREADS $REMOVE_DUPLICATIONS 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
      ./TRACES_mt_consensus.sh mtDNA.fa mt_aligned_sorted-$ORGAN_T.bam $ORGAN_T 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m"
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Storing files ...\e[0m"
      mkdir -p ../output_data/TRACES_mtdna_alignments;
      #rm -f ../output_data/TRACES_mtdna_alignments/*
      cp mtDNA.fa ../output_data/TRACES_mtdna_alignments/
      cp mtDNA.fa.fai ../output_data/TRACES_mtdna_alignments/
      mv mt_aligned_sorted-$ORGAN_T.bam ../output_data/TRACES_mtdna_alignments/
      mv mt_aligned_sorted-$ORGAN_T.bam.bai ../output_data/TRACES_mtdna_alignments/
      mkdir -p ../output_data/TRACES_mtdna_consensus;
      #rm -f ../output_data/TRACES_mtdna_consensus/*
      mv mt-consensus-$ORGAN_T.fa ../output_data/TRACES_mtdna_consensus/
      mkdir -p ../output_data/TRACES_mtdna_bed;
      #rm -f ../output_data/TRACES_mtdna_bed/*
      mv mt-calls-$ORGAN_T.bed ../output_data/TRACES_mtdna_bed/
      mv mt-coverage-$ORGAN_T.bed ../output_data/TRACES_mtdna_bed/
      mv mt-zero-coverage-$ORGAN_T.bed ../output_data/TRACES_mtdna_bed/
      mkdir -p ../output_data/TRACES_mtdna_statistics;
      echo -e "\e[34m[TRACESPipe]\e[32m Calculating coverage ...\e[0m";
      ./TRACES_overall_mtdna.sh $ORGAN_T
      C_BREADTH=`cat ../output_data/TRACES_mtdna_statistics/mt-total-horizontal-coverage-$ORGAN_T.txt`;
      C_DEPTH=`cat ../output_data/TRACES_mtdna_statistics/mt-total-depth-coverage-$ORGAN_T.txt`;
      echo -e "\e[34m[TRACESPipe]\e[35m Breadth (H) coverage: $C_BREADTH \e[0m";
      echo -e "\e[34m[TRACESPipe]\e[35m Depth-x (V) coverage: $C_DEPTH \e[0m";
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m"
      fi
    # ========================================================================
    # TRIM AND FILTER READS ADAPTERREMOVAL
    #
    if [[ "$RUN_MITO_DAMAGE_ON" -eq "1" ]];
      then
      #	
      CHECK_ADAPTERS_AR;
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Trimming, filtering, and collapsing with AdapterRemoval ...\e[0m";
      AdapterRemoval --threads $THREADS --file1 FW_READS.fq.gz --file2 RV_READS.fq.gz --outputcollapsed reads.fq --trimns --trimqualities --minlength 30 --collapse --adapter-list adapters_ar.fa 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
      echo -e "\e[34m[TRACESPipe]\e[32m Aligning data using bwa ...\e[0m";
      bwa index mtDNA.fa 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
      bwa mem -t $THREADS -I 0 -O 2 -N 0.02 -L 1024 -E 7 mtDNA.fa reads.fq > mt-$ORGAN_T.sam 2>> ../logs/Log-stderr-$ORGAN_T.txt;
      echo -e "\e[34m[TRACESPipe]\e[32m Adapting data with samtools ...\e[0m";
      samtools view -bSh mt-$ORGAN_T.sam > mt-$ORGAN_T.bam
      samtools view -bh -F4 mt-$ORGAN_T.bam > FIL-mt-$ORGAN_T.bam;
      echo -e "\e[34m[TRACESPipe]\e[32m Estimating the damage of mtDNA using mapDamage2 ...\e[0m";
      rm -fr ../output_data/TRACES_mtdna_damage_$ORGAN_T
      #mapDamage --rescale -d ../output_data/TRACES_mtdna_damage_$ORGAN_T -i FIL-$ORGAN_T.bam -r mtDNA.fa;
      mapDamage -d ../output_data/TRACES_mtdna_damage_$ORGAN_T -i FIL-mt-$ORGAN_T.bam -r mtDNA.fa 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m"
      #
      fi
    #
    # ========================================================================== 
    # CY VERYFICATION ALIGN/CONSENSUS
    #
    if [[ "$RUN_CY_ON" -eq "1" ]];
      then
      #
      CHECK_CY_DNA;
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Aligning reads to Y-chromosome ref with bowtie2 ...\e[0m";
      ./TRACES_cy_align_reads.sh cy.fa $ORGAN_T $THREADS 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
      ./TRACES_cy_consensus.sh cy.fa cy_aligned_sorted-$ORGAN_T.bam $ORGAN_T 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
      mkdir -p ../output_data/TRACES_cy_alignments;
      #rm -f ../output_data/TRACES_cy_alignments/*
      cp cy.fa ../output_data/TRACES_cy_alignments/
      cp cy.fa.fai ../output_data/TRACES_cy_alignments/
      mv cy_aligned_sorted-$ORGAN_T.bam ../output_data/TRACES_cy_alignments/
      mv cy_aligned_sorted-$ORGAN_T.bam.bai ../output_data/TRACES_cy_alignments/
      mkdir -p ../output_data/TRACES_cy_consensus;
      #rm -f ../output_data/TRACES_cy_consensus/*
      mv cy-consensus-$ORGAN_T.fa ../output_data/TRACES_cy_consensus/
      mkdir -p ../output_data/TRACES_cy_bed;
      #rm -f ../output_data/TRACES_cy_bed/*
      mv cy-calls-$ORGAN_T.bed ../output_data/TRACES_cy_bed/
      mv cy-coverage-$ORGAN_T.bed ../output_data/TRACES_cy_bed/
      mv cy-zero-coverage-$ORGAN_T.bed ../output_data/TRACES_cy_bed/
      mkdir -p ../output_data/TRACES_cy_statistics;
      ./TRACES_overall_cy.sh $ORGAN_T
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m"
      fi
    #
    # ==========================================================================
    # CY QUANTITY ESTIMATION
    #
    if [[ "$RUN_CY_QUANT_ON" -eq "1" ]];
      then
      #
      CHECK_CY_DNA;
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Estimating the quantity of Y-chromosome ...\e[0m";
      ./TRACES_estimate_cy_quantity.sh $ORGAN_T $THREADS 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
      CY_EST_CALUE=`grep "Dissimilarity" ../output_data/TRACES_results/REP_CY_$ORGAN_T.txt \
      | awk '{ print "BPS: "$6" ; NRC: "$16" "; }'`;
      echo -e "\e[34m[TRACESPipe]\e[36m Estimation: $CY_EST_CALUE \e[0m"; 
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      fi
    #
    # ==========================================================================
    # DE-NOVO ASSEMBLY
    #
    if [[ "$RUN_DE_NOVO_ASSEMBLY" -eq "1" ]];
      then
      echo -e "\e[34m[TRACESPipe]\e[32m Running do-novo DNA assembly with SPAdes ...\e[0m";
      ./TRACES_assemble_all.sh $ORGAN_T $THREADS 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      fi
    #
    #
    # ==========================================================================
    # HYBRID ASSEMBLY: BETWEEN ALIGNMENTS_CONSENSUS & SCAFFOLDS
    #
    # INPUT:
    # -> ALIGNMENTS_CONSENSUS_FASTA
    # -> SCAFFOLDS_MULTI-FASTA 
    # OUTPUT: 
    # ../output_data/TRACES_hybrid_consensus/<VIRUS>-consensus-<ORGAN>.fa
    #
    if [[ "$RUN_HYBRID" -eq "1" ]];
      then
      echo -e "\e[34m[TRACESPipe]\e[32m Running HYBRID assembly ...\e[0m";
      mkdir -p ../output_data/TRACES_hybrid/
      mkdir -p ../output_data/TRACES_hybrid_consensus/
      mkdir -p ../output_data/TRACES_hybrid_alignments/
      mkdir -p ../output_data/TRACES_hybrid_bed/
      SCAFFOLDS_PATH="../output_data/TRACES_denovo_$ORGAN_T/scaffolds.fasta";
      HYBRID_CON_PATH="../output_data/TRACES_hybrid_consensus";
      HYBRID_ALI_PATH="../output_data/TRACES_hybrid_alignments";
      HYBRID_BED_PATH="../output_data/TRACES_hybrid_bed";
      HYBRID_PATH="../output_data/TRACES_hybrid";
      CON_PATH="../output_data/TRACES_viral_consensus";
      #
      mkdir -p $HYBRID_CON_PATH;
      mkdir -p $HYBRID_BED_PATH;
      mkdir -p $HYBRID_ALI_PATH;
      #
      for VIRUS in "${VIRUSES[@]}"
        do
	#cp $CON_PATH/$VIRUS-consensus-$ORGAN_T.fa $ORGAN_T-$VIRUS.fa	
	#
	if [ -f ../output_data/TRACES_viral_alignments/$ORGAN_T-$VIRUS.fa ];
          then
          cp ../output_data/TRACES_viral_alignments/$ORGAN_T-$VIRUS.fa $ORGAN_T-$VIRUS.fa	
          ./TRACES_hybrid.sh $VIRUS $SCAFFOLDS_PATH $THREADS $ORGAN_T 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
	  ./TRACES_hybrid_consensus.sh $ORGAN_T-$VIRUS.fa $HYBRID_ALI_PATH/scaffolds_aligned_sorted_$VIRUS-$ORGAN_T.bam $ORGAN_T $VIRUS 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
          #
          mv $VIRUS-consensus-$ORGAN_T.fa $HYBRID_CON_PATH
	  #
          mv $VIRUS-calls-$ORGAN_T.bed $HYBRID_BED_PATH
          mv $VIRUS-coverage-$ORGAN_T.bed $HYBRID_BED_PATH
          mv $VIRUS-zero-coverage-$ORGAN_T.bed $HYBRID_BED_PATH
          #
	  mv scaffolds_aligned_sorted_$VIRUS-$ORGAN_T.bam $HYBRID_ALI_PATH
          mv scaffolds_aligned_sorted_$VIRUS-$ORGAN_T.bam.bai $HYBRID_ALI_PATH
	  mv $ORGAN_T-$VIRUS.fa $HYBRID_ALI_PATH
          mv $ORGAN_T-$VIRUS.fa.fai $HYBRID_ALI_PATH
          #
          fi
        done
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      fi
    #
    # ==========================================================================
    #
    done
    #
  # #
  #
  # ============================================================================ 
  #
  # BUILD COMPLETE VIRAL META TABLE FOR MULTIPLE ORGANS:
  #
  if [[ "$RUN_META_ON" -eq "1" ]];
    then
    ./TRACES_get_report_meta.sh
    fi
  #
  # ============================================================================
  # CLEAN DATA:
  rm -f FW_READS.fq.gz RV_READS.fq.gz
  #
  # ============================================================================
  fi
#
# ==============================================================================
################################################################################
