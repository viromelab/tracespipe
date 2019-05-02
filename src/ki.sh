#!/bin/bash
#
# 
#
#
mapfile -t READS < reads_info.txt
#
SHOW_HELP=0;
INSTALL=0;
BUILD_VDB=0;
GEN_ADAPTERS=0;
GET_PHIX=0;
GET_MITO=0;
RUN_ANALYSIS=0;
#
for i in "$@"
  do
  case $i in
    -h|--help)
      SHOW_HELP=1;
      shift
    ;;
    -i|--install)
      INSTALL=1;
      SHOW_HELP=0;
      shift
    ;;
    -vi|--build-viral)
      BUILD_VDB=1;
      SHOW_HELP=0;
      shift
    ;;
    -gad|--gen-adapters)
      GEN_ADAPTERS=1;
      SHOW_HELP=0;
      shift
    ;;
    -gp|--get-phix)
      GET_PHIX=1;
      SHOW_HELP=0;
      shift
    ;;
    -gm|--get-mito)
      GET_MITO=1;
      SHOW_HELP=0;
      shift
    ;;
    -ra|--run-analysis)
      RUN_ANALYSIS=1;
      SHOW_HELP=0;
      shift
    ;;
    -all|--run-all)
      INSTALL=1;
      BUILD_VDB=1;
      GEN_ADAPTERS=1;
      GET_PHIX=1;
      GET_MITO=1;
      RUN_ANALYSIS=1;
      SHOW_HELP=0;
      shift
    ;;
    *) # unknown option
    echo "Invalid arg "$1
    echo "For help, try: ./ki.sh -h"
    ;;
  esac
  done
#
# ==============================================================================
# HELP
#
if [ "$SHOW_HELP" -eq "1" ];
  then
    echo "                                                                "
    echo "                                                                "
    echo -e "\e[34m                       ██╗  ██╗ ██╗                             "
    echo "                       ██║ ██╔╝ ██║                             "
    echo "                       █████╔╝  ██║                             "
    echo "                       ██╔═██╗  ██║                             "
    echo "                       ██║  ██╗ ██║                             "
    echo -e "                       ╚═╝  ╚═╝ ╚═╝                             \e[0m"
    echo "                                                                "
    echo -e "\e[93m    Usage: ./ki.sh [options]                                    \e[0m"
    echo "                                                                "
    echo -e "    \e[32mAn automatic pipeline for viral genome identification\e[0m" 
    echo -e "    \e[32min the contexts of clinical virology and forensics\e[0m.         "
    echo "                                                                "
    echo "    -h,   --help           Show this help message and exit,     "
    echo "    -i,   --install        Installation of all the tools,       "
    echo "    -vdb, --build-viral    Build viral database,                "
    echo "    -gad, --gen-adapters   Generate FASTA file with adapters,   "
    echo "    -gp,  --get-phix       Downloads PhiX genomes,              "
    echo "    -gm,  --get-mito       Downloads human Mitochondrial genome,"
    echo "    -ra,  --run-analysis   Run data analysis,                   "
    echo "                                                                "
    echo "    -all, --run-all        Run all the options.                 "
    echo "                                                                "
    echo -e "\e[93m    Example: ./ki.sh -all                                         \e[0m"
    echo "                                                                "
    echo "    reads_info.txt -> 'name:readsf1:readsf1:readsr1:readsr2'    "
    echo "    The reads must be in the src/ folder.                       "
    echo "                                                                "
    exit 1
  fi
#
# ==============================================================================
#
if [[ "$INSTALL" -eq "1" ]];
  then
  ./ki_install.sh
  fi
#
# ==============================================================================
#
if [[ "$BUILD_VDB" -eq "1" ]];
  then
  #gto_build_dbs -vi
  gunzip VDB.fa.gz
  fi
#
# ==============================================================================
#
if [[ "$GEN_ADAPTERS" -eq "1" ]];
  then
  ./ki_generate_adapters.sh
  fi
#
# ==============================================================================
#
if [[ "$GET_PHIX" -eq "1" ]];
  then
  ./ki_get_phix.sh
  fi
#
# ==============================================================================
#
if [[ "$GET_MITO" -eq "1" ]];
  then
  ./ki_get_mito.sh
  fi
#
# ==============================================================================
#
if [[ "$RUN_ANALYSIS" -eq "1" ]];
  then
  for read in "${READS[@]}" # 
    do
    ORGAN_T=`echo $read | tr ':' '\t' | awk '{ print $1 }'`;
    SPL_R1A=`echo $read | tr ':' '\t' | awk '{ print $2 }'`;
    SPL_R2A=`echo $read | tr ':' '\t' | awk '{ print $3 }'`;
    SPL_R1B=`echo $read | tr ':' '\t' | awk '{ print $4 }'`;
    SPL_R2B=`echo $read | tr ':' '\t' | awk '{ print $5 }'`;
    echo -e "\e[34m[ki]\e[93m Running ORGAN=$ORGAN_T F1=$SPL_R1A R1=$SPL_R2A F2=$SPL_R1B R2=$SPL_R2B \e[0m";
    #
#   # MAKE RESULTS FOLDER & CLEAN
#   mkdir -p results;
#   rm -f results/*
    #
    # MERGE FILES
    rm -f FW_READS.fq.gz RV_READS.fq.gz
    echo -e "\e[34m[ki]\e[32m Mergging the files ...\e[0m";
    zcat $SPL_R1A $SPL_R1B | gzip > FW_READS.fq.gz
    zcat $SPL_R2A $SPL_R2B | gzip > RV_READS.fq.gz
    echo -e "\e[34m[ki]\e[32m Done!\e[0m";
    #
    echo -e "\e[34m[ki]\e[32m Trimming and filtering with Trimmomatic ...\e[0m";
    ./ki_trim_filter_reads.sh
    echo -e "\e[34m[ki]\e[32m Done!\e[0m";
    #
    echo -e "\e[34m[ki]\e[32m Removing PhiX from the samples with MAGNET ...\e[0m";
    ./ki_remove_phix.sh
    echo -e "\e[34m[ki]\e[32m Done!\e[0m";
    #
    echo -e "\e[34m[ki]\e[32m Running metagenomic analysis with FALCON ...\e[0m";
    ./ki_metagenomics.sh $ORGAN_T
    echo -e "\e[34m[ki]\e[32m Done!\e[0m";
    #
    echo -e "\e[34m[ki]\e[32m Building complexity profiles with gto ...\e[0m";
    cat NP-o_fw_pr.fq NP-o_fw_unpr.fq NP-o_rv_pr.fq NP-o_rv_unpr.fq > ki_sample_reads.fq
    ./ki_profiles.sh GIS-$ORGAN_T VDB.fa ki_sample_reads.fq $ORGAN_T
    echo -e "\e[34m[ki]\e[32m Done!\e[0m";
    #
    echo -e "\e[34m[ki]\e[32m Extracting mitochondrial reads with MAGNET ...\e[0m";
    ./ki_extract_mito.sh
    echo -e "\e[34m[ki]\e[32m Done!\e[0m";
    #
    echo -e "\e[34m[ki]\e[32m Running mitochondrial DNA assembly with SPAdes ...\e[0m";
    ./ki_assemble_mito.sh $ORGAN_T
    echo -e "\e[34m[ki]\e[32m Done!\e[0m";
    #
    done
  fi
#
# ==============================================================================
################################################################################
