#!/bin/bash
#
#
declare -a READS=("healthy-skin:MT18021B_S1_L001_R1_001.fastq.gz:MT18021B_S1_L001_R2_001.fastq.gz:MT18021B_S1_L002_R1_001.fastq.gz:MT18021B_S1_L002_R2_001.fastq.gz" "skin:MT18021C_S2_L001_R1_001.fastq.gz:MT18021C_S2_L001_R2_001.fastq.gz:MT18021C_S2_L002_R1_001.fastq.gz:MT18021C_S2_L002_R2_001.fastq.gz" "bone:MT18021D_S3_L001_R1_001.fastq.gz:MT18021D_S3_L001_R2_001.fastq.gz:MT18021D_S3_L002_R1_001.fastq.gz:MT18021D_S3_L002_R2_001.fastq.gz" "colon:MT18021E_S4_L001_R1_001.fastq.gz:MT18021E_S4_L001_R2_001.fastq.gz:MT18021E_S4_L002_R1_001.fastq.gz:MT18021E_S4_L002_R2_001.fastq.gz" "heart:MT18021F_S5_L001_R1_001.fastq.gz:MT18021F_S5_L001_R2_001.fastq.gz:MT18021F_S5_L002_R1_001.fastq.gz:MT18021F_S5_L002_R2_001.fastq.gz" "liver:MT18021G_S6_L001_R1_001.fastq.gz:MT18021G_S6_L001_R2_001.fastq.gz:MT18021G_S6_L002_R1_001.fastq.gz:MT18021G_S6_L002_R2_001.fastq.gz" "spleen:MT18022B_S7_L001_R1_001.fastq.gz:MT18022B_S7_L001_R2_001.fastq.gz:MT18022B_S7_L002_R1_001.fastq.gz:MT18022B_S7_L002_R2_001.fastq.gz" "kidney:MT18022C_S8_L001_R1_001.fastq.gz:MT18022C_S8_L001_R2_001.fastq.gz:MT18022C_S8_L002_R1_001.fastq.gz:MT18022C_S8_L002_R2_001.fastq.gz" "lung:MT18022D_S9_L001_R1_001.fastq.gz:MT18022D_S9_L001_R2_001.fastq.gz:MT18022D_S9_L002_R1_001.fastq.gz:MT18022D_S9_L002_R2_001.fastq.gz" "plasma:MT18022E_S10_L001_R1_001.fastq.gz:MT18022E_S10_L001_R2_001.fastq.gz:MT18022E_S10_L002_R1_001.fastq.gz:MT18022E_S10_L002_R2_001.fastq.gz" "blood:MT18022F_S11_L001_R1_001.fastq.gz:MT18022F_S11_L001_R2_001.fastq.gz:MT18022F_S11_L002_R1_001.fastq.gz:MT18022F_S11_L002_R2_001.fastq.gz" "bone-marrow:MT18022G_S12_L001_R1_001.fastq.gz:MT18022G_S12_L001_R2_001.fastq.gz:MT18022G_S12_L002_R1_001.fastq.gz:MT18022G_S12_L002_R2_001.fastq.gz" "teeth:MT18023B_S13_L001_R1_001.fastq.gz:MT18023B_S13_L001_R2_001.fastq.gz:MT18023B_S13_L002_R1_001.fastq.gz:MT18023B_S13_L002_R2_001.fastq.gz" "brain:MT18023C_S14_L001_R1_001.fastq.gz:MT18023C_S14_L001_R2_001.fastq.gz:MT18023C_S14_L002_R1_001.fastq.gz:MT18023C_S14_L002_R2_001.fastq.gz")
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
  #gto_build_dbs.sh -vi
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
    echo "RUN ORGAN=$ORGAN_T F1=$SPL_R1A R1=$SPL_R2A F2=$SPL_R1B R2=$SPL_R2B";
    #
    ## MERGE FILES
    echo "MERGGING FILES ...";
    rm -f FW_READS.fq.gz RV_READS.fq.gz
    zcat $SPL_R1A $SPL_R1B | gzip > FW_READS.fq.gz
    zcat $SPL_R2A $SPL_R2B | gzip > RV_READS.fq.gz
    #
    ./ki_trim_filter_reads.sh
    #
    ./ki_remove_phix.sh
    #
    ./ki_metagenomics.sh $ORGAN_T
    #
    ./ki_profiles.sh GIS-$ORGAN_T VDB.fa $ORGAN_T
    #
    ./ki_extract_mito.sh
    #
    ./ki_assemble_mito.sh $ORGAN_T
    #
    done
  fi
#
# ==============================================================================
################################################################################
