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
#
INSTALL=0;
BUILD_VDB_ALL=0;
BUILD_VDB_REF=0;
BUILD_UDB=0;
#
GEN_ADAPTERS=0;
GET_PHIX=0;
GET_MITO=0;
GET_CY=0;
GET_EXTRA=0;
#
RUN_ANALYSIS=0;
#
RUN_META_ON=0;
RUN_PROFILES_ON=0;
RUN_META_NON_VIRAL_ON=0;
#
RUN_MITO_ON=0;
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
RUN_POLY1=0;
RUN_POLY2=0;
RUN_POLY3=0;
RUN_POLY4=0;
RUN_POLY5=0;
RUN_POLY6=0;
RUN_POLY7=0;
RUN_POLY8=0;
RUN_POLY9=0;
RUN_POLY10=0;
RUN_POLY11=0;
RUN_POLY12=0;
RUN_POLY13=0;
RUN_POLY14=0;
#
RUN_TTV_ON=0;
RUN_HBOV1_ON=0;
RUN_HBOVNOT1_ON=0;
RUN_HBV_ON=0;
RUN_HPV_ON=0;
RUN_VARV_ON=0;
#
RUN_SPECIFIC=0;
#
RUN_CY_ON=0;
#
RUN_DE_NOVO_ASSEMBLY=0;
#
#
# ==============================================================================
# CHECK IF FILES EXIST
#
CHECK_META_INFO () {
  if [ ! -f meta_info.txt ];
    then
    echo -e "\e[31mERROR: meta_info.txt file not found!\e[0m"
    echo "Please create a meta information file before the run."
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
CHECK_MT () {
  if [ ! -f mtDNA.fa ];
    then
    echo -e "\e[31mERROR: reference mitochondrial DNA (mtDNA.fa) not found!\e[0m"
    echo "TIP: before this, run: ./TRACESPipe.sh --get-mito"
    echo "For addition information, see the instructions at the web page."
    exit 1;
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
for i in "$@"
  do
  case $i in
    -h|--help|?)
      SHOW_HELP=1;
      shift
    ;;
    -v|-V|--version)
      SHOW_VERSION=1;
      shift
    ;;
    -i|--install)
      INSTALL=1;
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
    -gy|--get-y-chromo)
      GET_CY=1;
      SHOW_HELP=0;
      shift
    ;;
    -gx|--get-extra-vir)
      GET_EXTRA=1;
      SHOW_HELP=0;
      shift
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
      RUN_POLY1=1;
      RUN_POLY2=1;
      RUN_POLY3=1;
      RUN_POLY4=1;
      RUN_POLY5=1;
      RUN_POLY6=1;
      RUN_POLY7=1;
      RUN_POLY8=1;
      RUN_POLY9=1;
      RUN_POLY10=1;
      RUN_POLY11=1;
      RUN_POLY12=1;
      RUN_POLY13=1;
      RUN_POLY14=1;
      RUN_TTV_ON=1;
      RUN_HBOV1_ON=1;
      RUN_HBOVNOT1_ON=1;
      RUN_HBV_ON=1;
      RUN_HPV_ON=1;
      RUN_VARV_ON=1;
      RUN_CY_ON=1;
      RUN_MITO_ON=1;
      RUN_DE_NOVO_ASSEMBLY=1;
      SHOW_HELP=0;
      shift
    ;;
    -rsr|--run-specific)
      RUN_ANALYSIS=1;
      RUN_SPECIFIC=1;
      SHOW_HELP=0;
      shift
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
      RUN_POLY1=1;
      RUN_POLY2=1;
      RUN_POLY3=1;
      RUN_POLY4=1;
      RUN_POLY5=1;
      RUN_POLY6=1;
      RUN_POLY7=1;
      RUN_POLY8=1;
      RUN_POLY9=1;
      RUN_POLY10=1;
      RUN_POLY11=1;
      RUN_POLY12=1;
      RUN_POLY13=1;
      RUN_POLY14=1;
      RUN_TTV_ON=1;
      RUN_HBOV1_ON=1;
      RUN_HBOVNOT1_ON=1;
      RUN_HBV_ON=1;
      RUN_HPV_ON=1;
      RUN_VARV_ON=1;
      shift
    ;;
    -rb|--run-b19)
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
    -rcy|--run-y-chromo)
      RUN_ANALYSIS=1;
      RUN_CY_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rda|--run-de-novo)
      RUN_ANALYSIS=1;
      RUN_DE_NOVO_ASSEMBLY=1;
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
      RUN_POLY1=1;
      RUN_POLY2=1;
      RUN_POLY3=1;
      RUN_POLY4=1;
      RUN_POLY5=1;
      RUN_POLY6=1;
      RUN_POLY7=1;
      RUN_POLY8=1;
      RUN_POLY9=1;
      RUN_POLY10=1;
      RUN_POLY11=1;
      RUN_POLY12=1;
      RUN_POLY13=1;
      RUN_POLY14=1;
      RUN_TTV_ON=1;
      RUN_HBOV1_ON=1;
      RUN_HBOVNOT1_ON=1;
      RUN_HBV_ON=1;
      RUN_HPV_ON=1;
      RUN_VARV_ON=1;
      RUN_MITO_ON=1;
      RUN_CY_ON=1;
      RUN_DE_NOVO_ASSEMBLY=1;
      SHOW_HELP=0;
      shift
    ;;
    *) # unknown option
    echo "Invalid arg "$1
    echo "For help, try: ./TRACESPipe.sh -h"
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
  echo -e "    \e[32mand analysis of viral and human-host genomes at a multi-organ level\e[0m."
  echo "                                                                "
  echo -e "\e[93m    Usage: ./TRACESPipe.sh [options]                             \e[0m"
  echo "                                                                "
  echo "    -h,    --help            Show this help message and exit,     "
  echo "    -v,    --version         Show the version and some information,  "
  echo "                                                                  "
  echo "    -i,    --install         Installation of all the tools,       "
  echo "                                                                  "
  echo "    -vdb,  --build-viral     Build viral database (all sequences), "
  echo "    -vdbr, --build-viral-r   Build viral database (references only),  "
  echo "    -udb,  --build-unviral   Build non viral database (control),  "
  echo "                                                                  "
  echo "    -gx,   --get-extra-vir   Downloads/appends (VDB) extra viral seq, "
  echo "    -gad,  --gen-adapters    Generate FASTA file with adapters,   "
  echo "    -gp,   --get-phix        Extracts PhiX genomes (Needs viral DB),  "
  echo "    -gm,   --get-mito        Downloads human Mitochondrial genome,"
  echo "    -gy,   --get-y-chromo    Downloads human Y-chromosome,        "
  echo "                                                                  "
  echo "    -rm,   --run-meta        Run viral metagenomic identification,    "
  echo "    -ro,   --run-meta-nv     Run NON-viral metagenomic identification,   "
  echo "                                                                  "
  echo "    -rava, --run-all-v-alig  Run all viral align/sort/consensus seqs,    "
  echo "                                                                 "
  echo "    -rb,   --run-b19         Run B19   align and consensus seq,    "
  echo "    -rh1,  --run-hv1         Run HV1   align and consensus seq,    "
  echo "    -rh2,  --run-hv2         Run HV2   align and consensus seq,    "
  echo "    -rh3,  --run-hv3         Run HV3   align and consensus seq,    "
  echo "    -rh4,  --run-hv4         Run HV4   align and consensus seq,    "
  echo "    -rh5,  --run-hv5         Run HV5   align and consensus seq,    "
  echo "    -rh6,  --run-hv6         Run HV6   align and consensus seq,    "
  echo "    -rh6a, --run-hv6a        Run HV6A  align and consensus seq,    "
  echo "    -rh6b, --run-hv6b        Run HV6B  align and consensus seq,    "
  echo "    -rh7,  --run-hv7         Run HV7   align and consensus seq,    "
  echo "    -rh8,  --run-hv8         Run HV8   align and consensus seq,    "
  echo "    -rh8,  --run-hv8         Run HV8   align and consensus seq,    "

  echo "    -rp1,  --run-poly1       Run Polyoma 1  align and consensus seq,  "
  echo "    -rp2,  --run-poly2       Run Polyoma 2  align and consensus seq,  "
  echo "    -rp3,  --run-poly3       Run Polyoma 3  align and consensus seq,  "
  echo "    -rp4,  --run-poly4       Run Polyoma 4  align and consensus seq,  "
  echo "    -rp5,  --run-poly5       Run Polyoma 5  align and consensus seq,  "
  echo "    -rp6,  --run-poly6       Run Polyoma 6  align and consensus seq,  "
  echo "    -rp7,  --run-poly7       Run Polyoma 7  align and consensus seq,  "
  echo "    -rp8,  --run-poly8       Run Polyoma 8  align and consensus seq,  "
  echo "    -rp9,  --run-poly9       Run Polyoma 9  align and consensus seq,  "
  echo "    -rp10, --run-poly10      Run Polyoma 10 align and consensus seq,  "
  echo "    -rp11, --run-poly11      Run Polyoma 11 align and consensus seq,  "
  echo "    -rp12, --run-poly12      Run Polyoma 12 align and consensus seq,  "
  echo "    -rp13, --run-poly13      Run Polyoma 13 align and consensus seq,  "
  echo "    -rp14, --run-poly14      Run Polyoma 14 align and consensus seq,  "

  echo "    -rtt,  --run-ttv         Run TTV   align and consensus seq,    "
  echo "    -rbv1, --run-hbov1       Run HBoV1 align and consensus seq,    "
  echo "    -rbv0, --run-hbovnot1    Run HBoV (2,3,...) align/consensus seq, "
  echo "    -rhbv, --run-hbv         Run HBV   align and consensus seq,    "
  echo "    -rhpv, --run-hpv         Run HPV   align and consensus seq,    "
  echo "    -rvar, --run-varv        Run VARV  align and consensus seq,    "
  echo "                                                                 "
  echo "    -rsr,  --run-specific    Run specific REF align/consensus seq, "
  echo "                                                                 "
  echo "    -rmt,  --run-mito        Run Mito  align and consensus seq,   "
  echo "    -rcy,  --run-y-chromo    Run CY    align and consensus seq,    "
  echo "                                                                  "
  echo "    -rda,  --run-de-novo     Run de-novo assembly,               "
  echo "                                                                 "
  echo "    -ra,   --run-analysis    Run data analysis,                   "
  echo "    -all,  --run-all         Run all the options.                 "
  echo "                                                                "
  echo -e "\e[93m    Example: ./TRACESPipe.sh --run-meta --run-b19 --run-mito \e[0m"
  echo "                                                                "
  echo "    meta_info.txt -> 'organ:reads_forward.fa.gz:reads_reverse.fa.gz'  "
  echo "    The reads and meta_info.txt must be in the src/ folder.     "
  echo "                                                                "
  exit 1
  fi
#
# ==============================================================================
# VERSION
#
if [ "$SHOW_VERSION" -eq "1" ];
  then
  echo "                                       ";
  echo "           TRACESPipe                  ";
  echo "                                       ";
  echo "         Version: 1.0.0                ";
  echo "                                       ";
  echo "   Department of Virology and          ";
  echo " Department of Forensic Medicine,      ";
  echo " University of Helsinki, Finland.      ";
  echo "                &                      ";
  echo "           IEETA/DETI,                 ";
  echo " University of Aveiro, Portugal.       ";
  echo "                                       ";
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
if [[ "$GET_CY" -eq "1" ]];
  then
  ./TRACES_get_cy.sh
  fi
#
# ==============================================================================
#
if [[ "$GET_EXTRA" -eq "1" ]];
  then
  ./TRACES_get_enriched_sequences.sh VDB.fa
  fi
#
# ==============================================================================
#
if [[ "$RUN_ANALYSIS" -eq "1" ]];
  then
  #
  CHECK_META_INFO;
  #
  mapfile -t READS < meta_info.txt
  #
  for read in "${READS[@]}" # 
    do
    #
    ORGAN_T=`echo $read | tr ':' '\t' | awk '{ print $1 }'`;
    SPL_Forward=`echo $read | tr ':' '\t' | awk '{ print $2 }'`;
    SPL_Reverse=`echo $read | tr ':' '\t' | awk '{ print $3 }'`;
    echo -e "\e[34m[TRACES]\e[93m Running: Organ=$ORGAN_T Forward=$SPL_Forward Reverse=$SPL_Reverse\e[0m";
    #
    rm -f FW_READS.fq.gz RV_READS.fq.gz
    echo -e "\e[34m[TRACES]\e[32m Copping an instance of the files ...\e[0m";
    cp $SPL_Forward FW_READS.fq.gz;
    cp $SPL_Reverse RV_READS.fq.gz;
    echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
    #
    # ==========================================================================
    # TRIM AND FILTER READS
    #
    CHECK_ADAPTERS;
    #
    echo -e "\e[34m[TRACES]\e[32m Trimming and filtering with Trimmomatic ...\e[0m";
    ./TRACES_trim_filter_reads.sh
    echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
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
      echo -e "\e[34m[TRACES]\e[32m Removing PhiX from the samples with MAGNET ...\e[0m";
      ./TRACES_remove_phix.sh # IT IS USED ONLY FOR FALCON
      #
      # fastq_pair test_R1.fastq test_R2.fastq: [needs adaptation]
      # IF you want to remove Phix also before assembly
      #
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[TRACES]\e[32m Running viral metagenomic analysis with FALCON ...\e[0m";
      ./TRACES_metagenomics_viral.sh $ORGAN_T VDB.fa 10000
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[TRACES]\e[32m Finding the best references ...\e[0m";
      ./TRACES_get_best_B19.sh $ORGAN_T > REPORT_META_VIRAL_$ORGAN_T.txt
      #
      ./TRACES_get_best_HV1.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_HV2.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_HV3.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_HV4.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_HV5.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_HV6.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_HV6A.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_HV6B.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_HV7.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_HV8.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      #
      ./TRACES_get_best_HBV.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      #
      ./TRACES_get_best_Poly1.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_Poly2.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_Poly3.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_Poly4.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_Poly5.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_Poly6.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_Poly7.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_Poly8.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_Poly9.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_Poly10.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_Poly11.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_Poly12.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_Poly13.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_Poly14.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      # 
      ./TRACES_get_best_HPV.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_TTV.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_HBoV1.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_HBoVnot1.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_VARV.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
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
      echo -e "\e[34m[TRACES]\e[32m Building complexity profiles with gto ...\e[0m";
      cat NP-o_fw_pr.fq NP-o_fw_unpr.fq NP-o_rv_pr.fq NP-o_rv_unpr.fq > P_TRACES_sample_reads.fq
      ./TRACES_profiles.sh GIS-$ORGAN_T VDB.fa P_TRACES_sample_reads.fq $ORGAN_T
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
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
      echo -e "\e[34m[TRACES]\e[32m Running NON viral metagenomic analysis with FALCON ...\e[0m";
      ./TRACES_metagenomics.sh $ORGAN_T DB.fa 10000 
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
      fi
    #
    # ==========================================================================
    # RUN SPECIFIC SPECIFIC ALIGN/CONSENSUS
    #
    if [[ "$RUN_SPECIFIC" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to B19 ref with bowtie2 ...\e[0m";
      #
      echo "Extracting sequence from VDB.fa ..."
      #
      CHECK_VDB;
      #
      ###gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > SPECIFIC-$IDN.fa
      echo "Aliggning ..."
      ###./TRACES_b19_align_reads.sh SPECIFIC-$IDN.fa $ORGAN_T
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
      ###./TRACES_b19_consensus.sh SPECIFIC-$IDN.fa b19_aligned_sorted-$ORGAN_T.bam $ORGAN_T $IDN
      fi
    #
    # ========================================================================== 
    # DETAILED VIRAL ALIGN/CONSENSUS
    #
    if [[ "$RUN_B19_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to B19 ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_B19.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]];
        then
        echo "Extracting sequence from VDB.fa ..."
	#
        CHECK_VDB;
        #
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-B19.fa
        echo "Aliggning ..."
	./TRACES_viral_align_reads.sh $ORGAN_T-B19.fa $ORGAN_T B19
        #./TRACES_b19_align_reads.sh $ORGAN_T-B19.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_viral_consensus.sh $ORGAN_T-B19.fa viral_aligned_sorted-$ORGAN_T-B19.bam $ORGAN_T B19
        fi
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
      fi
    #
    if [[ "$RUN_HV1_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to HV1 ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_HV1.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]];
        then
        echo "Extracting sequence from VDB.fa ..."
	#
	CHECK_VDB;
	#
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-HV1.fa
        echo "Aliggning ..."
        ./TRACES_hv1_align_reads.sh $ORGAN_T-HV1.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_viral_consensus.sh $ORGAN_T-HV1.fa hv1_aligned_sorted-$ORGAN_T.bam $ORGAN_T HV1
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
        fi
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
      fi
    #
    if [[ "$RUN_HV2_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to HV2 ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_HV2.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]];
        then
        echo "Extracting sequence from VDB.fa ..."
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-HV2.fa
        echo "Aliggning ..."
        ./TRACES_hv2_align_reads.sh $ORGAN_T-HV2.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_viral_consensus.sh $ORGAN_T-HV2.fa hv1_aligned_sorted-$ORGAN_T.bam $ORGAN_T HV2
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
        fi
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
      fi
    #
    if [[ "$RUN_HV3_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to HV3 ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_HV3.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]];
        then
        echo "Extracting sequence from VDB.fa ..."
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-HV3.fa
        echo "Aliggning ..."
        ./TRACES_hv3_align_reads.sh $ORGAN_T-HV3.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_viral_consensus.sh $ORGAN_T-HV3.fa hv3_aligned_sorted-$ORGAN_T.bam $ORGAN_T HV3
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
        fi
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
      fi
    #
    if [[ "$RUN_HV4_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to HV4 ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_HV4.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]];
        then
        echo "Extracting sequence from VDB.fa ..."
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-HV4.fa
        echo "Aliggning ..."
        ./TRACES_hv4_align_reads.sh $ORGAN_T-HV4.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_viral_consensus.sh $ORGAN_T-HV4.fa hv4_aligned_sorted-$ORGAN_T.bam $ORGAN_T HV4
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
        fi
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
      fi
    #
    if [[ "$RUN_HV5_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to HV5 ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_HV5.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]];
        then
        echo "Extracting sequence from VDB.fa ..."
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-HV5.fa
        echo "Aliggning ..."
        ./TRACES_hv5_align_reads.sh $ORGAN_T-HV5.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_viral_consensus.sh $ORGAN_T-HV5.fa hv5_aligned_sorted-$ORGAN_T.bam $ORGAN_T HV5
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
        fi
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
      fi
    #
    if [[ "$RUN_HV6_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to HV6 ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_HV6.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]];
        then
        echo "Extracting sequence from VDB.fa ..."
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-HV6.fa
        echo "Aliggning ..."
        ./TRACES_hv6_align_reads.sh $ORGAN_T-HV6.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_viral_consensus.sh $ORGAN_T-HV6.fa hv6_aligned_sorted-$ORGAN_T.bam $ORGAN_T HV6
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
        fi
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
      fi
    #
    if [[ "$RUN_HV6A_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to HV6A ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_HV6A.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]];
        then
        echo "Extracting sequence from VDB.fa ..."
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-HV6A.fa
        echo "Aliggning ..."
        ./TRACES_hv6A_align_reads.sh $ORGAN_T-HV6A.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_viral_consensus.sh $ORGAN_T-HV6A.fa hv6a_aligned_sorted-$ORGAN_T.bam $ORGAN_T HV6A
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
        fi
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
      fi
    #
    if [[ "$RUN_HV6B_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to HV6B ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_HV6B.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]];
        then
        echo "Extracting sequence from VDB.fa ..."
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-HV6B.fa
        echo "Aliggning ..."
        ./TRACES_hv6B_align_reads.sh $ORGAN_T-HV6B.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_viral_consensus.sh $ORGAN_T-HV6B.fa hv6b_aligned_sorted-$ORGAN_T.bam $ORGAN_T HV6B
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
        fi
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
      fi
    #
    if [[ "$RUN_HV7_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to HV7 ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_HV7.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]];
        then
        echo "Extracting sequence from VDB.fa ..."
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-HV7.fa
        echo "Aliggning ..."
        ./TRACES_hv7_align_reads.sh $ORGAN_T-HV7.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_viral_consensus.sh $ORGAN_T-HV7.fa hv7_aligned_sorted-$ORGAN_T.bam $ORGAN_T HV7
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
        fi
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
      fi
    #
    if [[ "$RUN_HV8_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to HV8 ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_HV8.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]];
        then
        echo "Extracting sequence from VDB.fa ..."
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-HV8.fa
        echo "Aliggning ..."
        ./TRACES_hv8_align_reads.sh $ORGAN_T-HV8.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_viral_consensus.sh $ORGAN_T-HV8.fa hv8_aligned_sorted-$ORGAN_T.bam $ORGAN_T HV8
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
        fi
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
      fi
    #

  ############################################################################################

    if [[ "$RUN_POLY1_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to Polyomavirus 1 ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_Poly1.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]]; 
        then 
        echo "Extracting sequence from VDB.fa ...";
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-POLY1.fa
        echo "Aliggning ..." 
        ./TRACES_poly1_align_reads.sh $ORGAN_T-POLY1.fa $ORGAN_T 12
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"; 
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m"; 
        ./TRACES_viral_consensus.sh $ORGAN_T-POLY1.fa poly1_aligned_sorted-$ORGAN_T.bam $ORGAN_T POLY1
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m" 
        fi 
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
      fi 
    #
    if [[ "$RUN_POLY2_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to Polyomavirus 2 ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_Poly2.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]]; 
        then 
        echo "Extracting sequence from VDB.fa ...";
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-POLY2.fa
        echo "Aliggning ..." 
        ./TRACES_poly2_align_reads.sh $ORGAN_T-POLY2.fa $ORGAN_T 12
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"; 
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m"; 
        ./TRACES_viral_consensus.sh $ORGAN_T-POLY2.fa poly2_aligned_sorted-$ORGAN_T.bam $ORGAN_T POLY2
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m" 
        fi 
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
      fi 
    #
    if [[ "$RUN_POLY3_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to Polyomavirus 3 ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_Poly3.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]]; 
        then 
        echo "Extracting sequence from VDB.fa ...";
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-POLY3.fa
        echo "Aliggning ..." 
        ./TRACES_poly3_align_reads.sh $ORGAN_T-POLY3.fa $ORGAN_T 12
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"; 
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m"; 
        ./TRACES_viral_consensus.sh $ORGAN_T-POLY3.fa poly3_aligned_sorted-$ORGAN_T.bam $ORGAN_T POLY3
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m" 
        fi 
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
      fi 
    #
    if [[ "$RUN_POLY4_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to Polyomavirus 4 ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_Poly4.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]]; 
        then 
        echo "Extracting sequence from VDB.fa ...";
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-POLY4.fa
        echo "Aliggning ..." 
        ./TRACES_poly4_align_reads.sh $ORGAN_T-POLY4.fa $ORGAN_T 12
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"; 
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m"; 
        ./TRACES_viral_consensus.sh $ORGAN_T-POLY4.fa poly4_aligned_sorted-$ORGAN_T.bam $ORGAN_T POLY4
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m" 
        fi 
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
      fi 
    #
    if [[ "$RUN_POLY5_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to Polyomavirus 5 ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_Poly5.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]]; 
        then 
        echo "Extracting sequence from VDB.fa ...";
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-POLY5.fa
        echo "Aliggning ..." 
        ./TRACES_poly5_align_reads.sh $ORGAN_T-POLY5.fa $ORGAN_T 12
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"; 
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m"; 
        ./TRACES_viral_consensus.sh $ORGAN_T-POLY5.fa poly5_aligned_sorted-$ORGAN_T.bam $ORGAN_T POLY5
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m" 
        fi 
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
      fi 
    #
    if [[ "$RUN_POLY6_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to Polyomavirus 6 ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_Poly6.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]]; 
        then 
        echo "Extracting sequence from VDB.fa ...";
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-POLY6.fa
        echo "Aliggning ..." 
        ./TRACES_poly6_align_reads.sh $ORGAN_T-POLY6.fa $ORGAN_T 12
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"; 
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m"; 
        ./TRACES_viral_consensus.sh $ORGAN_T-POLY6.fa poly6_aligned_sorted-$ORGAN_T.bam $ORGAN_T POLY6
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m" 
        fi 
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
      fi 
    #
    if [[ "$RUN_POLY7_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to Polyomavirus 7 ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_Poly7.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]]; 
        then 
        echo "Extracting sequence from VDB.fa ...";
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-POLY7.fa
        echo "Aliggning ..." 
        ./TRACES_poly7_align_reads.sh $ORGAN_T-POLY7.fa $ORGAN_T 12
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"; 
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m"; 
        ./TRACES_viral_consensus.sh $ORGAN_T-POLY7.fa poly7_aligned_sorted-$ORGAN_T.bam $ORGAN_T POLY7
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m" 
        fi 
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
      fi 
    #
    if [[ "$RUN_POLY8_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to Polyomavirus 8 ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_Poly8.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]]; 
        then 
        echo "Extracting sequence from VDB.fa ...";
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-POLY8.fa
        echo "Aliggning ..." 
        ./TRACES_poly8_align_reads.sh $ORGAN_T-POLY8.fa $ORGAN_T 12
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"; 
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m"; 
        ./TRACES_viral_consensus.sh $ORGAN_T-POLY8.fa poly8_aligned_sorted-$ORGAN_T.bam $ORGAN_T POLY8
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m" 
        fi 
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
      fi 
    #
    if [[ "$RUN_POLY9_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to Polyomavirus 9 ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_Poly9.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]]; 
        then 
        echo "Extracting sequence from VDB.fa ...";
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-POLY9.fa
        echo "Aliggning ..." 
        ./TRACES_poly9_align_reads.sh $ORGAN_T-POLY9.fa $ORGAN_T 12
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"; 
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m"; 
        ./TRACES_viral_consensus.sh $ORGAN_T-POLY9.fa poly9_aligned_sorted-$ORGAN_T.bam $ORGAN_T POLY9
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m" 
        fi 
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
      fi 
    #
    if [[ "$RUN_POLY10_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to Polyomavirus 10 ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_Poly10.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]]; 
        then 
        echo "Extracting sequence from VDB.fa ...";
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-POLY10.fa
        echo "Aliggning ..." 
        ./TRACES_poly10_align_reads.sh $ORGAN_T-POLY10.fa $ORGAN_T 12
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"; 
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m"; 
        ./TRACES_viral_consensus.sh $ORGAN_T-POLY10.fa poly10_aligned_sorted-$ORGAN_T.bam $ORGAN_T POLY10
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m" 
        fi 
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
      fi 
    #
    if [[ "$RUN_POLY11_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to Polyomavirus 11 ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_Poly11.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]]; 
        then 
        echo "Extracting sequence from VDB.fa ...";
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-POLY11.fa
        echo "Aliggning ..." 
        ./TRACES_poly11_align_reads.sh $ORGAN_T-POLY11.fa $ORGAN_T 12
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"; 
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m"; 
        ./TRACES_viral_consensus.sh $ORGAN_T-POLY11.fa poly11_aligned_sorted-$ORGAN_T.bam $ORGAN_T POLY11
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m" 
        fi 
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
      fi 
    #
    if [[ "$RUN_POLY12_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to Polyomavirus 12 ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_Poly12.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]]; 
        then 
        echo "Extracting sequence from VDB.fa ...";
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-POLY12.fa
        echo "Aliggning ..." 
        ./TRACES_poly12_align_reads.sh $ORGAN_T-POLY12.fa $ORGAN_T 12
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"; 
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m"; 
        ./TRACES_viral_consensus.sh $ORGAN_T-POLY12.fa poly12_aligned_sorted-$ORGAN_T.bam $ORGAN_T POLY12
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m" 
        fi 
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
      fi 
    #
    if [[ "$RUN_POLY13_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to Polyomavirus 13 ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_Poly13.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]]; 
        then 
        echo "Extracting sequence from VDB.fa ...";
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-POLY13.fa
        echo "Aliggning ..." 
        ./TRACES_poly13_align_reads.sh $ORGAN_T-POLY13.fa $ORGAN_T 12
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"; 
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m"; 
        ./TRACES_viral_consensus.sh $ORGAN_T-POLY13.fa poly13_aligned_sorted-$ORGAN_T.bam $ORGAN_T POLY13
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m" 
        fi 
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
      fi 
    #
    if [[ "$RUN_POLY14_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to Polyomavirus 14 ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_Poly14.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]]; 
        then 
        echo "Extracting sequence from VDB.fa ...";
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-POLY14.fa
        echo "Aliggning ..." 
        ./TRACES_poly14_align_reads.sh $ORGAN_T-POLY14.fa $ORGAN_T 12
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"; 
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m"; 
        ./TRACES_viral_consensus.sh $ORGAN_T-POLY14.fa poly14_aligned_sorted-$ORGAN_T.bam $ORGAN_T POLY14
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m" 
        fi 
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
      fi 
    #

  ############################################################################################

    #
    if [[ "$RUN_TTV_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to TTV ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_TTV.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]];
        then
        echo "Extracting sequence from VDB.fa ..."
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-TTV.fa
        echo "Aliggning ..."
        ./TRACES_ttv_align_reads.sh $ORGAN_T-TTV.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_viral_consensus.sh $ORGAN_T-TTV.fa ttv_aligned_sorted-$ORGAN_T.bam $ORGAN_T TTV
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
        fi
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
      fi
    #
    if [[ "$RUN_HBOV1_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to HBoV1 ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_HBoV1.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]];
        then
        echo "Extracting sequence from HBoV1.fa ..."
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-HBoV1.fa
        echo "Aliggning ..."
        ./TRACES_hbov1_align_reads.sh $ORGAN_T-HBoV1.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_viral_consensus.sh $ORGAN_T-HBoV1.fa hbov1_aligned_sorted-$ORGAN_T.bam $ORGAN_T HBoV1
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
        fi
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
      fi
    #
    if [[ "$RUN_HBOVNOT1_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to HBoV NOT 1 ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_HBoVnot1.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]];
        then
        echo "Extracting sequence from HBoVnot1.fa ..."
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-HBoVnot1.fa
        echo "Aliggning ..."
        ./TRACES_hbovnot1_align_reads.sh $ORGAN_T-HBoVnot1.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_viral_consensus.sh $ORGAN_T-HBoVnot1.fa hbovnot1_aligned_sorted-$ORGAN_T.bam $ORGAN_T HBoVnot1
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
        fi
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
      fi
    #
    if [[ "$RUN_HBV_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to HBV ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_HBV.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]];
        then
        echo "Extracting sequence from HBV.fa ..."
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-HBV.fa
        echo "Aliggning ..."
        ./TRACES_hbv_align_reads.sh $ORGAN_T-HBV.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_viral_consensus.sh $ORGAN_T-HBV.fa hbv_aligned_sorted-$ORGAN_T.bam $ORGAN_T HBV
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
        fi
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
      fi
    #
    if [[ "$RUN_HPV_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to HPV ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_HPV.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]];
        then
        echo "Extracting sequence from HPV.fa ..."
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-HPV.fa
        echo "Aliggning ..."
        ./TRACES_hpv_align_reads.sh $ORGAN_T-HPV.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_viral_consensus.sh $ORGAN_T-HPV.fa hpv_aligned_sorted-$ORGAN_T.bam $ORGAN_T HPV
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
        fi
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
      fi
    #
    if [[ "$RUN_VARV_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to VARV ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_VARV.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]];
        then
        echo "Extracting sequence from VARV.fa ..."
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-VARV.fa
        echo "Aliggning ..."
        ./TRACES_varv_align_reads.sh $ORGAN_T-VARV.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_viral_consensus.sh $ORGAN_T-VARV.fa varv_aligned_sorted-$ORGAN_T.bam $ORGAN_T VARV
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
        fi
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
      fi
    #
    # ==========================================================================
    # MITOCHONDRIAL GENOME ALIGN AND CONSENSUS
    #
    if [[ "$RUN_MITO_ON" -eq "1" ]];
      then
      #
      CHECK_MT;
      #
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to mitochondrial ref with bowtie2 ...\e[0m";
      ./TRACES_mt_align_reads.sh mtDNA.fa $ORGAN_T
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
      ./TRACES_mt_consensus.sh mtDNA.fa mt_aligned_sorted-$ORGAN_T.bam $ORGAN_T
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
      fi
    #
    # ========================================================================== 
    # CY VERYFICATION
    #
    if [[ "$RUN_CY_ON" -eq "1" ]];
      then
      #
      CHECK_CY;
      #
      echo -e "\e[34m[TRACES]\e[32m Searching for Y chromosome halotypes ...\e[0m";
      ./TRACES_cy_markers.sh $ORGAN_T
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
      fi
    #
    # ==========================================================================
    # DE-NOVO ASSEMBLY
    #
    if [[ "$RUN_DE_NOVO_ASSEMBLY" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Running do-novo DNA assembly with SPAdes ...\e[0m";
      ./TRACES_assemble_all.sh $ORGAN_T
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
      fi
    #
    # ==========================================================================
    done
  #
  # ============================================================================
  # REDIRECT RESULTS TO SPECIFIC FOLDERS
  # ============================================================================
  #
  # RESULTS WITH REPORTS AND IMAGES
  mkdir -p TRACES_results;
  rm -f TRACES_results/*
  mv *.pdf TRACES_results/
  mv *.svg TRACES_results/
  mv REPORT_META_VIRAL_*.txt TRACES_results/
  #
  # ============================================================================ 
  # MITOCHONDRIAL SEQUENCES
  #
  # CONSENSUS MITO FILES
  mkdir -p TRACES_mtdna_consensus;
  #rm -f TRACES_consensus/*
  mv mt-consensus-*.fa TRACES_mtdna_consensus/
  #
  # ALIGNMENT MITO FILES
  mkdir -p TRACES_mtdna_alignments;
  #rm -f TRACES_mtdna_alignments/*
  cp mtDNA.fa TRACES_mtdna_alignments/
  cp mtDNA.fa.fai TRACES_mtdna_alignments/
  mv mt_aligned_sorted-*.bam TRACES_mtdna_alignments/
  mv mt_aligned_sorted-*.bam.bai TRACES_mtdna_alignments/
  #
  # BED MITO FILES
  mkdir -p TRACES_mtdna_bed;
  #rm -f TRACES_mtdna_bed/*
  mv mt-calls-*.bed TRACES_mdna_bed/
  #
  # ============================================================================
  # VIRAL SEQUENCES 
  #
  # CONSENSUS VIRAL FILES
  mkdir -p TRACES_viral_consensus;
  #rm -f TRACES_viral_consensus/*
  mv *-consensus-*.fa TRACES_viral_consensus/
  #
  # ALIGNMENT VIRAL FILES
  mkdir -p TRACES_viral_alignments;
  #rm -f TRACES_viral_alignments/*
  cp *.fa TRACES_viral_alignments/
  cp *.fa.fai TRACES_viral_alignments/
  rm -f TRACES_viral_alignments/VDB.fa
  mv *_aligned_sorted-*.bam TRACES_viral_alignments/
  mv *_aligned_sorted-*.bam.bai TRACES_viral_alignments/
  #
  # BED VIRAL FILES
  mkdir -p TRACES_viral_bed;
  #rm -f TRACES_viral_bed/*
  mv *-calls-*.bed TRACES_viral_bed/
  #
  # BUILD COMPLETE VIRAL META TABLE FOR MULTIPLE ORGANS:
  ./TRACES_get_report_meta.sh
  #
  fi
#
# ==============================================================================
################################################################################
