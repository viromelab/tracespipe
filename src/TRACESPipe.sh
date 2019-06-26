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
GET_THREADS=0;
THREADS=0;
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
#
RUN_DECRYPT=0;
RUN_ENCRYPT=0;
#
RUN_SPECIFIC=0;
#
RUN_CY_ON=0;
RUN_CY_QUANT_ON=0;
#
RUN_DE_NOVO_ASSEMBLY=0;
#
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
# ==============================================================================
#
ALIGN_AND_CONSENSUS () {
  #
  V_TAG="$1";
  echo -e "\e[34m[TRACESPipe]\e[32m Aliggning reads to $V_TAG best reference with bowtie2 ...\e[0m";
  CHECK_TOP "$ORGAN_T";
  V_INFO=`./TRACES_get_best_$V_TAG.sh $ORGAN_T`;
  echo "Best match: $V_INFO";
  V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
  if [[ "$V_GID" != "-" ]];
    then
    echo "Extracting sequence from VDB.fa ..."
    CHECK_VDB;
    gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > $ORGAN_T-$V_TAG.fa
    echo "Aliggning ..."
    ./TRACES_viral_align_reads.sh $ORGAN_T-$V_TAG.fa $ORGAN_T $V_TAG $THREADS
    echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
    #
    echo -e "\e[34m[TRACESPipe]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
    ./TRACES_viral_consensus.sh $ORGAN_T-$V_TAG.fa viral_aligned_sorted-$ORGAN_T-$V_TAG.bam $ORGAN_T $V_TAG
    mkdir -p ../output_data/TRACES_viral_alignments;
    #rm -f ../output_data/TRACES_viral_alignments/*
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
    fi
  echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
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
      RUN_CY_ON=1;
      RUN_CY_QUANT_ON=1;
      RUN_MITO_ON=1;
      RUN_DE_NOVO_ASSEMBLY=1;
      SHOW_HELP=0;
      shift
    ;;
    -rsr|--run-specific)
      RUN_ANALYSIS=1;
      RUN_SPECIFIC=1;
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
    -rya|--run-y-align)
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
    -rda|--run-de-novo)
      RUN_ANALYSIS=1;
      RUN_DE_NOVO_ASSEMBLY=1;
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
      RUN_MITO_ON=1;
      RUN_CY_ON=1;
      RUN_CY_QUANT_ON=1;
      RUN_DE_NOVO_ASSEMBLY=1;
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
  echo -e "    \e[32mand analysis of viral and human-host genomes at a multi-organ level\e[0m."
  echo "                                                                "
  echo -e "\e[93m    Usage: ./TRACESPipe.sh [options]                             \e[0m"
  echo "                                                                "
  echo "    -h,    --help            Show this help message and exit,     "
  echo "    -v,    --version         Show the version and some information,  "
  echo "    -f,    --force           Force running and overwrite of files,  "
  echo "                                                                  "
  echo "    -gmt,  --get-max-threads Get the number of maximum machine threads, "
  echo "    -t <THREADS>, --threads <THREADS>                             "
  echo "                             Number of threads to use, "
  echo "                                                                  "
  echo "    -i,    --install         Installation of all the tools,       "
  echo "                                                                  "
  echo "    -dec,  --decrypt         Decrypt (all files in ../encrypted_data),  "
  echo "    -enc,  --encrypt         Encrypt (all files in ../to_encrypt_data),  "
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
  echo "    -gax,  --get-all-aux     Runs -gad -gp -gm -gy,               "
  echo "                                                                  "
  echo "    -rm,   --run-meta        Run viral metagenomic identification,    "
  echo "    -ro,   --run-meta-nv     Run NON-viral metagenomic identification,   "
  echo "                                                                  "
  echo "    -rava, --run-all-v-alig  Run all viral align/sort/consensus seqs,    "
  echo "                                                                 "
  echo "    -rb19, --run-b19         Run B19   align and consensus seq,    "
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
  echo "    -rsr <ID>, --run-specific <ID/PATTERN>                        "
  echo "                             Run specific reference align/consensus, "
  echo "                                                                 "
  echo "    -rmt,  --run-mito        Run Mito align and consensus seq,   "
  echo "                                                                 "
  echo "    -rya,  --run-cy-align    Run CY align and consensus seq,    "
  echo "    -ryq,  --run-cy-quant    Estimate the quantity of CY DNA,    "
  echo "                                                                  "
  echo "    -rda,  --run-de-novo     Run de-novo assembly,               "
  echo "                                                                 "
  echo "    -ra,   --run-analysis    Run data analysis,                   "
  echo "    -all,  --run-all         Run all the options.                 "
  echo "                                                                "
  echo -e "\e[93m    Example: ./TRACESPipe.sh --run-meta --run-b19 --run-mito \e[0m"
  echo "                                                                "
  echo "    Add the file meta_info.txt at ../meta_data/ folder. Example:      "
  echo "    meta_info.txt -> 'organ:reads_forward.fa.gz:reads_reverse.fa.gz'  "
  echo "    The reads must be in the ../input_data/ folder.                   "
  echo "                                                                      "
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
  echo "Running with $THREADS threads.";
  else
  echo "Running with $THREADS threads.";
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
  mapfile -t READS < ../meta_data/meta_info.txt
  #
  for read in "${READS[@]}" # 
    do
    #
    ORGAN_T=`echo $read | tr ':' '\t' | awk '{ print $1 }'`;
    SPL_Forward=`echo $read | tr ':' '\t' | awk '{ print $2 }'`;
    SPL_Reverse=`echo $read | tr ':' '\t' | awk '{ print $3 }'`;
    echo -e "\e[34m[TRACESPipe]\e[93m Running: Organ=$ORGAN_T Forward=$SPL_Forward Reverse=$SPL_Reverse\e[0m";
    #
    rm -f FW_READS.fq.gz RV_READS.fq.gz
    echo -e "\e[34m[TRACESPipe]\e[32m Copping an instance of the files ...\e[0m";
    cp ../input_data/$SPL_Forward FW_READS.fq.gz;
    cp ../input_data/$SPL_Reverse RV_READS.fq.gz;
    echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
    #
    # ==========================================================================
    # TRIM AND FILTER READS
    #
    CHECK_ADAPTERS;
    #
    echo -e "\e[34m[TRACESPipe]\e[32m Trimming and filtering with Trimmomatic ...\e[0m";
    ./TRACES_trim_filter_reads.sh $THREADS
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
      ./TRACES_remove_phix.sh $THREADS
      # IT IS USED ONLY FOR FALCON
      #
      # fastq_pair test_R1.fastq test_R2.fastq: [needs adaptation]
      # IF you want to remove Phix also before assembly
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Running viral metagenomic analysis with FALCON ...\e[0m";
      ./TRACES_metagenomics_viral.sh $ORGAN_T VDB.fa 10000 $THREADS
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Finding the best references ...\e[0m";
      ./TRACES_get_best_B19.sh $ORGAN_T > ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      #
      ./TRACES_get_best_HV1.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_HV2.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_HV3.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_HV4.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_HV5.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_HV6.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_HV6A.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_HV6B.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_HV7.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_HV8.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      #
      ./TRACES_get_best_POLY1.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_POLY2.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_POLY3.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_POLY4.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_POLY5.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_POLY6.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_POLY7.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_POLY8.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_POLY9.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_POLY10.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_POLY11.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_POLY12.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_POLY13.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_POLY14.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      # 
      ./TRACES_get_best_HBV.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_HPV.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_TTV.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_HBOV1.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_HBOVNOT1.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_VARV.sh $ORGAN_T >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
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
      echo -e "\e[34m[TRACESPipe]\e[32m Building complexity profiles with gto ...\e[0m";
      cat NP-o_fw_pr.fq NP-o_fw_unpr.fq NP-o_rv_pr.fq NP-o_rv_unpr.fq > P_TRACES_sample_reads.fq
      ./TRACES_profiles.sh GIS-$ORGAN_T VDB.fa P_TRACES_sample_reads.fq $ORGAN_T
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
      ./TRACES_metagenomics.sh $ORGAN_T DB.fa 12000 $THREADS 
      mkdir -p ../output_data/TRACES_results;
      #rm -f ../output_data/TRACES_results/*
      mv NV-$ORGAN_T.svg ../output_data/TRACES_results/
      mv NV-$ORGAN_T-HEAT.svg ../output_data/TRACES_results/
      mv REPORT_META_NON_VIRAL_$ORGAN_T.txt ../output_data/TRACES_results/
      cp top-non-viral-$ORGAN_T.csv  ../output_data/TRACES_results/
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      fi
    #
    # ==========================================================================
    # RUN SPECIFIC SPECIFIC ALIGN/CONSENSUS
    #
    if [[ "$RUN_SPECIFIC" -eq "1" ]];
      then
      echo -e "\e[34m[TRACESPipe]\e[32m Aliggning reads to specific viral ref(s) with pattern \"$SPECIFIC_ID\" using bowtie2 ...\e[0m";
      #
      echo "Extracting sequence with pattern \"$SPECIFIC_ID\" from VDB.fa ..."
      #
      CHECK_VDB;
      #
      gto_fasta_extract_read_by_pattern -p "$SPECIFIC_ID" < VDB.fa > SPECIFIC-$SPECIFIC_ID.fa
      echo "Aliggning ..."
      ./TRACES_viral_align_reads.sh SPECIFIC-$SPECIFIC_ID.fa $ORGAN_T $SPECIFIC_ID $THREADS
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
      ./TRACES_viral_consensus.sh SPECIFIC-$SPECIFIC_ID.fa viral_aligned_sorted-$ORGAN_T-$SPECIFIC_ID.bam $ORGAN_T $SPECIFIC_ID
      mkdir -p ../output_data/TRACES_specific_alignments;
      #rm -f ../output_data/TRACES_specific_alignments/*
      cp $ORGAN_T-$SPECIFIC_ID.fa ../output_data/TRACES_specific_alignments/
      cp $ORGAN_T-$SPECIFIC_ID.fa.fai ../output_data/TRACES_specific_alignments/
      mv viral_aligned_sorted-$ORGAN_T-$SPECIFIC_ID.bam ../output_data/TRACES_specific_alignments/
      mv viral_aligned_sorted-$ORGAN_T-$SPECIFIC_ID.bam.bai ../output_data/TRACES_specific_alignments/
      mkdir -p ../output_data/TRACES_specific_consensus;
      #rm -f ../output_data/TRACES_specific_consensus/*
      mv $SPECIFIC_ID-consensus-$ORGAN_T.fa ../output_data/TRACES_specific_consensus/
      mkdir -p ../output_data/TRACES_specific_bed;
      #rm -f ../output_data/TRACES_specific_bed/*
      mv $SPECIFIC_ID-calls-$ORGAN_T.bed ../output_data/TRACES_specific_bed/
      fi
    #
    # ========================================================================== 
    # DETAILED VIRAL ALIGN/CONSENSUS
    #
    if [[ "$RUN_B19_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "B19";
      fi
    #
    if [[ "$RUN_HV1_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HV1";
      fi
    #
    if [[ "$RUN_HV2_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HV2";
      fi
    #
    if [[ "$RUN_HV3_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HV3";
      fi
    #
    if [[ "$RUN_HV4_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HV4";
      fi
    #
    if [[ "$RUN_HV5_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HV5";
      fi
    #
    if [[ "$RUN_HV6_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HV6";
      fi
    #
    if [[ "$RUN_HV6A_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HV6A";
      fi
    #
    if [[ "$RUN_HV6B_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HV6B";
      fi
    #
    if [[ "$RUN_HV7_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HV7";
      fi
    #
    if [[ "$RUN_HV8_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HV8";
      fi
    #
    if [[ "$RUN_POLY1_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY1";
      fi
    #
    if [[ "$RUN_POLY2_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY2";
      fi
    #
    if [[ "$RUN_POLY2_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY2";
      fi
    #
    if [[ "$RUN_POLY3_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY3";
      fi
    #
    if [[ "$RUN_POLY4_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY4";
      fi
    #
    if [[ "$RUN_POLY5_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY5";
      fi
    #
    if [[ "$RUN_POLY6_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY6";
      fi
    #
    if [[ "$RUN_POLY7_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY7";
      fi
    #
    if [[ "$RUN_POLY8_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY8";
      fi
    #
    if [[ "$RUN_POLY9_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY9";
      fi
    #
    if [[ "$RUN_POLY10_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY10";
      fi
    #
    if [[ "$RUN_POLY11_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY11";
      fi
    #
    if [[ "$RUN_POLY12_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY12";
      fi
    #
    if [[ "$RUN_POLY13_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY13";
      fi
    #
    if [[ "$RUN_POLY14_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "POLY14";
      fi
    #
    if [[ "$RUN_TTV_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "TTV";
      fi
    #
    if [[ "$RUN_HBOV1_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HBOV1";
      fi
    #
    if [[ "$RUN_HBOVNOT1_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HBOVNOT1";
      fi
    #
    if [[ "$RUN_HBV_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HBV";
      fi
    #
    if [[ "$RUN_HPV_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "HPV";
      fi
    # 
    if [[ "$RUN_VARV_ON" -eq "1" ]];
      then
      ALIGN_AND_CONSENSUS "VARV";
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
      echo -e "\e[34m[TRACESPipe]\e[32m Aliggning reads to mitochondrial ref with bowtie2 ...\e[0m";
      ./TRACES_mt_align_reads.sh mtDNA.fa $ORGAN_T $THREADS
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
      ./TRACES_mt_consensus.sh mtDNA.fa mt_aligned_sorted-$ORGAN_T.bam $ORGAN_T
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
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m"
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
      echo -e "\e[34m[TRACESPipe]\e[32m Aliggning reads to Y-chromosome ref with bowtie2 ...\e[0m";
      ./TRACES_cy_align_reads.sh cy.fa $ORGAN_T $THREADS
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
      ./TRACES_cy_consensus.sh cy.fa cy_aligned_sorted-$ORGAN_T.bam $ORGAN_T
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
      ./TRACES_estimate_cy_quantity.sh $ORGAN_T $THREADS
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      fi
    #
    # ==========================================================================
    # DE-NOVO ASSEMBLY
    #
    if [[ "$RUN_DE_NOVO_ASSEMBLY" -eq "1" ]];
      then
      echo -e "\e[34m[TRACESPipe]\e[32m Running do-novo DNA assembly with SPAdes ...\e[0m";
      ./TRACES_assemble_all.sh $ORGAN_T $THREADS
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      fi
    #
    #
    # ==========================================================================
    done
    #
  # #
  #
  # ============================================================================ 
  #
  # BUILD COMPLETE VIRAL META TABLE FOR MULTIPLE ORGANS:
  #
  ./TRACES_get_report_meta.sh
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
