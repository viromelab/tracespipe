#!/bin/bash
#
##################################################################################
# ============================================================================== #
# =                                                                            = #
# =                                    KI                                      = #
# =                                                                            = #
# =          An automatic pipeline for viral genome identification             = #
# =           in the contexts of clinical virology and forensics.              = #
# =                                                                            = #
# ============================================================================== #
##################################################################################
#
SHOW_HELP=0;
INSTALL=0;
BUILD_VDB=0;
BUILD_UDB=0;
GEN_ADAPTERS=0;
GET_PHIX=0;
GET_MITO=0;
GET_CY=0;
RUN_ANALYSIS=0;
#
RUN_META_ON=0;
RUN_PROFILES_ON=0;
RUN_META_NON_VIRAL_ON=0;
RUN_MITO_ON=0;
RUN_MITO_CONSENSUS=0;
RUN_CY_ON=0;
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
    -i|--install)
      INSTALL=1;
      SHOW_HELP=0;
      shift
    ;;
    -vdb|--build-viral)
      BUILD_VDB=1;
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
    -ra|--run-analysis)
      RUN_ANALYSIS=1;
      RUN_META_ON=1;
      RUN_PROFILES_ON=1;
      RUN_META_NON_VIRAL_ON=1;
      RUN_MITO_ON=1;
      RUN_MITO_CONSENSUS=1;
      RUN_CY_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rm|--run-mito)
      RUN_ANALYSIS=1;
      RUN_MITO_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rmc|--run-mito-cons)
      RUN_ANALYSIS=1;
      RUN_MITO_ON=1;
      RUN_MITO_CONSENSUS=1;
      SHOW_HELP=0;
      shift
    ;;
    -all|--run-all)
      INSTALL=1;
      BUILD_VDB=1;
      BUILD_UDB=1;
      GEN_ADAPTERS=1;
      GET_PHIX=1;
      GET_MITO=1;
      GET_CY=1;
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
    echo "    -udb, --build-unviral  Build non viral database (control),  "
    echo "    -gad, --gen-adapters   Generate FASTA file with adapters,   "
    echo "    -gp,  --get-phix       Downloads PhiX genomes,              "
    echo "    -gm,  --get-mito       Downloads human Mitochondrial genome,"
    echo "    -gy,  --get-y-chromo   Downloads human Y-chromosome,        "
    echo "    -ra,  --run-analysis   Run data analysis,                   "
    echo "    -rm,  --run-mito       Run Mito align and sort (BAM),       "
    echo "    -rmc, --run-mito-cons  Run Mito align, sort and consensus seq,   "
    echo "                                                                "
    echo "    -all, --run-all        Run all the options.                 "
    echo "                                                                "
    echo -e "\e[93m    Example: ./ki.sh -all                                         \e[0m"
    echo "                                                                "
    echo "    meta_info.txt -> 'name:reads_forward.fa.gz:reads_reverse.fa.gz'  "
    echo "    The reads and meta_info.txt must be in the src/ folder.     "
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
  gto_build_dbs.sh -vi
  gunzip VDB.fa.gz
  fi
#
# ==============================================================================
#
if [[ "$BUILD_UDB" -eq "1" ]];
  then
  gto_build_dbs.sh -ba -ar -pr -fu -pl -in -mi -ps 
  # -vm -vo
  zcat BDB.fa.gz ADB.fa.gz PDB.fa.gz FDB.fa.gz TDB.fa.gz TDB.fa.gz IDB.fa.gz MTDB.fa.gz PLDB.fa.gz > DB.fa
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
# ==============================================================================
#
if [[ "$GET_CY" -eq "1" ]];
  then
  ./ki_get_cy.sh
  fi
#
# ==============================================================================
#
if [[ "$RUN_ANALYSIS" -eq "1" ]];
  then
  #
  mapfile -t READS < meta_info.txt
  #
  for read in "${READS[@]}" # 
    do
    #
    ORGAN_T=`echo $read | tr ':' '\t' | awk '{ print $1 }'`;
    SPL_Forward=`echo $read | tr ':' '\t' | awk '{ print $2 }'`;
    SPL_Reverse=`echo $read | tr ':' '\t' | awk '{ print $3 }'`;
    echo -e "\e[34m[ki]\e[93m Running: Organ=$ORGAN_T Forward=$SPL_Forward Reverse=$SPL_Reverse\e[0m";
    #
    rm -f FW_READS.fq.gz RV_READS.fq.gz
    echo -e "\e[34m[ki]\e[32m Copping an instance of the files ...\e[0m";
    cp $SPL_Forward FW_READS.fq.gz;
    cp $SPL_Reverse RV_READS.fq.gz;
    echo -e "\e[34m[ki]\e[32m Done!\e[0m";
    #
    echo -e "\e[34m[ki]\e[32m Trimming and filtering with Trimmomatic ...\e[0m";
    ./ki_trim_filter_reads.sh
    echo -e "\e[34m[ki]\e[32m Done!\e[0m";
    #
    # o_fw_pr.fq  o_fw_unpr.fq  o_rv_pr.fq  o_rv_unpr.fq
    #
    if [[ "$RUN_META_ON" -eq "1" ]];
      then
      echo -e "\e[34m[ki]\e[32m Removing PhiX from the samples with MAGNET ...\e[0m";
      ./ki_remove_phix.sh
      # CHANGE MAGNET: COMPATIBILITY PROBLEMS WITH SPADES AND PROBABLY BOWTIE2
      echo -e "\e[34m[ki]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[ki]\e[32m Running viral metagenomic analysis with FALCON ...\e[0m";
      ./ki_metagenomics.sh $ORGAN_T VDB.fa
      echo -e "\e[34m[ki]\e[32m Done!\e[0m";
      fi
    #
    if [[ "$RUN_PROFILES_ON" -eq "1" ]];
      then
      echo -e "\e[34m[ki]\e[32m Building complexity profiles with gto ...\e[0m";
      cat NP-o_fw_pr.fq NP-o_fw_unpr.fq NP-o_rv_pr.fq NP-o_rv_unpr.fq > ki_sample_reads.fq
      ./ki_profiles.sh GIS-$ORGAN_T VDB.fa ki_sample_reads.fq $ORGAN_T
      echo -e "\e[34m[ki]\e[32m Done!\e[0m";
      fi
    #
    if [[ "$RUN_META_NON_VIRAL_ON" -eq "1" ]];
      then
      echo -e "\e[34m[ki]\e[32m Running NON viral metagen. analysis with FALCON ...\e[0m";
      ./ki_metagenomics.sh $ORGAN_T-NON_VIRAL DB.fa 
      echo -e "\e[34m[ki]\e[32m Done!\e[0m";
      fi
    #
    if [[ "$RUN_MITO_ON" -eq "1" ]];
      then
      echo -e "\e[34m[ki]\e[32m Aliggning reads to mitochondrial ref with bowtie2 ...\e[0m";
      ./ki_align_reads.sh mtDNA.fa $ORGAN_T
      echo -e "\e[34m[ki]\e[32m Done!\e[0m";
      #
      if [[ "$RUN_MITO_CONSENSUS" -eq "1" ]];
        then
        echo -e "\e[34m[ki]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./ki_consensus.sh mtDNA.fa aligned_sorted-$ORGAN_T.bam $ORGAN_T
        echo -e "\e[34m[ki]\e[32m Done!\e[0m"
	fi
      #
   #   #
   #   echo -e "\e[34m[ki]\e[32m Extracting mitochondrial reads with MAGNET ...\e[0m";
   #   ./ki_extract_mito.sh
   #   echo -e "\e[34m[ki]\e[32m Done!\e[0m";
   #   #
   #   echo -e "\e[34m[ki]\e[32m Running mitochondrial DNA assembly with SPAdes ...\e[0m";
   #   ./ki_assemble_mito.sh $ORGAN_T
   #   echo -e "\e[34m[ki]\e[32m Done!\e[0m";
      fi
    #
    if [[ "$RUN_CY_ON" -eq "1" ]];
      then
      echo -e "\e[34m[ki]\e[32m Extracting Y chromosome reads with MAGNET ...\e[0m";
      ./ki_extract_classify_cy.sh $ORGAN_T
      echo -e "\e[34m[ki]\e[32m Done!\e[0m";
      fi
    #
    done
  #
  mkdir -p results;
  rm -f results/*
  mv *.pdf results/
  mv *.svg results/
  #
  fi
#
# ==============================================================================
################################################################################
