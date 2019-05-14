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
#
GEN_ADAPTERS=0;
GET_PHIX=0;
GET_MITO=0;
GET_B19=0;
GET_CY=0;
GET_HV3=0;
GET_HV4=0;
GET_HV7=0;
GET_AV=0;
GET_EXTRA=0;
#
RUN_ANALYSIS=0;
#
RUN_META_ON=0;
RUN_PROFILES_ON=0;
RUN_META_NON_VIRAL_ON=0;
RUN_MITO_ON=0;
RUN_B19_ON=0;
RUN_HV3_ON=0;
RUN_HV4_ON=0;
RUN_HV7_ON=0;
RUN_AV_ON=0;
RUN_CY_ON=0;
RUN_DE_NOVO_ASSEMBLY=0;
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
    -gb|--get-b19)
      GET_B19=1;
      SHOW_HELP=0;
      shift
    ;;
    -gh3|--get-hv3)
      GET_HV3=1;
      SHOW_HELP=0;
      shift
    ;;
    -gh4|--get-hv4)
      GET_HV4=1;
      SHOW_HELP=0;
      shift
    ;;
    -gh7|--get-hv7)
      GET_HV7=1;
      SHOW_HELP=0;
      shift
    ;;
    -gav|--get-av)
      GET_AV=1;
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
      RUN_MITO_ON=1;
      RUN_B19_ON=1;
      RUN_HV3_ON=1;
      RUN_HV4_ON=1;
      RUN_HV7_ON=1;
      RUN_AV_ON=1;
      RUN_CY_ON=1;
      #RUN_DE_NOVO_ASSEMBLY=1;
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
    -rb|--run-b19)
      RUN_ANALYSIS=1;
      RUN_B19_ON=1;
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
    -rh7|--run-hv7)
      RUN_ANALYSIS=1;
      RUN_HV7_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rav|--run-av)
      RUN_ANALYSIS=1;
      RUN_AV_ON=1;
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
      BUILD_VDB=1;
      BUILD_UDB=1;
      GEN_ADAPTERS=1;
      GET_PHIX=1;
      GET_MITO=1;
      GET_B19=1;
      GET_HV3=1;
      GET_HV4=1;
      GET_HV7=1;
      GET_AV=1;
      GET_CY=1;
      RUN_ANALYSIS=1;
      #
      RUN_META_ON=1;
      RUN_PROFILES_ON=1;
      RUN_META_NON_VIRAL_ON=1;
      RUN_MITO_ON=1;
      RUN_B19_ON=1;
      RUN_HV3_ON=1;
      RUN_HV4_ON=1;
      RUN_HV7_ON=1;
      RUN_AV_ON=1;
      RUN_CY_ON=1;
      RUN_DE_NOVO_ASSEMBLY=1;
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
    echo -e "    \e[32mAn automatic pipeline for viral genome identification\e[0m" 
    echo -e "    \e[32min the contexts of clinical virology and forensics\e[0m.         "
    echo "                                                                "
    echo -e "\e[93m    Usage: ./ki.sh [options]                                    \e[0m"
    echo "                                                                "
    echo "    -h,    --help           Show this help message and exit,     "
    echo "    -i,    --install        Installation of all the tools,       "
    echo "    -vdb,  --build-viral    Build viral database,                "
    echo "    -udb,  --build-unviral  Build non viral database (control),  "
    echo "                                                                 "
    echo "    -gad,  --gen-adapters   Generate FASTA file with adapters,   "
    echo "    -gp,   --get-phix       Extracts PhiX genomes (Needs viral DB),  "
    echo "    -gm,   --get-mito       Downloads human Mitochondrial genome,"
    echo "    -gb,   --get-b19        Extracts B19 genome (Needs viral DB),"
    echo "    -gh3,  --get-hv3        Extracts HV3 genome (Needs viral DB),"
    echo "    -gh4,  --get-hv4        Extracts HV4 genome (Needs viral DB),"
    echo "    -gh7,  --get-hv7        Extracts HV7 genome (Needs viral DB),"
    echo "    -gav,  --get-av         Extracts AV genome (Needs viral DB),"
    echo "    -gy,   --get-y-chromo   Downloads human Y-chromosome,        "
    echo "    -gx,   --get-extra-vir  Downloads/appends (VDB) extra viral seq, "
    echo "                                                                 "
    echo "    -ra,   --run-analysis   Run data analysis,                   "
    echo "                                                                 "
    echo "    -rm,   --run-meta       Run viral metagenomic identification,    "
    echo "    -ro,   --run-meta-nv    Run NON-viral metagenomic identification,   "
    echo "                                                                 "
    echo "    -rmt,  --run-mito       Run Mito align, sort and consensus seq,   "
    echo "    -rb,   --run-b19        Run B19 align, sort and consensus seq,    "
    echo "    -rh3,  --run-hv3        Run HV3 align, sort and consensus seq,    "
    echo "    -rh4,  --run-hv4        Run HV4 align, sort and consensus seq,    "
    echo "    -rh7,  --run-hv7        Run HV7 align, sort and consensus seq,    "
    echo "    -rav,  --run-av         Run AV align, sort and consensus seq,    "
    echo "                                                                "
    echo "    -rcy,  --run-y-chromo   Run Y align, sort and consensus seq,    "
    echo "    -rda,  --run-de-novo    Run de-novo assembly,               "
    echo "                                                                "
    echo "    -all,  --run-all        Run all the options.                 "
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
if [[ "$GET_B19" -eq "1" ]];
  then
  ./ki_get_b19.sh
  fi
# ==============================================================================
#
if [[ "$GET_HV3" -eq "1" ]];
  then
  ./ki_get_hv3.sh
  fi
# ==============================================================================
#
if [[ "$GET_HV4" -eq "1" ]];
  then
  ./ki_get_hv4.sh
  fi
# ==============================================================================
#
if [[ "$GET_HV7" -eq "1" ]];
  then
  ./ki_get_hv7.sh
  fi
# ==============================================================================
#
if [[ "$GET_AV" -eq "1" ]];
  then
  ./ki_get_av.sh
  fi
#
# ==============================================================================
#
if [[ "$GET_CY" -eq "1" ]];
  then
  ./ki_get_cy.sh
  fi
#
# ==============================================================================
#
if [[ "$GET_EXTRA" -eq "1" ]];
  then
  ./ki_get_enriched_sequences.sh VDB.fa
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
    # ==========================================================================
    # TRIM AND FILTER READS
    #
    echo -e "\e[34m[ki]\e[32m Trimming and filtering with Trimmomatic ...\e[0m";
    ./ki_trim_filter_reads.sh
    echo -e "\e[34m[ki]\e[32m Done!\e[0m";
    #
    # THE OUTPUT OF TRIMMING IS:
    # o_fw_pr.fq  o_fw_unpr.fq  o_rv_pr.fq  o_rv_unpr.fq
    #
    # ========================================================================== 
    # METAGENOMICS USING ONLY A VIRAL DATABASE
    #
    if [[ "$RUN_META_ON" -eq "1" ]];
      then
      echo -e "\e[34m[ki]\e[32m Removing PhiX from the samples with MAGNET ...\e[0m";
      ./ki_remove_phix.sh # IT IS USED ONLY FOR FALCON
      echo -e "\e[34m[ki]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[ki]\e[32m Running viral metagenomic analysis with FALCON ...\e[0m";
      ./ki_metagenomics.sh $ORGAN_T VDB.fa
      echo -e "\e[34m[ki]\e[32m Done!\e[0m";
      fi
    #
    # ==========================================================================
    # COMPLEXITY PROFILES FOR TOP IDENTIFIED VIRAL SEQUENCES
    #
    if [[ "$RUN_PROFILES_ON" -eq "1" ]];
      then
      echo -e "\e[34m[ki]\e[32m Building complexity profiles with gto ...\e[0m";
      cat NP-o_fw_pr.fq NP-o_fw_unpr.fq NP-o_rv_pr.fq NP-o_rv_unpr.fq > ki_sample_reads.fq
      ./ki_profiles.sh GIS-$ORGAN_T VDB.fa ki_sample_reads.fq $ORGAN_T
      echo -e "\e[34m[ki]\e[32m Done!\e[0m";
      fi
    #
    # ========================================================================== 
    # METAGENOMICS WITH ALL NON VIRAL
    #
    if [[ "$RUN_META_NON_VIRAL_ON" -eq "1" ]];
      then
      echo -e "\e[34m[ki]\e[32m Running NON viral metagen. analysis with FALCON ...\e[0m";
      ./ki_metagenomics.sh $ORGAN_T-NON_VIRAL DB.fa 
      echo -e "\e[34m[ki]\e[32m Done!\e[0m";
      fi
    #
    # ==========================================================================
    # MITOCHONDRIAL GENOME ALIGN AND CONSENSUS
    #
    if [[ "$RUN_MITO_ON" -eq "1" ]];
      then
      echo -e "\e[34m[ki]\e[32m Aliggning reads to mitochondrial ref with bowtie2 ...\e[0m";
      ./ki_mt_align_reads.sh mtDNA.fa $ORGAN_T
      echo -e "\e[34m[ki]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[ki]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
      ./ki_mt_consensus.sh mtDNA.fa mt_aligned_sorted-$ORGAN_T.bam $ORGAN_T
      echo -e "\e[34m[ki]\e[32m Done!\e[0m"
      fi
    #
    # ========================================================================== 
    # SPECIFIC VIRAL ALIGN AND CONSENSUS: B19, HV3, HV4, HV7, AV
    #
    if [[ "$RUN_B19_ON" -eq "1" ]];
      then
      echo -e "\e[34m[ki]\e[32m Aliggning reads to B19 ref with bowtie2 ...\e[0m";
      ./ki_b19_align_reads.sh B19.fa $ORGAN_T
      echo -e "\e[34m[ki]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[ki]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
      ./ki_b19_consensus.sh B19.fa b19_aligned_sorted-$ORGAN_T.bam $ORGAN_T
      echo -e "\e[34m[ki]\e[32m Done!\e[0m"
      fi
    #
    if [[ "$RUN_HV3_ON" -eq "1" ]];
      then
      echo -e "\e[34m[ki]\e[32m Aliggning reads to HV3 ref with bowtie2 ...\e[0m";
      ./ki_hv3_align_reads.sh HV3.fa $ORGAN_T
      echo -e "\e[34m[ki]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[ki]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
      ./ki_hv3_consensus.sh HV3.fa hv3_aligned_sorted-$ORGAN_T.bam $ORGAN_T
      echo -e "\e[34m[ki]\e[32m Done!\e[0m"
      fi
    #
    if [[ "$RUN_HV4_ON" -eq "1" ]];
      then
      echo -e "\e[34m[ki]\e[32m Aliggning reads to HV4 ref with bowtie2 ...\e[0m";
      ./ki_hv4_align_reads.sh HV4.fa $ORGAN_T
      echo -e "\e[34m[ki]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[ki]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
      ./ki_hv4_consensus.sh HV4.fa hv4_aligned_sorted-$ORGAN_T.bam $ORGAN_T
      echo -e "\e[34m[ki]\e[32m Done!\e[0m"
      fi
    #
    if [[ "$RUN_HV7_ON" -eq "1" ]];
      then
      echo -e "\e[34m[ki]\e[32m Aliggning reads to HV7 ref with bowtie2 ...\e[0m";
      ./ki_hv7_align_reads.sh HV7.fa $ORGAN_T
      echo -e "\e[34m[ki]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[ki]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
      ./ki_hv7_consensus.sh HV7.fa hv7_aligned_sorted-$ORGAN_T.bam $ORGAN_T
      echo -e "\e[34m[ki]\e[32m Done!\e[0m"
      fi	      
    #
    if [[ "$RUN_AV_ON" -eq "1" ]];
      then
      echo -e "\e[34m[ki]\e[32m Aliggning reads to AV ref with bowtie2 ...\e[0m";
      ./ki_av_align_reads.sh AV.fa $ORGAN_T
      echo -e "\e[34m[ki]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[ki]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
      ./ki_av_consensus.sh AV.fa av_aligned_sorted-$ORGAN_T.bam $ORGAN_T
      echo -e "\e[34m[ki]\e[32m Done!\e[0m"
      fi   
    #
    # ========================================================================== 
    # CY VERYFICATION
    #
    if [[ "$RUN_CY_ON" -eq "1" ]];
      then
      echo -e "\e[34m[ki]\e[32m Extracting Y chromosome reads with MAGNET ...\e[0m";
      ./ki_extract_classify_cy.sh $ORGAN_T
      echo -e "\e[34m[ki]\e[32m Done!\e[0m";
      fi
    #
    # ==========================================================================
    # DE-NOVO ASSEMBLY
    #
    if [[ "$RUN_DE_NOVO_ASSEMBLY" -eq "1" ]];
      then
      echo -e "\e[34m[ki]\e[32m Running do-novo DNA assembly with SPAdes ...\e[0m";
      ./ki_assemble_all.sh $ORGAN_T
      echo -e "\e[34m[ki]\e[32m Done!\e[0m";
      fi
    #
    # ==========================================================================
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
