#!/bin/bash
#
##################################################################################
# ============================================================================== #
# =                                                                            = #
# =                                 T R A C E S                                = #
# =                                                                            = #
# =       An automatic pipeline for viral and human genome identification      = #
# =           in the contexts of clinical virology and forensics.              = #
# =                                                                            = #
# ============================================================================== #
##################################################################################
#
SHOW_HELP=0;
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
RUN_TTV_ON=0;
RUN_JCV_ON=0;
RUN_MCV_ON=0;
RUN_BK_ON=0;
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
      RUN_MITO_ON=1;
      RUN_B19_ON=1;
      RUN_HV3_ON=1;
      RUN_HV4_ON=1;
      RUN_HV5_ON=1;
      RUN_HV6_ON=1;
      RUN_HV7_ON=1;
      RUN_TTV_ON=1;
      RUN_JCV_ON=1;
      RUN_MCV_ON=1;
      RUN_CY_ON=1;
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
      RUN_TTV_ON=1;
      RUN_JCV_ON=1;
      RUN_MCV_ON=1;
      RUN_BK_ON=1;
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
    -rtt|--run-ttv)
      RUN_ANALYSIS=1;
      RUN_TTV_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rjc|--run-jcv)
      RUN_ANALYSIS=1;
      RUN_JCV_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rmc|--run-mcv)
      RUN_ANALYSIS=1;
      RUN_MCV_ON=1;
      SHOW_HELP=0;
      shift
    ;;
    -rbk|--run-bk)
      RUN_ANALYSIS=1;
      RUN_BK_ON=1;
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
      RUN_TTV_ON=1;
      RUN_JCV_ON=1;
      RUN_MCV_ON=1;
      RUN_BK_ON=1;
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
    echo "For help, try: ./TRACES.sh -h"
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
    echo -e "\e[34m                                                              "
    echo "     ████████╗ ██████╗   █████╗   ██████╗ ███████╗ ███████╗         "
    echo "     ╚══██╔══╝ ██╔══██╗ ██╔══██╗ ██╔════╝ ██╔════╝ ██╔════╝         "
    echo "        ██║    ██████╔╝ ███████║ ██║      █████╗   ███████╗         "
    echo "        ██║    ██╔══██╗ ██╔══██║ ██║      ██╔══╝   ╚════██║         "
    echo "        ██║    ██║  ██║ ██║  ██║ ╚██████╗ ███████╗ ███████║         "
    echo "        ╚═╝    ╚═╝  ╚═╝ ╚═╝  ╚═╝  ╚═════╝ ╚══════╝ ╚══════╝         "
    echo "                                                                "
    echo -e "    \e[32mAn automatic pipeline for viral and human genome analysis\e[0m" 
    echo -e "    \e[32min the contexts of clinical virology and forensics\e[0m.         "
    echo "                                                                "
    echo -e "\e[93m    Usage: ./TRACES.sh [options]                                \e[0m"
    echo "                                                                "
    echo "    -h,    --help            Show this help message and exit,     "
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
    echo "    -ra,   --run-analysis    Run data analysis,                   "
    echo "                                                                  "
    echo "    -rm,   --run-meta        Run viral metagenomic identification,    "
    echo "    -ro,   --run-meta-nv     Run NON-viral metagenomic identification,   "
    echo "                                                                  "
    echo "    -rava, --run-all-v-alig  Run all viral align/sort/consensus seqs,    "
    echo "                                                                 "
    echo "    -rb,   --run-b19         Run B19 align, sort and consensus seq,    "
    echo "    -rh1,  --run-hv1         Run HV1 align, sort and consensus seq,    "
    echo "    -rh2,  --run-hv2         Run HV2 align, sort and consensus seq,    "
    echo "    -rh3,  --run-hv3         Run HV3 align, sort and consensus seq,    "
    echo "    -rh4,  --run-hv4         Run HV4 align, sort and consensus seq,    "
    echo "    -rh5,  --run-hv5         Run HV5 align, sort and consensus seq,    "
    echo "    -rh6,  --run-hv6         Run HV6 align, sort and consensus seq,    "
    echo "    -rh6a, --run-hv6a        Run HV6A align, sort and consensus seq,    "
    echo "    -rh6b, --run-hv6b        Run HV6B align, sort and consensus seq,    "
    echo "    -rh7,  --run-hv7         Run HV7 align, sort and consensus seq,    "
    echo "    -rh8,  --run-hv8         Run HV8 align, sort and consensus seq,    "
    echo "    -rtt,  --run-ttv         Run TTV align, sort and consensus seq,    "
    echo "    -rjc,  --run-jcv         Run JCV align, sort and consensus seq,    "
    echo "    -rmc,  --run-mcv         Run MCV align, sort and consensus seq,    "
    echo "    -rbk,  --run-bk          Run BK align, sort and consensus seq,    "
    echo "    -rbv1, --run-hbov1       Run HBoV1 align, sort and consensus seq,    "
    echo "    -rbv0, --run-hbovnot1    Run HBoV 2,3,... align, sort and consensus seq,    "
    echo "    -rhbv, --run-hbv         Run HBV align, sort and consensus seq,    "
    echo "    -rhpv, --run-hpv         Run HPV align, sort and consensus seq,    "
    echo "    -rvar, --run-varv        Run VARV align, sort and consensus seq,    "
    echo "                                                                 "
    echo "    -rsr,  --run-specific    Run specific REF align/consensus seq, "
    echo "                                                                 "
    echo "    -rmt,  --run-mito        Run Mito align, sort and consensus seq,   "
    echo "    -rcy,  --run-y-chromo    Run CY align, sort and consensus seq,    "
    echo "                                                                  "
    echo "    -rda,  --run-de-novo     Run de-novo assembly,               "
    echo "                                                                 "
    echo "    -all,  --run-all         Run all the options.                 "
    echo "                                                                "
    echo -e "\e[93m    Example: ./TRACES.sh --run-meta --run-mito           \e[0m"
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
      echo -e "\e[34m[TRACES]\e[32m Removing PhiX from the samples with MAGNET ...\e[0m";
      ./TRACES_remove_phix.sh # IT IS USED ONLY FOR FALCON
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[TRACES]\e[32m Running viral metagenomic analysis with FALCON ...\e[0m";
      ./TRACES_metagenomics_viral.sh $ORGAN_T VDB.fa 10000
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[TRACES]\e[32m Finding the best references ...\e[0m";
      ./TRACES_get_best_B19.sh $ORGAN_T > REPORT_META_VIRAL_$ORGAN_T.txt
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
      ./TRACES_get_best_HBV.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_MCV.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_JCV.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_HBoV1.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_BK.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_HPV.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
      ./TRACES_get_best_TTV.sh $ORGAN_T >> REPORT_META_VIRAL_$ORGAN_T.txt
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
      echo -e "\e[34m[TRACES]\e[32m Running NON viral metagenomic analysis with FALCON ...\e[0m";
      if [ ! -f DB.fa ]; 
        then
	echo -e "\e[31mERROR: Non-viral database FASTA file (DB.fa) not found!\e[0m"
	echo "TIP: first run ./TRACES --build-unviral"
	exit 1;
        fi
      ./TRACES_metagenomics.sh $ORGAN_T DB.fa 10000 
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
      fi
    #
    # ==========================================================================
    # MITOCHONDRIAL GENOME ALIGN AND CONSENSUS
    #
    if [[ "$RUN_MITO_ON" -eq "1" ]];
      then
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
    # RUN SPECIFIC SPECIFIC ALIGN/CONSENSUS
    #
    if [[ "$RUN_SPECIFIC" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to B19 ref with bowtie2 ...\e[0m";
      echo "Extracting sequence from VDB.fa ..."
      ###gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > SPECIFIC-$IDN.fa
      echo "Aliggning ..."
      ###./TRACES_b19_align_reads.sh SPECIFIC-$IDN.fa $ORGAN_T
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
      ###./TRACES_b19_consensus.sh SPECIFIC-$IDN.fa b19_aligned_sorted-$ORGAN_T.bam $ORGAN_T
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
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > B19.fa
        echo "Aliggning ..."
        ./TRACES_b19_align_reads.sh B19.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_b19_consensus.sh B19.fa b19_aligned_sorted-$ORGAN_T.bam $ORGAN_T
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
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > HV1.fa
        echo "Aliggning ..."
        ./TRACES_hv1_align_reads.sh HV1.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_hv1_consensus.sh HV1.fa hv1_aligned_sorted-$ORGAN_T.bam $ORGAN_T
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
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > HV2.fa
        echo "Aliggning ..."
        ./TRACES_hv2_align_reads.sh HV2.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_hv2_consensus.sh HV2.fa hv1_aligned_sorted-$ORGAN_T.bam $ORGAN_T
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
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > HV3.fa
        echo "Aliggning ..."
        ./TRACES_hv3_align_reads.sh HV3.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_hv3_consensus.sh HV3.fa hv3_aligned_sorted-$ORGAN_T.bam $ORGAN_T
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
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > HV4.fa
        echo "Aliggning ..."
        ./TRACES_hv4_align_reads.sh HV4.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_hv4_consensus.sh HV4.fa hv4_aligned_sorted-$ORGAN_T.bam $ORGAN_T
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
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > HV5.fa
        echo "Aliggning ..."
        ./TRACES_hv5_align_reads.sh HV5.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_hv4_consensus.sh HV5.fa hv5_aligned_sorted-$ORGAN_T.bam $ORGAN_T
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
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > HV6.fa
        echo "Aliggning ..."
        ./TRACES_hv6_align_reads.sh HV6.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_hv6_consensus.sh HV6.fa hv6_aligned_sorted-$ORGAN_T.bam $ORGAN_T
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
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > HV6A.fa
        echo "Aliggning ..."
        ./TRACES_hv6A_align_reads.sh HV6A.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_hv6a_consensus.sh HV6A.fa hv6a_aligned_sorted-$ORGAN_T.bam $ORGAN_T
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
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > HV6B.fa
        echo "Aliggning ..."
        ./TRACES_hv6B_align_reads.sh HV6B.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_hv6b_consensus.sh HV6B.fa hv6b_aligned_sorted-$ORGAN_T.bam $ORGAN_T
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
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > HV7.fa
        echo "Aliggning ..."
        ./TRACES_hv7_align_reads.sh HV7.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_hv7_consensus.sh HV7.fa hv7_aligned_sorted-$ORGAN_T.bam $ORGAN_T
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
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > HV8.fa
        echo "Aliggning ..."
        ./TRACES_hv8_align_reads.sh HV8.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_hv8_consensus.sh HV8.fa hv8_aligned_sorted-$ORGAN_T.bam $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
        fi
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
      fi
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
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > TTV.fa
        echo "Aliggning ..."
        ./TRACES_ttv_align_reads.sh TTV.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_ttv_consensus.sh TTV.fa ttv_aligned_sorted-$ORGAN_T.bam $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
        fi
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
      fi
    #
    if [[ "$RUN_JCV_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to JCV ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_JCV.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]];
        then
        echo "Extracting sequence from JCV.fa ..."
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > JCV.fa
        echo "Aliggning ..."
        ./TRACES_jcv_align_reads.sh JCV.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_jcv_consensus.sh JCV.fa jcv_aligned_sorted-$ORGAN_T.bam $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
        fi
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
      fi
    #
    if [[ "$RUN_MCV_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to MCV ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_MCV.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]];
        then
        echo "Extracting sequence from MCV.fa ..."
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > MCV.fa
        echo "Aliggning ..."
        ./TRACES_mcv_align_reads.sh MCV.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_mcv_consensus.sh MCV.fa mcv_aligned_sorted-$ORGAN_T.bam $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
        fi
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
      fi
    #
    if [[ "$RUN_BK_ON" -eq "1" ]];
      then
      echo -e "\e[34m[TRACES]\e[32m Aliggning reads to MK ref with bowtie2 ...\e[0m";
      V_INFO=`./TRACES_get_best_MK.sh $ORGAN_T`;
      echo "Best match: $V_INFO";
      V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
      if [[ "$V_GID" != "-" ]];
        then
        echo "Extracting sequence from MK.fa ..."
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > MK.fa
        echo "Aliggning ..."
        ./TRACES_bk_align_reads.sh MK.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_bk_consensus.sh BK.fa bk_aligned_sorted-$ORGAN_T.bam $ORGAN_T
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
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > HBoV1.fa
        echo "Aliggning ..."
        ./TRACES_hbov1_align_reads.sh HBoV1.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_hbov1_consensus.sh HBoV1.fa hbov1_aligned_sorted-$ORGAN_T.bam $ORGAN_T
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
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > HBoVnot1.fa
        echo "Aliggning ..."
        ./TRACES_hbovnot1_align_reads.sh HBoVnot1.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_hbovnot1_consensus.sh HBoVnot1.fa hbovnot1_aligned_sorted-$ORGAN_T.bam $ORGAN_T
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
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > HBV.fa
        echo "Aliggning ..."
        ./TRACES_hbv_align_reads.sh HBV.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_hbv_consensus.sh HBV.fa hbv_aligned_sorted-$ORGAN_T.bam $ORGAN_T
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
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > HPV.fa
        echo "Aliggning ..."
        ./TRACES_hpv_align_reads.sh HPV.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_hpv_consensus.sh HPV.fa hpv_aligned_sorted-$ORGAN_T.bam $ORGAN_T
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
        gto_fasta_extract_read_by_pattern -p "$V_GID" < VDB.fa > VARV.fa
        echo "Aliggning ..."
        ./TRACES_varv_align_reads.sh VARV.fa $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACES]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
        ./TRACES_varv_consensus.sh VARV.fa varv_aligned_sorted-$ORGAN_T.bam $ORGAN_T
        echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
        fi
      echo -e "\e[34m[TRACES]\e[32m Done!\e[0m"
      fi
    #
    # ========================================================================== 
    # CY VERYFICATION
    #
    if [[ "$RUN_CY_ON" -eq "1" ]];
      then
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
  # RESULTS WITH REPORTS AND IMAGES
  mkdir -p TRACES_results;
  rm -f TRACES_results/*
  mv *.pdf TRACES_results/
  mv *.svg TRACES_results/
  mv REPORT_META_VIRAL_*.txt TRACES_results/
  #
  # CONSENSUS FILES
  mkdir -p TRACES_consensus;
  rm -f TRACES_consensus/*
  mv *-consensus-*.fa TRACES_consensus/
  #
  # ALIGNMENT FILES
  mkdir -p TRACES_alignments;
  rm -f TRACES_alignments/*
  #mv *.fa TRACES_alignments/ #XXX: THIS WILL MOVE THE DATABASE!
  #mv *.fa.fai TRACES_alignments/
  mv *_aligned_sorted-*.bam TRACES_alignments/
  mv *_aligned_sorted-*.bam.bai TRACES_alignments/
  #
  # BED FILES
  mkdir -p TRACES_bed;
  rm -f TRACES_bed/*
  mv *-calls-*.bed TRACES_bed/
  #
  # BUILD COMPLETE VIRAL META TABLE FOR MULTIPLE ORGANS:
  ./TRACES_get_report_meta.sh
  #
  fi
#
# ==============================================================================
################################################################################
