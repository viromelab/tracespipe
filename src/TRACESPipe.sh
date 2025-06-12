#!/bin/bash
#
##################################################################################
# ============================================================================== #
# =                                                                            = #
# =                             T R A C E S P i p e                            = #
# =                                                                            = #
# =          A Next-Generation Sequencing pipeline for reconstruction          = #
# =        and analysis of viral and host genomes at multi-organ level.        = #
# =                                                                            = #
# ============================================================================== #
##################################################################################
#
SOURCE_DIR="$(dirname "$(readlink -f "$0")")"
#
SHOW_HELP=0;
SHOW_VERSION=0;
FORCE=0;
FLUSH_LOGS=0;
GET_THREADS=0;
THREADS=0;
#
INSTALL=0;
SHOW_PROG_VER=0;
UPDATE=0;
#
SAMPLE_TEST=0;
#
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
DOWNLOAD_MITO_SPECIES=0;
DOWNLOAD_MITO_POPULATION=0;
#
RUN_MITO_SPECIES=0;
RUN_MITO_POPULATION=0;
#
CREATE_BLAST_DB=0;
UPDATE_BLAST_DB=0;
SEARCH_BLAST_DB=0;
SEARCH_BLAST_REMOTE_DB=0;
BLAST_QUERY="";
#
FLUSH_OUTPUT=0;
RUN_ANALYSIS=0;
#
RUN_META_ON=0;
RUN_BEST_OF_BESTS=0;
RUN_PROFILES_ON=0;
RUN_META_NON_VIRAL_ON=0;
RUN_GID_COMPLEXITY_PROFILE=0;
GID_COMPLEXITY_PROFILE="";
COMPLEXITY_PROFILE_WINDOW="5";
COMPLEXITY_PROFILE_LEVEL="1";
#
RUN_MITO_ON=0;
RUN_MITO_DAMAGE_ON=0;
#
RUN_GID_DAMAGE_ANALYSIS=0;
GID_DAMAGE="";
#
VIEW_TOP=0;
#
REMOVE_DUPLICATIONS=0;
HIGH_SENSITIVITY=0;
#
RUN_ALL_VIRAL=0;
#
RUN_VISUAL_ALIGN=0;
#
RUN_COVERAGE_TABLE=0;
RUN_COVERAGE_TABLE_CSV=0;
RUN_COVERAGE_PROFILE=0;
COVERAGE_NAME="";
MAX_COVERAGE_PROFILE=0;
COVERAGE_MIN_X="0";
COVERAGE_LOG_SCALE="";
COVERAGE_WINDOW_SIZE="5";
COVERAGE_DROP="1";
#
RUN_CHANGE_MT=0;
#
RUN_DECRYPT=0;
RUN_ENCRYPT=0;
#
RUN_SPECIFIC=0;
RUN_DENOVO_SPECIFIC=0;
RUN_SPECIFIC_SENSITIVE=0;
#
RUN_MULTIORGAN_CONSENSUS=0;
#
RUN_CY_ON=0;
RUN_CY_QUANT_ON=0;
#
RUN_DE_NOVO_ASSEMBLY=0;
#
RUN_HYBRID=0;
#
RUN_DIFF=0;
RUN_SPECIFIC_DIFF=0;
RUN_DIFF_VIRUS="";
RUN_DIFF_ID="";
#
RUN_BLAST_RECONSTRUCTED=0;
#
B_ALT_VIRAL_DB=0;
VIRAL_DATABASE_METADATA="";
VIRAL_DATABASE_FILE="VDB.fa"
#
# ==============================================================================
#
# DEFAULT VALUES:
#
MINIMAL_SIMILARITY_VALUE=1.0;
TSIZE=2;
TOP_SIZE=0;
TOP_SIZE_VIR=0;
CACHE=70;
#
# ==============================================================================
# THESE ARE THE CURRENT FLAGGED VIRUSES OR VIRUSES GROUPS FOR ENHANCED ASSEMBLY:
# IF THIS LIST IS AUGMENTED THEN IT REQUIRES -> TRACES_get_best_<$pattern>.sh
# FILE TO SEARCH FOR SPECIFIC CHARACTERISTICS OF EACH VIRUS TYPE
# 
# HERV IS CURRENTLY BEING ADDRESSED AS HAPLOID IN ALIGNMENTS -> THIS WILL 
# REQUIRE ADAPTATION IN THE FUTURE.
#
declare -a VIRUSES=("B19" "HV1" "HV2" "HV3" "HV4" "HV5" "HV6" "HV6A" "HV6B" 
                    "HV7" "HV8" "POLY1" "POLY2" "POLY3" "POLY4" "POLY5" 
		    "POLY6" "POLY7" "POLY8" "POLY9" "POLY10" "POLY11" "POLY12" 
		    "POLY13" "POLY14" "HBV" "HPV" "TTV" "HBOV1" "HBOVNOT1" 
		    "VARV" "SV40" "CUTA" "HERV");
#
# ==============================================================================
# CHECK INTERNAL SYSTEM FILES
#
CHECK_FILTERING_SYSTEM_FILES () {
    if [ -n "$VIRAL_DATABASE_METADATA" ]; then
        if [ ! -f TRACES_get_best_by_meta.sh ]; then
            echo -e "\e[31mERROR: TRACES_get_best_by_meta.sh file not found!\e[0m"
            echo "This file may have been deleted acidentally."
            echo "System is corrupted!"
            exit 1;
        fi
    else
        for VIRUS in "${VIRUSES[@]}"; do
            if [ ! -f TRACES_get_best_$VIRUS.sh ]; then
                echo -e "\e[31mERROR: TRACES_get_best_$VIRUS.sh file not found!\e[0m"
                echo "This file may have been deleted acidentally or VIRAL array changed."
                echo "System is corrupted!"
                exit 1;
                fi
        done
    fi
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
  FILE_GZ1=`file -L ../input_data/$1 | grep gzip | wc -l`;
  FILE_GZ2=`file -L ../input_data/$2 | grep gzip | wc -l`;
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
  if [ ! -f "$VIRAL_DATABASE_FILE" ];
    then
    echo -e "\e[31mERROR: viral database ("$VIRAL_DATABASE_FILE") not found!\e[0m"
    if [ "$B_ALT_VIRAL_DB" -eq 0 ]; then
        echo "TIP: before this, run: ./TRACESPipe.sh --build-viral"
        echo " OR"
        echo "Use the --alt-viral-db option to specify an alternative VDB location"
        echo "For addition information, see the instructions at the web page."
    fi
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
PROGRAM_EXISTS () {
  printf "Checking $1 ... ";
  if ! [ -x "$(command -v $1)" ];
    then
    echo -e "\e[41mERROR\e[49m: $1 is not installed." >&2;
    echo -e "\e[42mTIP\e[49m: Try: ./TRACESPipe.sh --install" >&2;
    # exit 1;
    else
    echo -e "\e[42mSUCCESS!\e[49m";
    fi
  }
#
#
CHECK_PROGRAMS () {
  PROGRAM_EXISTS "trimmomatic";
  PROGRAM_EXISTS "cryfa";
  PROGRAM_EXISTS "MAGNET";
  PROGRAM_EXISTS "FALCON";
  PROGRAM_EXISTS "gto";
  PROGRAM_EXISTS "spades.py";
  PROGRAM_EXISTS "igv";
  PROGRAM_EXISTS "bowtie2";
  PROGRAM_EXISTS "samtools";
  PROGRAM_EXISTS "bcftools";
  PROGRAM_EXISTS "bedops";
  PROGRAM_EXISTS "bedtools";
  PROGRAM_EXISTS "fastq_pair";
  PROGRAM_EXISTS "efetch";
  PROGRAM_EXISTS "mapDamage";
  PROGRAM_EXISTS "tabix";
  PROGRAM_EXISTS "AdapterRemoval";
  PROGRAM_EXISTS "bwa";
  PROGRAM_EXISTS "art_illumina";
  PROGRAM_EXISTS "blastn";
  PROGRAM_EXISTS "dnadiff";
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
  ORGAN="$8";
  HIGH="$7";
  DUPL="$6";
  THREADS="$5";
  V_TAG="$1";
  IDX_TAG="$4";
  #
  V_GID="X";
  #
  echo -e "\e[34m[TRACESPipe]\e[32m Assessing $V_TAG best reference ...\e[0m";
  #
  CHECK_TOP "$ORGAN";
  #
  if [ -n "$VIRAL_DATABASE_METADATA" ]; then
      V_INFO="$(./TRACES_get_best_by_meta.sh "$ORGAN" "$V_TAG" "$VIRAL_DATABASE_METADATA")";
    else
        V_INFO=`./TRACES_get_best_$V_TAG.sh $ORGAN`;
    fi
  V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
  V_VAL=`echo "$V_INFO" | awk '{ print $1; }'`;
  # 
  echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
  #
  if [[ "$V_GID" != "-" ]] && [[ $(bc <<< "$V_VAL > $2") -eq 1  ]];
    then
    #
    if [[ "$3" -eq "1" ]]
      then	    
      echo -e "\e[34m[TRACESPipe]\e[96m Best reference: $V_GID\e[0m";
      V_GID=`sed "${IDX_TAG}q;d" ../output_data/TRACES_results/REPORT_META_VIRAL_BESTS_$ORGAN.txt | awk '{ print $2; }'`;
      echo -e "\e[34m[TRACESPipe]\e[96m Using best of bests reference: $V_GID\e[0m";
      fi
    #
    echo -e "\e[34m[TRACESPipe]\e[96m Similarity best match: $V_INFO\e[0m";
    echo -e "\e[34m[TRACESPipe]\e[32m Extracting sequence from "$VIRAL_DATABASE_FILE"\e[0m";
    CHECK_VDB;
    gto_fasta_extract_read_by_pattern -p "$V_GID" < "$VIRAL_DATABASE_FILE" | awk "/^>/ {n++} n>1 {exit} 1" > $ORGAN-$V_TAG.fa 2>> ../logs/Log-$ORGAN.txt;
    echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
    echo -e "\e[34m[TRACESPipe]\e[32m Aligning reads to $V_TAG best reference with bowtie2 ...\e[0m";
    ./TRACES_viral_align_reads.sh $ORGAN-$V_TAG.fa $ORGAN $V_TAG $THREADS $DUPL $HIGH 1>> ../logs/Log-stdout-$ORGAN.txt 2>> ../logs/Log-stderr-$ORGAN.txt;
    echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
    #
    echo -e "\e[34m[TRACESPipe]\e[32m Generate a consensus sequence with bcftools ...\e[0m";
    ./TRACES_viral_consensus.sh $ORGAN-$V_TAG.fa viral_aligned_sorted-$ORGAN-$V_TAG.bam $ORGAN $V_TAG 1>> ../logs/Log-stdout-$ORGAN.txt 2>> ../logs/Log-stderr-$ORGAN.txt;
    #
    mkdir -p ../output_data/TRACES_viral_alignments;
    cp $ORGAN-$V_TAG.fa ../output_data/TRACES_viral_alignments/
    cp $ORGAN-$V_TAG.fa.fai ../output_data/TRACES_viral_alignments/
    mv viral_aligned_sorted-$ORGAN-$V_TAG.bam ../output_data/TRACES_viral_alignments/
    mv viral_aligned_sorted-$ORGAN-$V_TAG.bam.bai ../output_data/TRACES_viral_alignments/
    mkdir -p ../output_data/TRACES_viral_consensus;
    #rm -f ../output_data/TRACES_viral_consensus/*
    mv $V_TAG-consensus-$ORGAN.fa ../output_data/TRACES_viral_consensus/
    mkdir -p ../output_data/TRACES_viral_bed;
    #rm -f ../output_data/TRACES_viral_bed/*
    mv $V_TAG-calls-$ORGAN.bed ../output_data/TRACES_viral_bed/
    mv $V_TAG-coverage-$ORGAN.bed ../output_data/TRACES_viral_bed/
    mv $V_TAG-zero-coverage-$ORGAN.bed ../output_data/TRACES_viral_bed/
    mv $V_TAG-$ORGAN-calls.vcf.gz ../output_data/TRACES_viral_bed/
    mkdir -p ../output_data/TRACES_viral_statistics;
    echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
    echo -e "\e[34m[TRACESPipe]\e[32m Calculating coverage ...\e[0m";
    ./TRACES_overall.sh viral $V_TAG $ORGAN
    C_BREADTH=`cat ../output_data/TRACES_viral_statistics/$V_TAG-total-horizontal-coverage-$ORGAN.txt`;
    C_DEPTH=`cat ../output_data/TRACES_viral_statistics/$V_TAG-total-depth-coverage-$ORGAN.txt`;
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
    -avdb|--alt-viral-db)
        B_ALT_VIRAL_DB=1;
        VIRAL_DATABASE_FILE="$2" 
        shift 2;
    ;;
    -flog|--flush-logs)
      FLUSH_LOGS=1;
      shift
    ;;
    -fout|--flush-output)
      FLUSH_OUTPUT=1;
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
    -spv|--show-prog-ver)
      SHOW_PROG_VER=1;
      SHOW_HELP=0;
      shift
    ;;
    -st|--sample)
      SAMPLE_TEST=1;
      SHOW_HELP=0;
      shift
    ;;
    -vdb|--build-viral)
      BUILD_VDB_ALL=1;
      SHOW_HELP=0;
      shift
    ;;
    -vdbm|--viral-db-metadata)
        VIRAL_DATABASE_METADATA="$2";
        if [ -s "$VIRAL_DATABASE_METADATA" ]; then
            #Parse the unique virus groups in the metadata (ignore lines without at least two columns
            readarray -t VIRUSES < <(awk -F '\\t' '
                (FNR > 1 && $1 && $2){VirusSet[$2]=1;}
                END {n=asorti(VirusSet,vl); for(i=1;i<=n;i++){print vl[i]}}
            ' "$VIRAL_DATABASE_METADATA")
            #Ensure at least one virus group is present
            if [ "${#VIRUSES[@]}" -lt 1 ]; then
                >&2 echo -e "\e[31mERROR: Provided Viral database Metadata file did not contain any valid accession-virus pairs\e[0m" 
                exit 1;
            fi
        else
            >&2 echo -e "\e[31mERROR: Provided Viral database Metadata file is non-existant or empty\e[0m" 
            exit 1;
        fi
        shift 2
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
    -sfrs|--search-blast-remote-db)
      SEARCH_BLAST_REMOTE_DB=1;
      BLAST_QUERY="$2";
      SHOW_HELP=0;
      shift 2
    ;;
    -rdup|--remove-dup)
      REMOVE_DUPLICATIONS=1;
      SHOW_HELP=0;
      shift
    ;;
    -gbb|--best-of-bests)
      RUN_BEST_OF_BESTS=1;
      SHOW_HELP=0;
      shift
    ;;
    -vhs|--very-sensitive)
      HIGH_SENSITIVITY=1;
      RUN_SPECIFIC_SENSITIVE=1;
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
    -dwms|--download-mito-species)
      DOWNLOAD_MITO_SPECIES=1;
      SHOW_HELP=0;
      shift; 
    ;;
    -dwmp|--download-mito-population)
      DOWNLOAD_MITO_POPULATION=1;
      SHOW_HELP=0;
      shift;
    ;;
    -aums|--auth-mito-species)
      RUN_MITO_SPECIES=1;
      SHOW_HELP=0;
      shift;
    ;;
    -aump|--auth-mito-population)
      RUN_MITO_POPULATION=1;
      SHOW_HELP=0;
      shift;
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
    -rpgi|--run-gid-complexity-profile)
      RUN_GID_COMPLEXITY_PROFILE=1;
      GID_COMPLEXITY_PROFILE="$2";
      SHOW_HELP=0;
      shift 2
    ;;
    -cpwi|--complexity-profile-window)
      COMPLEXITY_PROFILE_WINDOW="$2";
      SHOW_HELP=0;
      shift 2
    ;;
    -cple|--complexity-profile-level)
      COMPLEXITY_PROFILE_LEVEL="$2";
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
      RUN_HERV_ON=1;
      RUN_ALL_VIRAL=1;
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
    -rsd|--run-de-novo-specific)
      RUN_ANALYSIS=1;
      RUN_DENOVO_SPECIFIC=1;
      SPECIFIC_DENOVO_ID="$2";
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
    -clog|--coverage-log-scale)
      COVERAGE_LOG_SCALE="$2";
      SHOW_HELP=0;
      shift 2;
    ;;
    -cwis|--coverage-window-size)
      COVERAGE_WINDOW_SIZE="$2";
      SHOW_HELP=0;
      shift 2;
    ;;
    -cdro|--coverage-drop)
      COVERAGE_DROP="$2";
      SHOW_HELP=0;
      shift 2;
    ;;
    -covm|--coverage-min-x)
      COVERAGE_MIN_X="$2";
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
    -rgid|--run-gid-damage)
      RUN_GID_DAMAGE_ANALYSIS=1;
      GID_DAMAGE="$2";
      SHOW_HELP=0;
      shift 2;
    ;;
    -rava|--run-all-v-alig)
      RUN_ANALYSIS=1;
      RUN_ALL_VIRAL=1;
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
    -rmhc|--run-multiorgan-consensus)
      RUN_ANALYSIS=1;
      RUN_MULTIORGAN_CONSENSUS=1;
      SHOW_HELP=0;
      shift
    ;;
    -diff|--run-diff)
      RUN_DIFF=1;
      RUN_ANALYSIS=1;
      SHOW_HELP=0;
      shift
    ;;
    -sdiff|--run-specific-diff)
      RUN_SPECIFIC_DIFF=1;
      RUN_DIFF_VIRUS="$2";
      RUN_DIFF_ID="$3";
      SHOW_HELP=0;
      shift 3
    ;;
    -brec|--blast-reconstructed)
      RUN_BLAST_RECONSTRUCTED=1;
      RUN_ANALYSIS=1;
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
      RUN_MITO_ON=1;
      RUN_MITO_DAMAGE_ON=1;
      RUN_ALL_VIRAL=1;
      RUN_CY_ON=1;
      RUN_CY_QUANT_ON=1;
      RUN_DE_NOVO_ASSEMBLY=1;
      RUN_HYBRID=1;
      RUN_MULTIORGAN_CONSENSUS=1;
      RUN_DIFF=1;
      SHOW_HELP=0;
      shift
    ;;
    -c|--cache)
      CACHE="$2";
      SHOW_HELP=0;
      shift 2;
    ;;
    -ts|--top-size)
      TOP_SIZE="$2";
      SHOW_HELP=0;
      shift 2;
    ;;    
    -tsv|--top-size-virus)
      TOP_SIZE_VIR="$2";
      SHOW_HELP=0;
      shift 2;
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
  echo "                                                                       "
  echo -e "\e[34m                                                              "
  echo "         ████████╗ ██████╗   █████╗   ██████╗ ███████╗ ███████╗        "
  echo "         ╚══██╔══╝ ██╔══██╗ ██╔══██╗ ██╔════╝ ██╔════╝ ██╔════╝        "
  echo "            ██║    ██████╔╝ ███████║ ██║      █████╗   ███████╗        "
  echo "            ██║    ██╔══██╗ ██╔══██║ ██║      ██╔══╝   ╚════██║        "
  echo "            ██║    ██║  ██║ ██║  ██║ ╚██████╗ ███████╗ ███████║        "
  echo "            ╚═╝    ╚═╝  ╚═╝ ╚═╝  ╚═╝  ╚═════╝ ╚══════╝ ╚══════╝        "
  echo "                                                                       "
  echo "                             P I P E L I N E                           "
  echo "                                                                       "
  echo -e "    \e[32m      |  A hybrid pipeline for reconstruction & analysis  | \e[0m" 
  echo -e "    \e[32m      |  of viral and host genomes at multi-organ level.  | \e[0m"
  echo "                                                                "
  echo -e "\e[93m    Usage: ./TRACESPipe.sh [options]                     \e[0m"
  echo "                                                                       "
  echo -e "\e[34m ===========            GENERAL OPTIONS              ==========\e[0m"
  echo "                                                                       "
  echo "    -h,     --help            Show this help message and exit,         "
  echo "    -v,     --version         Show the version and some information,   "
  echo "    -f,     --force           Force running and overwrite of files,    "
  echo "    -flog,  --flush-logs      Flush logs (delete logs),                "
  echo "    -fout,  --flush-output    Flush output data (delete all output_data), "
  echo "    -t <THREADS>, --threads <THREADS> Number of threads to use,        "
  echo "                                                                       "
  echo -e "\e[34m ===========            SETUP COMMANDS               ==========\e[0m"
  echo "                                                                       "
  echo "    -i,     --install         Installation of all the tools,           "
  echo "    -up,    --update          Update all the tools in TRACESPipe,      "
  echo "    -spv,   --show-prog-ver   Show included programs versions,         "
  echo "                                                                       "
  echo "    -st,    --sample          Creates human ref. VDB and sample organ, "
  echo "                                                                       "
  echo "    -gmt,   --get-max-threads Get the number of maximum machine threads,"
  echo "                                                                       "
  echo "    -dec,   --decrypt         Decrypt (all files in ../encrypted_data), "
  echo "    -enc,   --encrypt         Encrypt (all files in ../to_encrypt_data),"
  echo "                                                                       "
  echo "    -vdb,   --build-viral     Build viral database (all) [Recommended], "
  echo "    -vdbr,  --build-viral-r   Build viral database (references only),  "
  echo "    -udb,   --build-unviral   Build non viral database (control),      "
  echo "                                                                       "
  echo "    -afs <FASTA>, --add-fasta <FASTA>                                  "
  echo "                              Add a FASTA sequence to the "$VIRAL_DATABASE_FILE",      "
  echo "    -aes <ID>, --add-extra-seq <ID>                                    "
  echo "                              Add extra sequence to the "$VIRAL_DATABASE_FILE",        "
  echo "    -gx,    --get-extra-vir   Downloads/appends (VDB) extra viral seq, "
  echo "                                                                       "
  echo "    -gad,   --gen-adapters    Generate FASTA file with adapters,       "
  echo "    -gp,    --get-phix        Extracts PhiX genomes (Needs viral DB),  "
  echo "    -gm,    --get-mito        Downloads human Mitochondrial genome,    "
  echo "                                                                       "
  echo "    -dwms,  --download-mito-species                                    "
  echo "                              Downloads the complete NCBI mitogenomes  "
  echo "                              database containing the existing species,"
  echo "                                                                       "
  echo "    -dwmp,  --download-mito-population                                 "
  echo "                              Downloads two complete mitogenome databases "
  echo "                              with healthy and pathogenic sequences,   "
  echo "                                                                       "
  echo "    -aums,  --auth-mito-species                                        "
  echo "                              Autheticate the mitogenome species,      "
  echo "                                                                       "
  echo "    -aump,  --auth-mito-population                                     "
  echo "                              Authenticate closest population,         "
  echo "                                                                       "
  echo "    -cmt <ID>, --change-mito <ID>                                      "
  echo "                              Set any Mitochondrial genome by ID,      "
  echo "                                                                       "
  echo "    -gy,    --get-y-chromo    Downloads human Y-chromosome,            "
  echo "    -gax,   --get-all-aux     Runs -gad -gp -gm -gy,                   "
  echo "                                                                       "
  echo "    -cbn,   --create-blast-db It creates a nucleotide blast database,  "
  echo "    -ubn,   --update-blast-db It updates a nucleotide blast database,  "
  echo "                                                                       "
  echo -e "\e[34m ===========           ANALYSIS COMMANDS             ==========\e[0m"
  echo "                                                                       "
  echo "    -ra,    --run-analysis    Run data analysis (core),                      "
  echo "    -all,   --run-all         Run all the options (excluding the specific).  "
  echo "                                                                       "
  echo "    -sfs <FASTA>, --search-blast-db <FASTA>                            "
  echo "                              It blasts the nucleotide (nt) blast DB,  "
  echo "    -sfrs <FASTA>, --search-blast-remote-db <FASTA>                    "
  echo "                              It blasts remotly thenucleotide (nt) blast "
  echo "                              database (it requires internet connection), "
  echo "                                                                       "
  echo "    -gbb,   --best-of-bests   Identifies the best of bests references  "
  echo "                              between multiple organs [similar reference], "
  echo "                                                                       "
  echo "    -rpro,  --run-profiles    Run complexity and relative profiles (control), "
  echo "                                                                       "
  echo "    -rpgi <ID>,  --run-gid-complexity-profile <ID>                     "
  echo "                              Run complexity profiles by GID,          "
  echo "    -rm,    --run-meta        Run viral metagenomic identification,    "
  echo "    -ro,    --run-meta-nv     Run NON-viral metagenomic identification,"
  echo "    -rava,  --run-all-v-alig  Run all viral align/sort/consensus seqs  "
  echo "                              from a specific list,                    "
  echo "                                                                       "
  echo "    -rsd <ID>, --run-de-novo-specific <ID/PATTERN>                     "
  echo "                              Run specific alignments of the de-novo   "
  echo "                              to the reference genome,                 "
  echo "    -rsr <ID>, --run-specific <ID/PATTERN>                             "
  echo "                              Run specific reference align/consensus,  "
  echo "                                                                       "
  echo "    -rsx <ID>, --run-extreme <ID/PATTERN>                              "
  echo "                              Run specific reference align/consensus   "
  echo "                              using extreme sensitivity;               "
  echo "                              Retained for backwards compatibility;    "
  echo "                              Now an alias for -vhs -rsr <ID/PATTERN>, "
  echo "                                                                       "
  echo "    -rmt,   --run-mito        Run Mito align and consensus seq,        "
  echo "    -rmtd,  --run-mito-dam    Run Mito damage only,                    "
  echo "                                                                       "
  echo "    -rgid <ID>, --run-gid-damage <ID>                                  "
  echo "                              Run damage pattern analysis by GID,      "
  echo "                                                                       "
  echo "    -rya,   --run-cy-align    Run CY align and consensus seq,          "
  echo "    -ryq,   --run-cy-quant    Estimate the quantity of CY DNA,         "
  echo "                                                                       "
  echo "    -rda,   --run-de-novo     Run de-novo assembly,                    "
  echo "                                                                       "
  echo "    -rhyb,  --run-hybrid      Run hybrid assembly (align/de-novo),     "
  echo "                                                                       "
  echo "    -rmhc,  --run-multiorgan-consensus                                 "
  echo "                              Run alignments/consensus between all the "
  echo "                              reconstructed organ sequences,           "
  echo "                                                                       "
  echo "    -vis,   --visual-align    Run Visualization tool for alignments,   "
  echo "    -covl,  --coverage-latex  Run coverage table in Latex format,      "
  echo "    -covc,  --coverage-csv    Run coverage table in CSV format,        "
  echo "                                                                       "
  echo "    -covp <NAME>,  --coverage-profile <BED_NAME_FILE>                   "
  echo "                              Run coverage profile for specific BED file, "
  echo "                                                                       "
  echo "    -diff,  --run-diff        Run diff -> reference and hybrid (ident/SNPs), "
  echo "                                                                       "
  echo "    -sdiff <V_NAME> <ID/PATTERN>, --run-specific-diff <V_NAME> <ID/PATTERN>  "
  echo "                              Run specific diff of reconstructed to a virus  "
  echo "                              pattern of ID. Example: -sdiff B19 AY386330.1, "
  echo "                                                                       "
  echo "    -brec,  --blast-reconstructed                                      "
  echo "                              Run local blast over reconstructed genomes, "
  echo "                                                                             "
  echo "                                                                       "
  echo -e "\e[34m ===========           ANALYSIS OPTIONS              ==========\e[0m"
  echo "                                                                       "
  echo "    -avdb <FASTA>, --alt-viral-db <FASTA>                              "
  echo "                              Specify a path to fasta file containing  "
  echo "                               viral sequences                         "
  echo "                              Sequence names must include the accession"
  echo "                               as the first field, either whitespace or"
  echo "                               underscore (_) delimited (NC_ handled)  "
  echo "    -vdbm <PATH>, --viral-db-metadata <PATH>                           "
  echo "                              Specify a path to a tab del file which   "
  echo "                               has sequence GID/ACC in the first column"
  echo "                               and a Name representing a virus label   "
  echo "                               in the second                           "
  echo "                              This changes the meta-analysis behaviour:"
  echo "                               rather than using internal virus names  "
  echo "                               and inclusion criteria, the relationships"
  echo "                               defined in the provided file are used.  "
  echo "                              This allows the user to define groupings "
  echo "                               of interest in the default, or user     "
  echo "                               provided viral database.                "
  echo "                                                                       "
  echo "    -rdup,  --remove-dup      Remove duplications (e.g. PCR dup),      "
  echo "    -vhs,   --very-sensitive  Aligns with very high sensitivity (slower),  "
  echo "                                                                       "
  echo "    -iss <SIZE>, --inter-sim-size <SIZE>                               "
  echo "                              Inter-genome similarity top size (control), "
  echo "                                                                       "
  echo "    -cpwi <VALUE>, --complexity-profile-window <VALUE>                 "
  echo "                              Complexity profile window size,          "
  echo "    -cple <VALUE>, --complexity-profile-level <VALUE>                  "
  echo "                              Complexity profile compression level [1;10], "
  echo "                                                                       "
  echo "    -mis <VALUE>, --min-similarity <VALUE>                             "
  echo "                              Minimum similarity value to consider the "
  echo "                              sequence for alignment-consensus (filter), "
  echo "                                                                       "
  echo "    -top <VALUE>, --view-top <VALUE>                                   "
  echo "                              Display the top <VALUE> with the highest "
  echo "                              similarity (by descending order),        "
  echo "                                                                       "
  echo "    -c <VALUE>,   --cache <VALUE>                                      "
  echo "                              Cache to be used by FALCON-meta,         "
  echo "    -tsv <VALUE>,   --top-size-virus <VALUE>                           "
  echo "                              Top size to be used by FALCON-meta when  "
  echo "                              using TRACES_metagenomic_viral.sh;       "
  echo "                              default:0 -> seq count in viral db       "
  echo "    -ts <VALUE>,   --top-size <VALUE>                                  "
  echo "                              Top size to be used by FALCON-meta when  "
  echo "                              using TRACES_metagenomic.sh;             "
  echo "                              default:0 -> seq count in  non-viral db  "
  echo "                                                                       "
  echo "    -cmax <MAX>,   --max-coverage <MAX_COVERAGE>                       "
  echo "                              Maximum depth coverage (depth normalization), "
  echo "    -clog <VALUE>, --coverage-log-scale <VALUE>                        "
  echo "                              Coverage profile logarithmic scale VALUE=Base, "
  echo "    -cwis <VALUE>, --coverage-window-size <VALUE>                      "
  echo "                              Coverage window size for low-pass filter, "
  echo "    -cdro <VALUE>, --coverage-drop <VALUE>                             "
  echo "                              Coverage drop size (sampling),           "
  echo "    -covm <VALUE>, --coverage-min-x <VALUE>                             "
  echo "                              Coverage minimum value for x-axis        "
  echo "                                                                       "
  echo -e "\e[34m ===========                EXAMPLES                 ==========\e[0m"
  echo "                                                                       "
  echo -e "\e[93m    Ex: ./TRACESPipe.sh --flush-output --flush-logs --run-mito --run-meta \e[0m"
  echo -e "\e[93m    --remove-dup --run-de-novo --run-hybrid --min-similarity 1 --run-diff \e[0m" 
  echo -e "\e[93m    --very-sensitive --best-of-bests --run-multiorgan-consensus \e[0m"
  echo "                                                                       "
  echo "    Add the file meta_info.txt at ../meta_data/ folder. Example:       "
  echo "    meta_info.txt -> 'organ:reads_forward.fa.gz:reads_reverse.fa.gz'   "
  echo "    The reads must be GZIPed in the ../input_data/ folder.             "
  echo "    The output results are at ../output_data/ folder.                  "
  echo "                                                                       "
  echo -e "\e[32m    Contact: tracespipe@gmail.com                        \e[0m"
  echo "                                                                       "
  exit 1;
  fi
#
# ==============================================================================
# VERSION
#
if [ "$SHOW_VERSION" -eq "1" ]; then
    version="1.1.3";
    versionFile="$SOURCE_DIR/../Version.txt"
    [ -f "$versionFile" ] && version=$(cat "$versionFile")
  echo "                                                                      ";
  echo "                              TRACESPipe                              ";
  echo "                                                                      ";
  echo "                            Version: $version                            ";
  echo "                                                                      ";
  echo "                      Department of Virology and                      ";
  echo "                   Department of Forensic Medicine,                   ";
  echo "                   University of Helsinki, Finland.                   ";
  echo "                                  &                                   ";
  echo "                             IEETA/DETI,                              ";
  echo "                    University of Aveiro, Portugal.                   ";
  echo "                                                                      ";
  echo "                       diogo.pratas@helsinki.fi                       ";
  echo "                      zachery.dickson@helsinki.fi                     ";
  echo "                                                                      ";
  exit 0;
  fi
#
# ==============================================================================
#
# MAKE SURE FOLDERS EXIST
#
mkdir -p ../logs/
mkdir -p ../output_data/
#
# MAKE SURE PROGRAMS EXIST
#
CHECK_PROGRAMS;
#
# DELETE OLD RUN FILES
#
if [[ "$RUN_DIFF" -eq "1" ]];
  then
  rm -f ../output_data/TRACES_diff/Viral_Diff.txt
  rm -f ../output_data/TRACES_diff/Viral_Diff_after_blastn.txt
  rm -f ../output_data/TRACES_diff/mtDNA_Diff.txt
  fi
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
if [[ "$FLUSH_OUTPUT" -eq "1" ]];
  then
  echo -e "\e[34m[TRACESPipe]\e[32m Flushing output data ...\e[0m";
  rm -fr ../output_data/*
  echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
  fi
#
# ==============================================================================
#
if [[ "$SHOW_PROG_VER" -eq "1" ]];
  then
  ./TRACES_get_program_versions.sh
  exit 0;
  fi
#
# ==============================================================================
#
if [[ "$SAMPLE_TEST" -eq "1" ]];
  then
  ./TRACES_simulate_simple_sample.sh
  exit 0;
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
if [[ "$DOWNLOAD_MITO_SPECIES" -eq "1" ]];
  then
  echo -e "\e[34m[TRACESPipe]\e[32m Downloading species mitogenome database ...\e[0m";
  ./TRACES_download_base_species.sh
  echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
  fi
#
# ==============================================================================
#
if [[ "$DOWNLOAD_MITO_POPULATION" -eq "1" ]];
  then
  echo -e "\e[34m[TRACESPipe]\e[32m Downloading two mitogenomes population databases ...\e[0m";
  ./TRACES_download_base_population.sh
  echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
  fi
#
# ==============================================================================
#
if [[ "$CREATE_BLAST_DB" -eq "1" ]];
  then
  echo -e "\e[34m[TRACESPipe]\e[32m Creating Blast nt database ...\e[0m";
  echo -e "\e[34m[TRACESPipe]\e[32m (This may take a while...)\e[0m";
  ./TRACES_blastn_create_n_db.sh 1>> ../logs/Log-stdout-system.txt 2>> ../logs/Log-stderr-system.txt;
  echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
  exit 0;
  fi
#
# ==============================================================================
#
if [[ "$UPDATE_BLAST_DB" -eq "1" ]];
  then
  echo -e "\e[34m[TRACESPipe]\e[32m Updating Blast nt database ...\e[0m";
  echo -e "\e[34m[TRACESPipe]\e[32m (This may take a while...)\e[0m";
  ./TRACES_blastn_update_n_db.sh 1>> ../logs/Log-stdout-system.txt 2>> ../logs/Log-stderr-system.txt;
  echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
  exit 0;
  fi
#
# ==============================================================================
#
if [[ "$SEARCH_BLAST_DB" -eq "1" ]];
  then
  echo -e "\e[34m[TRACESPipe]\e[32m Searching $BLAST_QUERY in local Blast nt database ...\e[0m";
  BLASTS_OUTPUT="../output_data/TRACES_blasts/";
  mkdir -p $BLASTS_OUTPUT;
  ./TRACES_blastn_n_db.sh $BLAST_QUERY > $BLASTS_OUTPUT/$BLAST_QUERY.txt 2>> ../logs/Log-stderr-system.txt;
  head -n 10 $BLASTS_OUTPUT/$BLAST_QUERY.txt;
  echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
  exit 0;
  fi
#
# ==============================================================================
#
if [[ "$SEARCH_BLAST_REMOTE_DB" -eq "1" ]];
  then
  echo -e "\e[34m[TRACESPipe]\e[32m Searching $BLAST_QUERY in remote Blast nt database ...\e[0m";
  BLASTS_OUTPUT="../output_data/TRACES_blasts/";
  mkdir -p $BLASTS_OUTPUT;
  blastn -db nt -task blastn-short -query $BLAST_QUERY -remote > $BLASTS_OUTPUT/x.txt 2>> ../logs/Log-stderr-system.txt;
  echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
  exit 0;
  fi
#
# ==============================================================================
#
if [[ "$RUN_GID_COMPLEXITY_PROFILE" -eq "1" ]];
  then
  efetch -db nucleotide -format fasta -id "$GID_COMPLEXITY_PROFILE" > $GID_COMPLEXITY_PROFILE.fa 
  CHECK_E_FILE $GID_COMPLEXITY_PROFILE.fa
  #
  echo -e "\e[34m[TRACESPipe]\e[32m Building complexity profiles for $GID_COMPLEXITY_PROFILE ...\e[0m";
  ./TRACES_complexity_profile.sh $GID_COMPLEXITY_PROFILE.fa $COMPLEXITY_PROFILE_LEVEL $COMPLEXITY_PROFILE_WINDOW 1>> ../logs/Log-stdout.txt 2>> ../logs/Log-stderr.txt;
  mkdir -p ../output_data/TRACES_results/profiles/
  cp complexity-profile-$GID_COMPLEXITY_PROFILE.fa.pdf ../output_data/TRACES_results/profiles/
  echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
  #
  fi
#
# ==============================================================================
#
if [[ "$RUN_COVERAGE_PROFILE" -eq "1" ]];
  then
  #
  CHECK_E_FILE $COVERAGE_NAME
  #
  rm -f x.projected.profile;
  ./TRACES_project_coordinates.sh $COVERAGE_NAME $COVERAGE_MAX | gto_filter -w $COVERAGE_WINDOW_SIZE -d $COVERAGE_DROP > x.projected.profile;
  #
  if [[ "$COVERAGE_LOG_SCALE" -eq "" ]];
    then
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
    set yrange [$COVERAGE_MIN_X:]
    set xrange [:]
    set xtics auto
    set grid
    set ylabel "Depth"
    set xlabel "Position"
    set border linewidth 1.5
    set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 5 ps 0.4 # --- blue
    set style line 2 lc rgb '#000000' lt 1 lw 2 pt 6 ps 0.4 # --- green
    set style line 3 lc rgb '#dd181f' lt 1 lw 4 pt 7 ps 0.4 # --- ?
    set style line 4 lc rgb '#4d1811' lt 1 lw 4 pt 8 ps 0.4 # --- ?
    set style line 5 lc rgb '#1d121f' lt 1 lw 4 pt 9 ps 0.4 # --- ?
    plot "x.projected.profile" using 1:2 t "Depth coverage" with lines ls 2
EOF
    else
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
    set yrange [$COVERAGE_MIN_X:]
    set xrange [:]
    set xtics auto
    set logscale y $COVERAGE_LOG_SCALE
    set grid
    set ylabel "Depth"
    set xlabel "Position"
    set border linewidth 1.5
    set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 5 ps 0.4 # --- blue
    set style line 2 lc rgb '#000000' lt 1 lw 2 pt 6 ps 0.4 # --- green
    set style line 3 lc rgb '#dd181f' lt 1 lw 4 pt 7 ps 0.4 # --- ?
    set style line 4 lc rgb '#4d1811' lt 1 lw 4 pt 8 ps 0.4 # --- ?
    set style line 5 lc rgb '#1d121f' lt 1 lw 4 pt 9 ps 0.4 # --- ?
    plot "x.projected.profile" using 1:2 t "Depth coverage" with lines ls 2
EOF
    fi
  rm -f x.projected.profile;
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
if [[ "$HIGH_SENSITIVITY" -eq "1" ]];
  then
  echo -e "\e[34m[TRACESPipe]\e[32m Running with High sensitivity.\e[0m";
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
    if [ "$B_ALT_VIRAL_DB" -eq 1 ]; then
       >&2 echo -e "\e[31mERROR: Cannot add fasta when using an alternative viral DB\e[0m" 
       exit 1
    fi
  CHECK_VDB;
  ./TRACES_add_fasta_to_database.sh "$VIRAL_DATABASE_FILE" $NEW_FASTA
  fi
#
# ==============================================================================
#
if [[ "$GET_EXTRA" -eq "1" ]];
  then
  if [ "$B_ALT_VIRAL_DB" -eq 1 ]; then
     >&2 echo -e "\e[31mERROR: Cannot add extra viral seq when using an alternative viral DB\e[0m" 
     exit 1
  fi
  CHECK_ENRICH;
  ./TRACES_get_enriched_sequences.sh "$VIRAL_DATABASE_FILE"
  fi
#
# ==============================================================================
#
if [[ "$ADD_EXTRA_SEQ" -eq "1" ]];
  then
    if [ "$B_ALT_VIRAL_DB" -eq 1 ]; then
       >&2 echo -e "\e[31mERROR: Cannot add extra seq when using an alternative viral DB\e[0m" 
       exit 1
    fi
  CHECK_VDB;
  ./TRACES_get_extra_seq.sh "$VIRAL_DATABASE_FILE" $NEW_SEQ_ID
  fi
#
#
# ==============================================================================
#
if [[ "$RUN_DIFF" -eq "1" ]];
   then
   mkdir -p ../output_data/TRACES_diff
   printf "Organ\tVirus\tID\tSimilarity\tAligned bases\tidentity\tSNPs\tBreadth\tDepth\n" > ../output_data/TRACES_diff/Viral_Diff.txt
   printf "Organ\tID\tSimilarity\tAligned bases\tidentity\tSNPs\tBreadth\tDepth\n" > ../output_data/TRACES_diff/mtDNA_Diff.txt
   fi
#
# ==============================================================================
#
if [[ "$RUN_SPECIFIC_DIFF" -eq "1" ]];
   then
   #
   echo -e "\e[34m[TRACESPipe]\e[32m Running specific diff [VIRUS: $RUN_DIFF_VIRUS ; ID: $RUN_DIFF_ID] ...\e[0m";
   if [[ ! "${VIRUSES[@]}" =~ "$RUN_DIFF_VIRUS" ]]; 
     then
     echo -e "\e[34m[TRACESPipe]\e[31m ERROR: virus label does not exist!\e[0m";
     echo -e "\e[34m[TRACESPipe]\e[33m TIP -> existing labels: ${VIRUSES[*]}\e[0m";
     exit 1;
     fi
   #
   ./TRACES_run_specific_diff.sh $RUN_DIFF_VIRUS $RUN_DIFF_ID $THREADS ${VIRUSES[@]}
   echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
   fi
#
# ==============================================================================
#
if [[ "$RUN_MITO_POPULATION" -eq "1" ]];
  then
  #
  echo -e "\e[34m[TRACESPipe]\e[32m Running population mitogenome authentication ...\e[0m";
  ./TRACES_auth_population.sh $THREADS
  echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
  #
  fi
#
# ==============================================================================
#
if [[ "$RUN_MITO_SPECIES" -eq "1" ]];
  then
  #       
  echo -e "\e[34m[TRACESPipe]\e[32m Running species mitogenome authentication ...\e[0m";  
  ./TRACES_auth_species.sh $THREADS
  echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
  #
  fi
#
# ==============================================================================
#
if [[ "$RUN_GID_DAMAGE_ANALYSIS" -eq "1" ]];
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
    echo -e "\e[34m[TRACESPipe]\e[93m Running organ=$ORGAN_T\e[0m";
    #
    rm -f FW_READS.fq.gz RV_READS.fq.gz
    echo -e "\e[34m[TRACESPipe]\e[32m Copying an instance of the files ...\e[0m";
    cp ../input_data/$SPL_Forward FW_READS.fq.gz;
    cp ../input_data/$SPL_Reverse RV_READS.fq.gz;
    echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
    #
    ./TRACES_gid_damage_patterns.sh $GID_DAMAGE $ORGAN_T $THREADS
    #
    done
  #
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
  # ============================================================================
  # 
  if [[ "$RUN_BEST_OF_BESTS" -eq "1" ]];
    then
    #
    echo -e "\e[34m[TRACESPipe]\e[36m Running best of bests ...\e[0m";
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
      echo -e "\e[34m[TRACESPipe]\e[93m Running organ=$ORGAN_T\e[0m";
      #
      rm -f FW_READS.fq.gz RV_READS.fq.gz
      echo -e "\e[34m[TRACESPipe]\e[32m Copying an instance of the files ...\e[0m";
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
      if [[ "$RUN_META_ON" -eq "1" ]];
        then
	#
        CHECK_VDB;
        CHECK_PHIX;
        [ $TOP_SIZE_VIR -eq 0 ] && TOP_SIZE_VIR=$(grep -c '^>' "$VIRAL_DATABASE_FILE");
        #
        mv o_fw_pr.fq NP-o_fw_pr.fq;
        mv o_fw_unpr.fq NP-o_fw_unpr.fq;
        mv o_rv_pr.fq NP-o_rv_pr.fq;
        mv o_rv_unpr.fq NP-o_rv_unpr.fq;
        #
        echo -e "\e[34m[TRACESPipe]\e[32m Running viral metagenomic analysis with FALCON-meta ...\e[0m";
        mkdir -p ../output_data/TRACES_results
        ./TRACES_metagenomics_viral.sh $ORGAN_T "$VIRAL_DATABASE_FILE" $TOP_SIZE_VIR $THREADS $TSIZE $CACHE 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
        echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
        #
        echo -e "\e[34m[TRACESPipe]\e[32m Finding the best references ...\e[0m";
        #
        cp top-$ORGAN_T.csv ../output_data/TRACES_results/
        #
        rm -f ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt;
        for VIRUS in "${VIRUSES[@]}"
          do
            if [ -n "$VIRAL_DATABASE_METADATA" ]; then
                ./TRACES_get_best_by_meta.sh "$ORGAN_T" "$VIRUS" "$VIRAL_DATABASE_METADATA"
            else
                ./TRACES_get_best_$VIRUS.sh $ORGAN_T 
            fi >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
          done
        echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
        #
        fi
      #	
      done
    #
    IDX=1;
    rm -f TMP_VIRAL_GENERAL.txt;
    for VIRUS in "${VIRUSES[@]}"
      do
      #
      rm -f V_F_STRINGS;
      for read in "${READS[@]}";  
        do
        ORGAN_T=`echo $read | tr ':' '\t' | awk '{ print $1 }'`;
        sed "${IDX}q;d" ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt | awk '{ print $2;}' >> V_F_STRINGS;
        done
      ((++IDX));
      #
      BEST_OF_BESTS=`grep -v "-" V_F_STRINGS | awk '{++a[$0]}END{for(i in a)if(a[i]>max){max=a[i];k=i}print k}'`;
      if [[ -z "$BEST_OF_BESTS" ]];
        then
        printf "%s\t%s\n" "-" "-" >> TMP_VIRAL_GENERAL.txt;
        else
        printf "Best\t$BEST_OF_BESTS\n" >> TMP_VIRAL_GENERAL.txt;
      fi
      #
      done
    #
    for read in "${READS[@]}";
      do
      ORGAN_T=`echo $read | tr ':' '\t' | awk '{ print $1 }'`;
      cp TMP_VIRAL_GENERAL.txt ../output_data/TRACES_results/REPORT_META_VIRAL_BESTS_$ORGAN_T.txt
      done
    #
    fi
  #  
  # ============================================================================
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
    echo -e "\e[34m[TRACESPipe]\e[32m Copying an instance of the files ...\e[0m";
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
    if [[ "$RUN_META_ON" -eq "1" && "$RUN_BEST_OF_BESTS" -ne "1" ]];
      then
      CHECK_VDB;
      CHECK_PHIX;
      [ $TOP_SIZE_VIR -eq 0 ] && TOP_SIZE_VIR=$(grep -c '^>' "$VIRAL_DATABASE_FILE");
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Removing PhiX from the samples with MAGNET ...\e[0m";
      #./TRACES_remove_phix.sh $THREADS 
      # XXX: REMOVE FROM THE DB PHIX OR CHANGE MAGNET FOR PAIRED-END READS 
      # XXX: ALTERNATIVE - USE BWA OR OTHER ALIGNER TO EXTRACT THOSE READS 
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
      ./TRACES_metagenomics_viral.sh $ORGAN_T "$VIRAL_DATABASE_FILE" $TOP_SIZE_VIR $THREADS $TSIZE $CACHE 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Finding the best references ...\e[0m";
      #
      cp top-$ORGAN_T.csv ../output_data/TRACES_results/
      #
      rm -f ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt;
      for VIRUS in "${VIRUSES[@]}"
        do
        if [ -n "$VIRAL_DATABASE_METADATA" ]; then
                ./TRACES_get_best_by_meta.sh "$ORGAN_T" "$VIRUS" "$VIRAL_DATABASE_METADATA"
            else
                ./TRACES_get_best_$VIRUS.sh $ORGAN_T 
            fi >> ../output_data/TRACES_results/REPORT_META_VIRAL_$ORGAN_T.txt
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
      ./TRACES_profiles.sh GIS-$ORGAN_T "$VIRAL_DATABASE_FILE" P_TRACES_sample_reads.fq $ORGAN_T 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
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
      [ $TOP_SIZE -eq 0 ] && TOP_SIZE=$(grep -c '^>' DB.fa);
      #	
      echo -e "\e[34m[TRACESPipe]\e[32m Running NON viral metagenomic analysis with FALCON ...\e[0m";
      ./TRACES_metagenomics.sh $ORGAN_T DB.fa $TOP_SIZE $THREADS $TSIZE $CACHE 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
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
	echo -e "\e[34m[TRACESPipe]\e[32m Possibly computer ran out of RAM!\e[0m";
	else
        head -n $VTOP_SIZE ../output_data/TRACES_results/top-$ORGAN_T.csv
        fi
      fi
    #
    # ==========================================================================
    # RUN SPECIFIC ALIGN/CONSENSUS
    #
    if [[ "$RUN_SPECIFIC" -eq "1" ]];
      then
      echo -e "\e[34m[TRACESPipe]\e[32m Aligning reads to specific viral ref(s) with pattern \"$SPECIFIC_ID\" using bowtie2 ...\e[0m";
      #
      CHECK_VDB;
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Extracting sequence with pattern \"$SPECIFIC_ID\" from "$VIRAL_DATABASE_FILE" ...\e[0m";
      gto_fasta_extract_read_by_pattern -p "$SPECIFIC_ID" < "$VIRAL_DATABASE_FILE" | awk "/^>/ {n++} n>1 {exit} 1" > SPECIFIC-$SPECIFIC_ID.fa
      echo -e "\e[34m[TRACESPipe]\e[32m Aligning ... \e[0m";
      ./TRACES_viral_align_reads.sh SPECIFIC-$SPECIFIC_ID.fa $ORGAN_T $SPECIFIC_ID $THREADS $REMOVE_DUPLICATIONS $RUN_SPECIFIC_SENSITIVE 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
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
      mv $SPECIFIC_ID-consensus-$ORGAN_T.fa ../output_data/TRACES_specific_consensus/
      mkdir -p ../output_data/TRACES_specific_bed;
      mv $SPECIFIC_ID-calls-$ORGAN_T.bed ../output_data/TRACES_specific_bed/
      mv $SPECIFIC_ID-coverage-$ORGAN_T.bed ../output_data/TRACES_specific_bed/
      mv $SPECIFIC_ID-zero-coverage-$ORGAN_T.bed ../output_data/TRACES_specific_bed/
      mv $SPECIFIC_ID-$ORGAN_T-calls.vcf.gz ../output_data/TRACES_specific_bed/
      mkdir -p ../output_data/TRACES_specific_statistics;
      echo -e "\e[34m[TRACESPipe]\e[32m Calculating coverage for alignment-based assembly ...\e[0m";
      ./TRACES_overall.sh specific $SPECIFIC_ID $ORGAN_T
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
    if [[ "$RUN_ALL_VIRAL" -eq "1" ]];
      then
      IDX=1;
      for VIRUS in "${VIRUSES[@]}"
        do
        ALIGN_AND_CONSENSUS "$VIRUS" "$MINIMAL_SIMILARITY_VALUE" "$RUN_BEST_OF_BESTS" "$IDX" "$THREADS" "$REMOVE_DUPLICATIONS" "$HIGH_SENSITIVITY" "$ORGAN_T"
	((++IDX));
	done
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
      ./TRACES_mt_align_reads.sh mtDNA.fa $ORGAN_T $THREADS $REMOVE_DUPLICATIONS $HIGH_SENSITIVITY 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
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
      echo -e "\e[34m[TRACESPipe]\e[32m Calculating coverage for alignment-based assembly ...\e[0m";
      ./TRACES_overall.sh mtdna mt $ORGAN_T
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
      #
      rm -f mtDNA.fa.amb mtDNA.fa.ann mtDNA.fa.bwt mtDNA.fa.pac mtDNA.fa.sa;
      #
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
      ./TRACES_cy_align_reads.sh cy.fa $ORGAN_T $THREADS $HIGH_SENSITIVITY 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
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
      ./TRACES_overall.sh cy cy $ORGAN_T
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
      echo -e "\e[34m[TRACESPipe]\e[32m Running do-novo DNA assembly with metaSPAdes ...\e[0m";
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
      mkdir -p ../output_data/TRACES_viral_consensus/
      mkdir -p ../output_data/TRACES_hybrid/
      mkdir -p ../output_data/TRACES_hybrid_consensus/
      mkdir -p ../output_data/TRACES_hybrid_alignments/
      mkdir -p ../output_data/TRACES_hybrid_bed/
      SCAFFOLDS_PATH="../output_data/TRACES_denovo_$ORGAN_T/scaffolds.fasta";
      ALIGNM_CON_PATH="../output_data/TRACES_viral_consensus";
      HYBRID_CON_PATH="../output_data/TRACES_hybrid_consensus";
      HYBRID_ALI_PATH="../output_data/TRACES_hybrid_alignments";
      HYBRID_BED_PATH="../output_data/TRACES_hybrid_bed";
      HYBRID_PATH="../output_data/TRACES_hybrid";
      CON_PATH="../output_data/TRACES_viral_consensus";
      #
      mkdir -p $HYBRID_CON_PATH;
      mkdir -p $HYBRID_BED_PATH;
      mkdir -p $HYBRID_ALI_PATH;
      mkdir -p ../output_data/TRACES_hybrid_R2_alignments/
      mkdir -p ../output_data/TRACES_hybrid_R2_consensus/
      mkdir -p ../output_data/TRACES_hybrid_R2_bed/
      #
      for VIRUS in "${VIRUSES[@]}"
        do
	#
	if [ -f ../output_data/TRACES_viral_alignments/$ORGAN_T-$VIRUS.fa ];
          then
          echo -e "\e[34m[TRACESPipe]\e[32m Hybrid $VIRUS\e[0m";
          echo -e "\e[34m[TRACESPipe]\e[32m Mode 1\e[0m";
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
          echo -e "\e[34m[TRACESPipe]\e[32m Mode 2\e[0m";
          cp ../output_data/TRACES_hybrid_consensus/$VIRUS-consensus-$ORGAN.fa R2-$ORGAN_T-$VIRUS.fa
          ./TRACES_hybrid_R2.sh R2-$ORGAN_T-$VIRUS.fa $SCAFFOLDS_PATH $VIRUS $ORGAN_T $THREADS 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
          #
          # ROUND 3 and 4
	  #
          echo -e "\e[34m[TRACESPipe]\e[32m Mode 3\e[0m";
	  mkdir -p ../output_data/TRACES_hybrid_R3_alignments/
          mkdir -p ../output_data/TRACES_hybrid_R3_consensus/
          mkdir -p ../output_data/TRACES_hybrid_R3_bed/
          mkdir -p ../output_data/TRACES_hybrid_R4_alignments/
          mkdir -p ../output_data/TRACES_hybrid_R4_consensus/
          mkdir -p ../output_data/TRACES_hybrid_R4_bed/
          #
	  FALCON -F -t 500 -x top-$VIRUS-$ORGAN_T.txt ../output_data/TRACES_hybrid_R2_consensus/$VIRUS-consensus-$ORGAN_T.fa ../output_data/TRACES_denovo_$ORGAN_T/scaffolds.fasta 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
	  BEST_TOP_SIM=`head -n 1 top-$VIRUS-$ORGAN_T.txt | awk '{ print $3;}'`; 
	  #TOP_NLINES_SCAFFOLDS=`wc -l top-$VIRUS-$ORGAN_T.txt`;
	  if [[ "$BEST_TOP_SIM" > "1.0" ]];
	    then
	    cp top-$VIRUS-$ORGAN_T.txt ../output_data/TRACES_hybrid_R3_consensus/
	    FIL_NAME=`cat top-$VIRUS-$ORGAN_T.txt | head -n 1 | awk '{ print $4;}'`; 
	    gto_fasta_extract_read_by_pattern -p "$FIL_NAME" < ../output_data/TRACES_denovo_$ORGAN_T/scaffolds.fasta > $VIRUS-$ORGAN_T-SCAFFOLD.fa 2>> ../logs/Log-stderr-$ORGAN_T.txt;
	    cp $VIRUS-$ORGAN_T-SCAFFOLD.fa ../output_data/TRACES_hybrid_R3_consensus/$VIRUS-consensus-$ORGAN_T.fa
	    #
            echo -e "\e[34m[TRACESPipe]\e[32m Mode 4\e[0m";
            ./TRACES_hybrid_R4.sh $VIRUS-$ORGAN_T-SCAFFOLD.fa ../output_data/TRACES_hybrid_consensus/$VIRUS-consensus-$ORGAN.fa $VIRUS $ORGAN_T $THREADS 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
            fi
          #
          PC_R0=`grep -v ">" ../output_data/TRACES_viral_consensus/$VIRUS-consensus-$ORGAN_T.fa | tr -d -c "ACGTacgt" | tr "acgt" "ACGT" | gto_info | grep "Number of symbols" | awk '{ print $5; }'`;
          if [[ "$PC_R0" == "" ]]; then PC_R0=0; fi

	  PC_R1=`grep -v ">" ../output_data/TRACES_hybrid_consensus/$VIRUS-consensus-$ORGAN_T.fa | tr -d -c "ACGTacgt" | tr "acgt" "ACGT" | gto_info | grep "Number of symbols" | awk '{ print $5; }'`;
	  if [[ "$PC_R1" == "" ]]; then PC_R1=0; fi

	  PC_R2=`grep -v ">" ../output_data/TRACES_hybrid_R2_consensus/$VIRUS-consensus-$ORGAN_T.fa | tr -d -c "ACGTacgt" | tr "acgt" "ACGT" | gto_info | grep "Number of symbols" | awk '{ print $5; }'`;
	  if [[ "$PC_R2" == "" ]]; then PC_R2=0; fi

          if [ ! -f ../output_data/TRACES_hybrid_R3_consensus/$VIRUS-consensus-$ORGAN_T.fa ]; then
            echo "R3 not found, using R2 results."
            cp ../output_data/TRACES_hybrid_R2_consensus/$VIRUS-consensus-$ORGAN_T.fa ../output_data/TRACES_hybrid_R3_consensus/
            
          fi
	  PC_R3=`grep -v ">" ../output_data/TRACES_hybrid_R3_consensus/$VIRUS-consensus-$ORGAN_T.fa | tr -d -c "ACGTacgt" | tr "acgt" "ACGT" | gto_info | grep "Number of symbols" | awk '{ print $5; }'`;
	  if [[ "$PC_R3" == "" ]]; then PC_R3=0; fi
          
          if [ ! -f ../output_data/TRACES_hybrid_R4_consensus/$VIRUS-consensus-$ORGAN_T.fa ]; then
            echo "R4 not found, using R2 results."
            cp ../output_data/TRACES_hybrid_R2_consensus/$VIRUS-consensus-$ORGAN_T.fa ../output_data/TRACES_hybrid_R4_consensus/
            
          fi
	  PC_R4=`grep -v ">" ../output_data/TRACES_hybrid_R4_consensus/$VIRUS-consensus-$ORGAN_T.fa | tr -d -c "ACGTacgt" | tr "acgt" "ACGT" | gto_info | grep "Number of symbols" | awk '{ print $5; }'`;
	  if [[ "$PC_R4" == "" ]]; then PC_R4=0; fi
	  #
	  echo -e "\e[34m[TRACESPipe]\e[32m Bases: {$PC_R0, $PC_R1, $PC_R2, $PC_R3, $PC_R4}\e[0m";
          N_SIZES=("$PC_R0" "$PC_R1" "$PC_R2" "$PC_R3" "$PC_R4")
	  MAX=0;
	  MAX_POSITION=0;
	  IDX_POSITION=0;
	  for i in "${N_SIZES[@]}"; 
	    do 
	    if [[ "$i" -gt "$MAX" ]];
	      then
	      MAX=$i;
	      MAX_POSITION=$IDX_POSITION;
	      fi
	    ((++IDX_POSITION));
            done
	  OUT_R5_PATH="../output_data/TRACES_hybrid_R5_consensus/"
	  mkdir -p $OUT_R5_PATH;
	  echo -e "\e[34m[TRACESPipe]\e[32m Highest index value: $MAX_POSITION\e[0m";
          case "$MAX_POSITION" in
            0)
	    echo -e "\e[34m[TRACESPipe]\e[32m Choosing mode 0 ...\e[0m";
            cp ../output_data/TRACES_viral_consensus/$VIRUS-consensus-$ORGAN_T.fa $OUT_R5_PATH
            ;;
            1)
	    echo -e "\e[34m[TRACESPipe]\e[32m Choosing mode 1 ...\e[0m";
            cp ../output_data/TRACES_hybrid_consensus/$VIRUS-consensus-$ORGAN_T.fa $OUT_R5_PATH
            ;;
            2)
	    echo -e "\e[34m[TRACESPipe]\e[32m Choosing mode 2 ...\e[0m";
	    cp ../output_data/TRACES_hybrid_R2_consensus/$VIRUS-consensus-$ORGAN_T.fa $OUT_R5_PATH
            ;;
            3)
	    echo -e "\e[34m[TRACESPipe]\e[32m Choosing mode 3 ...\e[0m";
	    cp ../output_data/TRACES_hybrid_R3_consensus/$VIRUS-consensus-$ORGAN_T.fa $OUT_R5_PATH
            ;;
            4)
	    echo -e "\e[34m[TRACESPipe]\e[32m Choosing mode 4 ...\e[0m";
	    cp ../output_data/TRACES_hybrid_R4_consensus/$VIRUS-consensus-$ORGAN_T.fa $OUT_R5_PATH
            ;;
            *)
	    echo "ERROR: unknown switch round!"
	    exit;
            esac
            #
          fi
        done
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      fi
    # 
    # DISCLAIMER:
    # THE PREVIOUS CODE IS A HACK TO MERGE THE BEST OF ALIGNMENT-FREE AND 
    # DE-NOVO APPROACHES. IT HAS ADVANTAGES AND DISADVANTAGES. STILL WAITING FOR
    # AN EFFICIENT PROGRAM TO BE DEVELOPED FOR THIS... PERHAPS WE WILL NEED TO 
    # DEPELOP IT.
    #
    # ==========================================================================
    #
    if [[ "$RUN_DIFF" -eq "1" ]];
      then
      #
      # VIRAL DNA
      #
      printf "\n" 1>> ../output_data/TRACES_diff/Viral_Diff.txt;
      echo -e "\e[34m[TRACESPipe]\e[32m Running dnadiff between references and reconstructed ...\e[0m";
      IDX_V=1;
      for VIRUS in "${VIRUSES[@]}"
        do
        if [ -f ../output_data/TRACES_viral_alignments/$ORGAN_T-$VIRUS.fa ];
          then
          if [ -f ../output_data/TRACES_hybrid_R5_consensus/$VIRUS-consensus-$ORGAN_T.fa ];
            then
            # Sanitize dnadiff input with gto
            gto_fasta_to_seq < ../output_data/TRACES_viral_alignments/$ORGAN_T-$VIRUS.fa | gto_fasta_from_seq -l 60 -n AtoCompare > $ORGAN_T-$VIRUS-G_A.fa;
            gto_fasta_to_seq < ../output_data/TRACES_hybrid_R5_consensus/$VIRUS-consensus-$ORGAN_T.fa | gto_fasta_from_seq -l 60 -n BtoBeCompared > $ORGAN_T-$VIRUS-G_B.fa
            #
            dnadiff $ORGAN_T-$VIRUS-G_A.fa $ORGAN_T-$VIRUS-G_B.fa 2>> ../logs/Log-stderr-$ORGAN_T.txt;
	    IDEN=`cat out.report | grep "AvgIdentity "  | head -n 1 | awk '{ print $2;}'`;
            ALBA=`cat out.report | grep "AlignedBases " | head -n 1 | awk '{ print $2;}'`;
            SNPS=`cat out.report | grep TotalSNPs | awk '{ print $2;}'`;
	    XBREADTH=`awk '{ print $3;}' ../output_data/TRACES_viral_statistics/$VIRUS-total-horizontal-coverage-$ORGAN_T.txt`;
            XDEPTH=`awk '{ print $3;}' ../output_data/TRACES_viral_statistics/$VIRUS-total-depth-coverage-$ORGAN_T.txt`;
            if [ -n "$VIRAL_DATABASE_METADATA" ]; then
                V_INFO="$(./TRACES_get_best_by_meta.sh "$ORGAN" "$VIRUS" "$VIRAL_DATABASE_METADATA")";
            else
                V_INFO=`./TRACES_get_best_$VIRUS.sh $ORGAN`;
            fi
            V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
            V_VAL=`echo "$V_INFO" | awk '{ print $1; }'`;
            #
            if [[ "$V_GID" != "-" ]] && [[ "$V_VAL" > "0" ]];
              then
              if [[ "$RUN_BEST_OF_BESTS" -eq "1" ]]
                then
                V_GID=`sed "${IDX_V}q;d" ../output_data/TRACES_results/REPORT_META_VIRAL_BESTS_$ORGAN_T.txt | awk '{ print $2; }'`;
                fi
              fi
	    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$ORGAN_T" "$VIRUS" "$V_GID" "$V_VAL" "$ALBA" "$IDEN" "$SNPS" "$XBREADTH" "$XDEPTH" 1>> ../output_data/TRACES_diff/Viral_Diff.txt
            rm -f $ORGAN_T-$VIRUS-G_A.fa $ORGAN_T-$VIRUS-G_B.fa ;
            fi
          fi
	  ((++IDX_V));
        done
      #
      # MT DNA
      #
      if [ -f ../output_data/TRACES_mtdna_alignments/mtDNA.fa ];
        then
        if [ -f ../output_data/TRACES_mtdna_consensus/mt-consensus-$ORGAN_T.fa ];
          then
	  # Sanitize dnadiff input
	  gto_fasta_to_seq < ../output_data/TRACES_mtdna_alignments/mtDNA.fa | gto_fasta_from_seq -l 60 -n AtoCompare > MT-G_A.fa;
          gto_fasta_to_seq < ../output_data/TRACES_mtdna_consensus/mt-consensus-$ORGAN_T.fa | gto_fasta_from_seq -l 60 -n BtoBeCompared > $ORGAN_T-MT-G_B.fa;
          #
	  cp ../output_data/TRACES_mtdna_alignments/mtDNA.fa MT-G_A.fa;
          cp ../output_data/TRACES_mtdna_consensus/mt-consensus-$ORGAN_T.fa $ORGAN_T-MT-G_B.fa;
          dnadiff MT-G_A.fa $ORGAN_T-MT-G_B.fa 2>> ../logs/Log-stderr-$ORGAN_T.txt;
          IDEN=`cat out.report | grep "AvgIdentity " | head -n 1 | awk '{ print $2;}'`;
          ALBA=`cat out.report | grep "AlignedBases " | head -n 1 | awk '{ print $2;}'`;
          SNPS=`cat out.report | grep TotalSNPs | awk '{ print $2;}'`;
	  XBREADTH=`awk '{ print $3;}' ../output_data/TRACES_mtdna_statistics/mt-total-horizontal-coverage-$ORGAN_T.txt`;
          XDEPTH=`awk '{ print $3;}' ../output_data/TRACES_mtdna_statistics/mt-total-depth-coverage-$ORGAN_T.txt`;
	  printf "%s\tMT\t-\t%s\t%s\t%s\t%s\t%s\n" "$ORGAN_T" "$ALBA" "$IDEN" "$SNPS" "$XBREADTH" "$XDEPTH" 1>> ../output_data/TRACES_diff/mtDNA_Diff.txt
          rm -f MT-G_A.fa $ORGAN_T-MT-G_B.fa ;
          fi
        fi
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      fi
    #
    # ==========================================================================
    # RUN SPECIFIC ALIGN/CONSENSUS USING DE-NOVO
    #
    if [[ "$RUN_DENOVO_SPECIFIC" -eq "1" ]];
      then
      echo -e "\e[34m[TRACESPipe]\e[32m Aligning de-novo scaffolds to specific viral ref(s) with pattern \"$SPECIFIC_DENOVO_ID\" using bowtie2 ...\e[0m";
      #
      CHECK_VDB;
      #
      echo -e "\e[34m[TRACESPipe]\e[32m Extracting sequence with pattern \"$SPECIFIC_DENOVO_ID\" from "$VIRAL_DATABASE_FILE" ...\e[0m";
      gto_fasta_extract_read_by_pattern -p "$SPECIFIC_DENOVO_ID" < "$VIRAL_DATABASE_FILE" | awk "/^>/ {n++} n>1 {exit} 1" > SPECIFIC-$SPECIFIC_DENOVO_ID.fa
      echo -e "\e[34m[TRACESPipe]\e[32m Aligning and creating consensus ... \e[0m";
      SCAFFOLDS_PATH="../output_data/TRACES_denovo_$ORGAN_T/scaffolds.fasta";
      ./TRACES_denovo_specific.sh SPECIFIC-$SPECIFIC_DENOVO_ID.fa $SCAFFOLDS_PATH $SPECIFIC_DENOVO_ID $ORGAN_T $THREADS 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      #
      mkdir -p ../output_data/TRACES_specific_denovo_alignments;
      cp SPECIFIC-$SPECIFIC_DENOVO_ID.fa ../output_data/TRACES_specific_denovo_alignments/
      cp SPECIFIC-$SPECIFIC_DENOVO_ID.fa.fai ../output_data/TRACES_specific_denovo_alignments/
      mv specific_aligned_sorted_$SPECIFIC_DENOVO_ID-$ORGAN_T.bam ../output_data/TRACES_specific_denovo_alignments/
      mv specific_aligned_sorted_$SPECIFIC_DENOVO_ID-$ORGAN_T.bam.bai ../output_data/TRACES_specific_denovo_alignments/
      mkdir -p ../output_data/TRACES_specific_denovo_consensus;
      mv $SPECIFIC_DENOVO_ID-consensus-$ORGAN_T.fa ../output_data/TRACES_specific_denovo_consensus/
      mkdir -p ../output_data/TRACES_specific_denovo_bed;
      mv $SPECIFIC_DENOVO_ID-calls-$ORGAN_T.bed ../output_data/TRACES_specific_denovo_bed/
      mv $SPECIFIC_DENOVO_ID-coverage-$ORGAN_T.bed ../output_data/TRACES_specific_denovo_bed/
      mv $SPECIFIC_DENOVO_ID-zero-coverage-$ORGAN_T.bed ../output_data/TRACES_specific_denovo_bed/
      echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
      fi
    #
    # ==========================================================================
    #
    done
    #
  ###
  #
  # ==========================================================================
  #
  if [[ "$RUN_MULTIORGAN_CONSENSUS" -eq "1" ]];
    then
    echo -e "\e[34m[TRACESPipe]\e[32m Running multi-organ (MO) consensus ...\e[0m";
    mkdir -p ../output_data/TRACES_multiorgan_alignments/
    mkdir -p ../output_data/TRACES_multiorgan_consensus/
    R5PATH="../output_data/TRACES_hybrid_R5_consensus";
    for VIRUS in "${VIRUSES[@]}"
      do
      echo -e "\e[34m[TRACESPipe]\e[32m Running MO consensus for $VIRUS ...\e[0m";
      #
      MAX_NB=0;
      MAX_RF="organ";
      #
      for read in "${READS[@]}" #
        do
        #
        ORGAN=`echo $read | tr ':' '\t' | awk '{ print $1 }'`;
	#
	if [ -f $R5PATH/$VIRUS-consensus-$ORGAN.fa ];
          then
          NS_V=`gto_fasta_to_seq < $R5PATH/$VIRUS-consensus-$ORGAN.fa \
	      | tr -d "N" \
              | gto_info \
              | grep "Number of symbols" \
	      | awk '{ print $5; }'`;
	  #
          if [[ "$NS_V" == "" ]]; 
	    then 
	    NS_V=0; 
            fi
	  #
          if [[ "$NS_V" -gt "$MAX_NB" ]];
	    then
            MAX_NB=$NS_V;
	    MAX_RF=$ORGAN;
	    fi
	  fi
        done
      #
      if [ -f $R5PATH/$VIRUS-consensus-$MAX_RF.fa ];
        then
        cp $R5PATH/$VIRUS-consensus-$MAX_RF.fa $VIRUS.fa
        #
	NSEQS=`ls  $R5PATH/$VIRUS-consensus-*.fa | wc -l`;
        if [[ "$NSEQS" -gt "1" ]]
          then         
          #
          rm -f $VIRUS-multiorgans.fa;
	  #
	  for read in "${READS[@]}" #
            do
            #
            ORGAN=`echo $read | tr ':' '\t' | awk '{ print $1 }'`;
            #
            if [[ "$ORGAN" != "$MAX_RF" ]];
              then
	      if [ -f $R5PATH/$VIRUS-consensus-$ORGAN.fa ];
	        then
                cat $R5PATH/$VIRUS-consensus-$ORGAN.fa >> $VIRUS-multiorgans.fa;
	        fi
              fi
            #
            done
	  #
          if [ -f $VIRUS-multiorgans.fa ];
	    then
            ./TRACES_multiorgan_consensus.sh $VIRUS $VIRUS.fa $VIRUS-multiorgans.fa $THREADS 1>> ../logs/Log-stdout-$ORGAN_T.txt 2>> ../logs/Log-stderr-$ORGAN_T.txt;
            cp $VIRUS.fa ../output_data/TRACES_multiorgan_alignments/$VIRUS.fa
            cp $VIRUS-data_aligned_sorted.bam.bai ../output_data/TRACES_multiorgan_alignments/
            cp $VIRUS-data_aligned_sorted.bam ../output_data/TRACES_multiorgan_alignments/
            cp $VIRUS-multiorgan-consensus.fa ../output_data/TRACES_multiorgan_consensus/
            fi
	  #
          fi
        fi
      done
    echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
    fi
  #
  # ==========================================================================
  #
  if [[ "$RUN_BLAST_RECONSTRUCTED" -eq "1" ]];
    then
    echo -e "\e[34m[TRACESPipe]\e[32m Running local blast for reconstructed data ...\e[0m";
    #
    BLASTS_OUTPUT="../output_data/TRACES_blasts/";
    mkdir -p $BLASTS_OUTPUT;
    R5PATH="../output_data/TRACES_hybrid_R5_consensus";
    #
    if [[ "$RUN_DIFF" -eq "1" ]];
      then
      printf "\n" 1>> ../output_data/TRACES_diff/Viral_Diff_after_blastn.txt;
      fi
    #
    for VIRUS in "${VIRUSES[@]}"
      do
      #
      for read in "${READS[@]}" #
        do
        ORGAN=`echo $read | tr ':' '\t' | awk '{ print $1 }'`;
        if [ -f $R5PATH/$VIRUS-consensus-$ORGAN.fa ];
          then
          ./TRACES_blastn_n_db.sh $R5PATH/$VIRUS-consensus-$ORGAN.fa > $BLASTS_OUTPUT/SINGLE-LIST-$VIRUS-$ORGAN.txt 2>> ../logs/Log-stderr-system.txt;
          head -n 1 $BLASTS_OUTPUT/SINGLE-LIST-$VIRUS-$ORGAN.txt > $BLASTS_OUTPUT/SINGLE-BEST-$VIRUS-$ORGAN.txt;
          #
          if [[ "$RUN_DIFF" -eq "1" ]];
            then
            IDBLAST=`head -n 1 $BLASTS_OUTPUT/SINGLE-LIST-$VIRUS-$ORGAN.txt | awk '{ print $2}' | tr '|' '\t' | awk '{ print $4;}'`;
            #
            # Try to find genome in "$VIRAL_DATABASE_FILE"
            CHECK_VDB;
            gto_fasta_extract_read_by_pattern -p "$IDBLAST" < "$VIRAL_DATABASE_FILE" \
            | awk "/^>/ {n++} n>1 {exit} 1" > TMP-BEST-BLAST.fa
            NVALSEQ=`wc -l TMP-BEST-BLAST.fa | awk '{ print $1;}'`;
            if [[ "$NVALSEQ" -eq "0" ]];
              then
              efetch -db nucleotide -format fasta -id "$IDBLAST" > TMP-BEST-BLAST.fa
              fi
            # Sanitize inputs to dnadiff     
            gto_fasta_to_seq < TMP-BEST-BLAST.fa \
            | gto_fasta_from_seq -l 60 -n AtoCompare > $ORGAN-$VIRUS-G_A.fa 2>> ../logs/Log-$ORGAN.txt;
            #
            gto_fasta_to_seq < $R5PATH/$VIRUS-consensus-$ORGAN.fa \
            | gto_fasta_from_seq -l 60 -n BtoBeCompared > $ORGAN-$VIRUS-G_B.fa
            #
            dnadiff $ORGAN-$VIRUS-G_A.fa $ORGAN-$VIRUS-G_B.fa 2>> ../logs/Log-stderr-$ORGAN.txt;
            #
            IDEN=`cat out.report | grep "AvgIdentity "  | head -n 1 | awk '{ print $2;}'`;
            ALBA=`cat out.report | grep "AlignedBases " | head -n 1 | awk '{ print $2;}'`;
            SNPS=`cat out.report | grep TotalSNPs | awk '{ print $2;}'`;
            XBREADTH=`awk '{ print $3;}' ../output_data/TRACES_viral_statistics/$VIRUS-total-horizontal-coverage-$ORGAN.txt`;
            XDEPTH=`awk '{ print $3;}' ../output_data/TRACES_viral_statistics/$VIRUS-total-depth-coverage-$ORGAN.txt`;
            if [ -n "$VIRAL_DATABASE_METADATA" ]; then
                V_INFO="$(./TRACES_get_best_by_meta.sh "$ORGAN" "$VIRUS" "$VIRAL_DATABASE_METADATA")";
            else
                V_INFO=`./TRACES_get_best_$VIRUS.sh $ORGAN`;
            fi
            V_GID=`echo "$V_INFO" | awk '{ print $2; }'`;
            V_VAL=`echo "$V_INFO" | awk '{ print $1; }'`;
            #
            if [[ "$V_GID" != "-" ]] && [[ "$V_VAL" > "0" ]];
              then
              if [[ "$RUN_BEST_OF_BESTS" -eq "1" ]]
                then
                V_GID=`sed "${IDX_V}q;d" ../output_data/TRACES_results/REPORT_META_VIRAL_BESTS_$ORGAN.txt | awk '{ print $2; }'`;
                else
                V_GID=$IDBLAST;
                fi
              printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$ORGAN" "$VIRUS" "$V_GID" "$V_VAL" "$ALBA" "$IDEN" "$SNPS" "$XBREADTH" "$XDEPTH" 1>> ../output_data/TRACES_diff/Viral_Diff_after_blastn.txt
              fi
            rm -f $ORGAN-$VIRUS-G_A.fa $ORGAN-$VIRUS-G_B.fa ;
            #
            fi
          #
          fi
        done
      #    
      if [ -f ../output_data/TRACES_multiorgan_consensus/$VIRUS-multiorgans.fa ];
        then
        ./TRACES_blastn_n_db.sh ../output_data/TRACES_multiorgan_consensus/$VIRUS-multiorgans.fa > $BLASTS_OUTPUT/MULTI-LIST-$VIRUS.txt 2>> ../logs/Log-stderr-system.txt;
        head -n 1 $BLASTS_OUTPUT/MULTI-LIST-$VIRUS.txt > $BLASTS_OUTPUT/MULTI-BEST-$VIRUS.txt;
        fi
      #
      done
    echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
    fi
  #
  # ==========================================================================
  #
  # BUILD COMPLETE VIRAL META TABLE FOR MULTIPLE ORGANS:
  #
  if [[ "$RUN_META_ON" -eq "1" ]];
    then
    ./TRACES_get_report_meta.sh
    fi
  #
  # ==============================================================================
  # CLEAN DATA:
  rm -f FW_READS.fq.gz RV_READS.fq.gz
  #
  # ============================================================================
  fi
#
# ==============================================================================
################################################################################
