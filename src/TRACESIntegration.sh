#!/bin/bash
#
HELP="0";
INSTALL="0";
OUTPUT="G";
RUN="0";
THREADS="8";
READS1="";
READS2="";
HOST="";
VIRUS="";
#
################################################################################
#
SHOW_MENU () {
  echo " ------------------------------------------------------- ";
  echo "                                                         ";
  echo " TRACESIntegration.sh >> viral integration analysis     ";
  echo "                                                         ";
  echo " Program options --------------------------------------- ";
  echo "                                                         ";
  echo " -h, --help                    Show this,                ";
  echo " -i, --install                 Installation (w/ conda),  ";
  echo "                                                         ";
  echo " -n <STR>,  --name <STR>       TAG name for uniqueness,  ";
  echo " -t <INT>,  --threads <INT>    Number of threads,        ";
  echo "                                                         ";
  echo " Input options ----------------------------------------- ";
  echo "                                                         ";
  echo " -v <ID>,    --virus <ID>      Virus ID (NCBI>efetch),   ";
  echo "                                                         ";
  echo " -hs <FILE>, --host <FILE>     Host FASTA filename,      ";
  echo "                                                         ";
  echo " -f1 <FILE>, --forward <FILE>  FASTQ forward filename,   ";
  echo " -f2 <FILE>, --reverse <FILE>  FASTQ reverse filename.   ";
  echo "                                                         ";
  echo " Example ----------------------------------------------- ";
  echo "                                                         ";
  echo " ./TRACESIntegration.sh --threads 8 --name out1  \\     ";
  echo "   --virus U31781.1 --host hs.fna --forward a1.fq \\     ";
  echo "   --reverse a2.fq                                       ";
  echo "                                                         ";
  echo " ------------------------------------------------------- ";
  }
#
################################################################################
#
CHECK_INPUT () {
  FILE=$1
  if [ -f "$FILE" ];
    then
    echo "Input filename: $FILE"
    else
    echo -e "\e[31mERROR: input file not found ($FILE)!\e[0m";
    SHOW_MENU;
    exit;
    fi
  }
#
################################################################################
#
PROGRAM_EXISTS () {
  if ! [ -x "$(command -v $1)" ];
    then
    echo -e "\e[41mERROR\e[49m: $1 is not installed." >&2;
    echo -e "\e[42mTIP\e[49m: Try: TRACESIntegration.sh --install" >&2;
    exit 1;
    fi
  }
#
################################################################################
#
CHECK_PROGRAMS () {
  PROGRAM_EXISTS "trimmomatic";
  PROGRAM_EXISTS "bwa";
  PROGRAM_EXISTS "samtools";
  PROGRAM_EXISTS "efetch";
  }
#
################################################################################
#
if [[ "$#" -lt 1 ]];
  then
  HELP=1;
  fi
#
POSITIONAL=();
#
while [[ $# -gt 0 ]]
  do
  i="$1";
  case $i in
    -h|--help|?)
      HELP=1;
      shift
    ;;
    -i|--install|--compile)
      INSTALL=1;
      shift
    ;;
    -n|--name|--tag)
      OUTPUT="$2";
      shift 2;
    ;;
    -t|--threads)
      THREADS="$2";
      shift 2;
    ;;
    -v|--virus)
      VIRUS="$2";
      shift 2;
    ;;
    -hs|--host)
      HOST="$2";
      shift 2;
    ;;
    -f1|--forward)
      READS1="$2";
      RUN=1;
      shift 2;
    ;;
    -f2|--reverse)
      READS2="$2";
      RUN=1;
      shift 2;
    ;;
    -*) # unknown option with small
    echo "Invalid arg ($1)!";
    echo "For help, try: TRACESIntegration.sh -h"
    exit 1;
    ;;
  esac
  done
#
set -- "${POSITIONAL[@]}" # restore positional parameters
#
################################################################################
#
if [[ "$HELP" -eq "1" ]];
  then
  SHOW_MENU;
  exit;
  fi
#
################################################################################
#
if [[ "$INSTALL" -eq "1" ]];
  then
  # ACTIVATE CONDA
  conda create --name python python=2.7 --yes
  conda activate python
  python -m pip install pyfaidx
  python -m pip install pysam
  export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib"
  #
  conda config --add channels conda-forge
  conda config --add channels bioconda
  conda config --add channels defaults
  #
  conda install -c bioconda trimmomatic --yes
  conda install -c bioconda bwa --yes
  conda install -c bioconda samtools=1.9 --yes
  conda install -c bioconda entrez-direct --yes
  #
  rm -fr master.zip SurVirus-master/
  wget https://github.com/kensung-lab/SurVirus/archive/refs/heads/master.zip
  unzip master.zip
  cd SurVirus-master/
  ./build_libs.sh
  cmake -DCMAKE_BUILD_TYPE=Release . && make
  cd ..
  #
  CHECK_PROGRAMS
  #
  printf ">PrefixPE/1\nTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n>PrefixPE/2\nGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT\n>PE1\nTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n>PE1_rc\nAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA\n>PE2\nGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT\n>PE2_rc\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC\n" > adapters.fa
  echo "Done!";
  exit;
  fi
#
################################################################################
#
if [[ "$RUN" -eq "1" ]];
  then
  #
  cp $READS1 SurVirus-master/
  cp $READS2 SurVirus-master/
  cp $HOST SurVirus-master/
  cp adapters.fa SurVirus-master/
  #
  cd SurVirus-master/
  #
  echo "RUNNING..."
  #
  CHECK_INPUT $READS1
  CHECK_INPUT $READS2
  CHECK_INPUT adapters.fa
  CHECK_INPUT $HOST
  #
  CHECK_PROGRAMS
  #
  efetch -db nucleotide -format FASTA -id $VIRUS > $VIRUS.fa
  CHECK_INPUT $VIRUS.fa
  #
  rm ../$OUTPUT/* -fr
  mkdir -p ../$OUTPUT/
  #
  ## RUN TRIMMING ==============================================================
  #
  rm -f o_fw_pr.fq o_fw_unpr.fq o_rv_pr.fq o_rv_unpr.fq;
  #
  trimmomatic PE -threads $THREADS -phred33 $READS1 $READS2 o_fw_pr.fq o_fw_unpr.fq o_rv_pr.fq o_rv_unpr.fq ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30
  #
  # ----------------------------------------------------------------------------
  #
  rm -fr $HOST+$VIRUS.fa;
  #
  bwa index -b 10000000 $HOST
  samtools faidx $HOST
  #
  bwa index $VIRUS.fa
  samtools faidx $VIRUS.fa
  #
  cat $HOST $VIRUS.fa > $HOST+$VIRUS.fa
  bwa index -b 10000000 $HOST+$VIRUS.fa
  samtools faidx $HOST+$VIRUS.fa
  #
  # ----------------------------------------------------------------------------
  #
  python surveyor.py o_fw_pr.fq,o_rv_pr.fq ../$OUTPUT $HOST $VIRUS.fa $HOST+$VIRUS.fa --fq
  #
  cp $VIRUS.fa* ../$OUTPUT/
  cp $HOST.* ../$OUTPUT/
  cp o_fw_pr.fq ../$OUTPUT/
  cp o_rv_pr.fq ../$OUTPUT/
  cp $HOST+$VIRUS.* ../$OUTPUT/
  #
  cd ..
  #
  fi
#
#===============================================================================
#
