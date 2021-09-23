#!/bin/bash
#
RUN="0";
INSTALL="0";
THREADS="8";
READS1="";
READS2="";
GID="NC_012920.1"
OUTPUT="damage_analysis";
#
################################################################################
#
SHOW_MENU () {
  echo " -------------------------------------------------------- ";
  echo "                                                          ";
  echo " TRACESPipe_damage_analysis.sh : damage pattern analysis  ";
  echo "                                                          ";
  echo " This subpipeline of TRACESPipe automatically generates   ";
  echo " plots for analysis of damage patterns using mapDamage2.  ";
  echo " Any reference sequence can be used through the GID.      ";
  echo "                                                          ";
  echo " Program options ---------------------------------------- ";
  echo "                                                          ";
  echo " -h, --help                  Show this,                   ";
  echo " -i, --install               Installation (w/ conda),     ";
  echo "                                                          ";
  echo " -t  <INT>, --threads <INT>  Number of threads,           ";
  echo " -g  <STR>, --gid <STR>      GID of genome/sequence,      ";
  echo " -o  <STR>, --output <STR>   output folder name,          ";
  echo "                                                          ";
  echo " -r1 <STR>, --reads1 <STR>   FASTQ reads (forward),       ";
  echo " -r2 <STR>, --reads2 <STR>   FASTQ reads (reverse).       ";
  echo "                                                          ";
  echo " Example -----------------------------------------------  ";
  echo "                                                          ";
  echo " TRACESPipe_damage_analysis.sh --gid NC_012920.1 \\       ";
  echo " --reads1 reads_fw.fq.gz --reads2 reads_rv.fq.gz \\       ";
  echo " --threads 8 --output damage_analysis                     ";
  echo "                                                          ";
  echo " -------------------------------------------------------  ";
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
  printf "Checking $1 ... ";
  if ! [ -x "$(command -v $1)" ];
    then
    echo -e "\e[41mERROR\e[49m: $1 is not installed." >&2;
    echo -e "\e[42mTIP\e[49m: Try: TRACESPipe_damage_analysis.sh --install" >&2;
    exit 1;
    else
    echo -e "\e[42mSUCCESS!\e[49m";
    fi
  }
#
################################################################################
#
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
    -c|--install|--compile)
      INSTALL=1;
      shift
    ;;
    -t|--threads)
      THREADS="$2";
      shift 2;
    ;;
    -r1|-r|--input1|--reads|--reads1)
      READS1="$2";
      RUN=1;
      shift 2;
    ;;
    -r2|--input2|--reads2)
      READS2="$2";
      RUN=1;
      shift 2;
    ;;
    -g|--gid)
      GID="$2";
      shift 2;
    ;;
    -o|--output)
      OUTPUT="$2";
      shift 2;
    ;;
    -*) # unknown option with small
    echo "Invalid arg ($1)!";
    echo "For help, try: TRACESPipe_damage_analysis.sh -h"
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
  #
  conda install -c bioconda AdapterRemoval -y
  conda install -c bioconda entrez-direct -y
  conda install -c bioconda bwa -y
  conda install -c bioconda mapdamage2=2.2.1 -y
  conda install -c bioconda samtools -y
  #
  PROGRAM_EXISTS "AdapterRemoval";
  PROGRAM_EXISTS "efetch";
  PROGRAM_EXISTS "bwa";
  PROGRAM_EXISTS "samtools";
  PROGRAM_EXISTS "mapDamage";
  #
  echo "Generating adapters for AdapterRemoval ...";
  echo "TACACTCTTTCCCTACACGACGCTCTTCCGATCT      AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA" >  adapters_ar.fa;
  echo "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT      AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" >> adapters_ar.fa;
  echo "Done!";
  exit;
  #
  fi
#
################################################################################
#
if [[ "$RUN" -eq "1" ]];
  then
  #
  CHECK_INPUT $READS1
  CHECK_INPUT adapters_ar.fa
  #
  rm -fr $GID.fa* $GID-READS.sam $GID-READS.bam FIL-$GID-READS.bam \
  SORT-FIL-$GID-READS.bam SORT-FIL-$GID-READS.bam.bai \
  reads-tracespipe-run-tmp.fq $OUTPUT/
  #
  efetch -db nucleotide -format fasta -id "$GID" > $GID.fa
  #
  if [[ $READS2 = *[!\ ]* ]];
    then
    AdapterRemoval --threads $THREADS --file1 $READS1 --file2 $READS2 \
    --outputcollapsed reads-tracespipe-run-tmp.fq --trimns --trimqualities \
    --minlength 20 --collapse --adapter-list adapters_ar.fa;
    else
    AdapterRemoval --threads $THREADS --file1 $READS1 \
    --outputcollapsed reads-tracespipe-run-tmp.fq --trimns --trimqualities \
    --minlength 20 --collapse --adapter-list adapters_ar.fa;
    fi
  #
  bwa index $GID.fa
  #
  bwa aln -l 1000 -n 0.01 $GID.fa reads-tracespipe-run-tmp.fq > $GID-READS.sai
  bwa samse $GID.fa $GID-READS.sai reads-tracespipe-run-tmp.fq > $GID-READS.sam
  #
  samtools view -bSh $GID-READS.sam > $GID-READS.bam;
  samtools view -bh -F4 $GID-READS.bam > FIL-$GID-READS.bam;
  samtools sort -o SORT-FIL-$GID-READS.bam FIL-$GID-READS.bam;
  samtools rmdup -s SORT-FIL-$GID-READS.bam RD-SORT-FIL-$GID-READS.bam;
  samtools view -h RD-SORT-FIL-$GID-READS.bam \
        | grep -v 'XT:A:R'\
        | grep -v 'XA:Z' \
        | grep -v 'XT:A:M' \
        | awk '{if($0~/X1:i:0/||$0~/^@/  )print $0}' \
        | samtools view -bS - > UNK-RD-SORT-FIL-$GID-READS.bam;
  #
  samtools index -b UNK-RD-SORT-FIL-$GID-READS.bam UNK-RD-SORT-FIL-$GID-READS.bam.bai
  #
  mapDamage --rescale -d $OUTPUT -i UNK-RD-SORT-FIL-$GID-READS.bam -r $GID.fa
  #
  rm -f $GID.fa* $GID-READS.sam $GID-READS.bam FIL-$GID-READS.bam \
  SORT-FIL-$GID-READS.bam SORT-FIL-$GID-READS.bam.bai \ 
  RD-SORT-FIL-$GID-READS.bam reads-tracespipe-run-tmp.fq $GID.fa.amb \
  $GID.fa.ann $GID.fa.bwt $GID.fa.fai $GID.fa.pac $GID.fa.sa;
  #
  fi
#
################################################################################
#

