#!/bin/bash
#
THREADS=8;
OUTPUT="../src";
#
DOWNLOAD=1;
RUN=1;
#
# ==============================================================================
#
CHECK_FILE () {
  if [ ! -f $1 ];
    then
    echo -e "\e[31mERROR: file $1 not found!\e[0m"
    echo "For addition information, see the instructions at the web page."
    exit 1;
    fi
  }
#
# ==============================================================================
#
if [[ "$DOWNLOAD" -eq "1" ]];
  then
  echo -e "\e[34m[TRACESPipe]\e[36m Downloading species mitogenomes database ...\e[0m";
  #
  rm -f mitochondrion.1.1.genomic.fna.gz mitochondrion.2.1.genomic.fna.gz MTs.fa
  #
  wget https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.1.genomic.fna.gz
  wget https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.2.1.genomic.fna.gz
  zcat mitochondrion.1.1.genomic.fna.gz mitochondrion.2.1.genomic.fna.gz > SPECIES-MTs.fa
  #
  echo -e "\e[34m[TRACESPipe]\e[36m Done!\e[0m";
  fi
#
# ==============================================================================
#
if [[ "$RUN" -eq "1" ]];
  then
  #
  echo -e "\e[34m[TRACESPipe]\e[36m Running species mitogenome authentication ...\e[0m";
  #
  CHECK_FILE "../meta_data/meta_info.txt";
  CHECK_FILE "SPECIES-MTs.fa";
  #
  mapfile -t READS < ../meta_data/meta_info.txt
  #
  for read in "${READS[@]}" #
    do
    #
    ORGAN=`echo $read | tr ':' '\t' | awk '{ print $1 }'`;
    SPL_Forward=`echo $read | tr ':' '\t' | awk '{ print $2 }'`;
    SPL_Reverse=`echo $read | tr ':' '\t' | awk '{ print $3 }'`;
    #
    echo -e "\e[34m[TRACESPipe]\e[36m Running $ORGAN sample ...\e[0m";
    zcat ../input_data/$SPL_Forward ../input_data/$SPL_Reverse > reads_human_auth.fq
    FALCON -c 30 -l 47 -t 200 -x $OUTPUT/MT-SPECIES-$ORGAN.txt -n $THREADS reads_human_auth.fq SPECIES-MTs.fa
    done
  echo -e "\e[34m[TRACESPipe]\e[36m Done!\e[0m";
  fi
#
# ==============================================================================
#
