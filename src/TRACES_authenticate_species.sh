#!/bin/bash
#
#
THREADS=8;
OUTPUT="../src";
#
DOWNLOAD=1;
RUN=1;
#
# ==============================================================================
#
if [[ "$DOWNLOAD" -eq "1" ]];
  then
  echo -e "\e[34m[TRACESPipe]\e[36m Downloading human mitogenomes database ...\e[0m";
  #
  rm -f alg_tot.fasta alg_pa_tot.fasta alg_tot.fasta alg_pa_tot.fasta
  #
  wget https://www.hmtdb.uniba.it/download_file?dataset=alg_tot.zip -O alg_tot.zip
  wget https://www.hmtdb.uniba.it/download_file?dataset=alg_pa_tot.zip -O alg_pa_tot.zip
  #
  unzip alg_tot.zip
  unzip alg_pa_tot.zip
  #
  echo -e "\e[34m[TRACESPipe]\e[36m Done!\e[0m";
  fi
#
# ==============================================================================
#
if [[ "$RUN" -eq "1" ]];
  then
  #
  echo -e "\e[34m[TRACESPipe]\e[36m Running human mitogenome authentication ...\e[0m";
  #
  if [ ! -f ../meta_data/meta_info.txt ];
    then
    echo -e "\e[31mERROR: meta_info.txt file not found!\e[0m"
    echo "Please create a meta information file before the run."
    echo "For addition information, see the instructions at the web page."
    exit 1;
    fi
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
    FALCON -c 30 -l 47 -t 200 -x $OUTPUT/MT-HEALT-$ORGAN.txt -n $THREADS reads_human_auth.fq alg_tot.fasta
    FALCON -c 30 -l 47 -t 200 -x $OUTPUT/MT-PATHO-$ORGAN.txt -n $THREADS reads_human_auth.fq alg_pa_tot.fasta
    done
  echo -e "\e[34m[TRACESPipe]\e[36m Done!\e[0m";
  fi
#
# ==============================================================================
#
