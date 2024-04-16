#!/bin/bash
#
# 1 -> TURNS THE RESPECTIVE COMPUTATION ON
# 0 -> IGNORES THE RESPECTIVE ACTION
#
GET_UTILS=0;
INSTALL_PIPELINE=1;
#
# ==============================================================================
# ################################# FUNCTIONS ##################################
# ==============================================================================
#
Program_installed () {
  printf "Checking $1 ... ";
  if ! [ -x "$(command -v $1)" ];
    then
    echo -e "\e[41mERROR\e[49m: $1 is not installed." >&2;
#    exit 1;
    else
    echo -e "\e[42mSUCCESS!\e[49m";
    fi
  }
#
# ==============================================================================
# ================================= GET UTILS ==================================
if [[ "$GET_UTILS" -eq "1" ]];
  then
  sudo apt-get install librsvg2-bin cmake git python3-pip curl
  pip3 install conda
  fi
#
# ==============================================================================
# ============================== INSTALL PIPELINE ==============================
if [[ "$INSTALL_PIPELINE" -eq "1" ]];
  then
  conda config --add channels conda-forge
  conda config --add channels bioconda
  conda config --add channels defaults
  #
  conda install -c bioconda trimmomatic --yes
  conda install -c seyedmorteza cryfa --yes
  conda install -c cobilab magnet --yes
  conda install -c cobilab falcon --yes
  conda install -c cobilab gto --yes
  conda install -c bioconda spades --yes
  conda install -c bioconda igv --yes
  conda install -c bioconda bowtie2 --yes
  conda install -c bioconda samtools=1.9 --yes
  conda install -c bioconda bcftools=1.9 --yes
  conda install -c bioconda bedops --yes
  conda install -c bioconda bedtools --yes
  conda install -c bioconda fastq-pair --yes
  conda install -c bioconda entrez-direct --yes
  conda install -c bioconda/label/cf201901 entrez-direct --yes
  conda install -c bioconda mapdamage2 --yes
 # conda install -c bioconda mapdamage2=2.1.1=pyr36_1
  conda install -c bioconda tabix --yes
  conda install -c bioconda adapterremoval --yes
  conda install -c bioconda bwa --yes
  conda install -c bioconda art --yes
  conda install -c bioconda blast --yes
  conda install -c bioconda mummer4 --yes
#  conda install -c bioconda ivar --yes
  #
  Program_installed "trimmomatic";
  Program_installed "cryfa";
  Program_installed "MAGNET";
  Program_installed "FALCON";
  Program_installed "gto";
  Program_installed "spades.py";
  Program_installed "igv";
  Program_installed "bowtie2";
  Program_installed "samtools";
  Program_installed "bcftools";
  Program_installed "bedops";
  Program_installed "bedtools";
  Program_installed "fastq_pair";
  Program_installed "efetch";
  Program_installed "mapDamage";
  Program_installed "tabix";
  Program_installed "AdapterRemoval";
  Program_installed "bwa";
  Program_installed "art_illumina";
  Program_installed "blastn";
  Program_installed "dnadiff";
#  Program_installed "ivar";
  fi
#
# ==============================================================================

