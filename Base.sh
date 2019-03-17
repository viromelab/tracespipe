#!/bin/bash
#
GET_TOOLS=1;
#
GET_EXTERNAL_B_TOOLS=1;
#
GET_FALCON=1;
GET_GULL=1;
GET_GOOSE=1;
GET_MAGNET=1;
GET_GECO=1;
GET_AC=1;
GET_GECO=1;
GET_CHESTER=1;
#
GET_CRYFA=1;
GET_SMASHPP=1;
GET_FRUIT=1;
#
#==============================================================================
#==============================================================================
#==============================================================================
# GET TOOLS
#
if [[ "$GET_TOOLS" -eq "1" ]]; then
  sudo apt-get update
  sudo apt-get install vim terminator git cmake conda
fi
#
#==============================================================================
#==============================================================================
#==============================================================================
#
if [[ "$GET_EXTERNAL_B_TOOLS" -eq "1" ]]; then
  #
  conda config --add channels defaults
  conda config --add channels bioconda
  conda config --add channels conda-forge
  #
  conda install -c bioconda bwa
  conda install -c bioconda bowtie
  conda install -c bioconda bowtie2 
  conda install -c bioconda mapDamage2 
  conda install -c bioconda centrifuge
  conda install -c bioconda spades 
  conda install -c bioconda bamtools 
  conda install -c bioconda samtools 
  conda install -c bioconda vcftools
  conda install -c bioconda fastqc 
  conda install -c bioconda fqzcomp 
  conda install trimmomatic
  conda install -c maxibor falcon
fi
#==============================================================================
# GET FALCON
#
if [[ "$GET_FALCON" -eq "1" ]]; then
  rm -fr falcon FALCON FALCON-* *.pl
  git clone https://github.com/pratas/falcon.git
  cd falcon/src/
  cmake .
  make
  cp FALCON ../../
  cp FALCON-FILTER ../../
  cp FALCON-EYE ../../
  cd ../../
  #cp falcon/scripts/*.pl .
fi
#==============================================================================
# GET GULL
#
if [[ "$GET_GULL" -eq "1" ]]; then
  rm -fr GULL/ GULL-map GULL-visual
  git clone https://github.com/pratas/GULL.git
  cd GULL/src/
  cmake .
  make
  cp GULL-map ../../
  cp GULL-visual ../../
  cd ../../
fi
#==============================================================================
# GET GOOSE
#
if [[ "$GET_GOOSE" -eq "1" ]]; then
  rm -fr goose/ goose-*
  git clone https://github.com/pratas/goose.git
  cd goose/src/
  make
  cp goose-* ../../
  cd ../../
fi
#==============================================================================
# GET KESTREL
#
if [[ "$GET_MAGNET" -eq "1" ]]; then
  rm -rf Magnet/ MAGNET
  git clone https://github.com/pratas/magnet.git
  cd Magnet/src/
  cmake .
  make
  cp Magnet ../../
  cd ../../
fi
#==============================================================================
# GET CHESTER
#
if [[ "$GET_CHESTER" -eq "1" ]]; then
  rm -fr chester*
  git clone https://github.com/pratas/chester.git
  cd chester/src/
  cmake .
  make
  cp CHESTER-filter ../../
  cp CHESTER-visual ../../
  cd ../../
fi
#==============================================================================
# GET AC
#
if [[ "$GET_AC" -eq "1" ]]; then
  git clone https://github.com/pratas/ac.git
  cd ac/src/
  cmake .
  make
  cp AC ../../
  cp AD ../../
  cd ../../
fi
#==============================================================================
# GET GECO
#
if [[ "$GET_GECO" -eq "1" ]]; then
  git clone https://github.com/pratas/geco.git
  cd geco/src/
  cmake .
  make
  cp GeCo ../../
  cp GeDe ../../
  cd ../../
fi
#==============================================================================
# GET CRYFA
#
if [[ "$GET_CRYFA" -eq "1" ]]; then
  git clone https://github.com/pratas/cryfa.git
  cd cryfa
  cmake .
  make
  cd ../
  mv cryfa cryfa_git
  cp cryfa_git/cryfa .
fi
#==============================================================================
# GET SMASHPP
#
if [[ "$GET_SMASHPP" -eq "1" ]]; then
  git clone https://pratas@github.com/smortezah/smashpp.git
  cd smashpp/
  cmake .
  make -j4
  cp smashpp ..
  cd ..
fi
#==============================================================================
# GET FRUIT
#
if [[ "$GET_FRUIT" -eq "1" ]]; then
  git clone https://pratas@github.com/smortezah/fruit.git
  cp fruit-map ..
  cp fruit-filter ..
  cp fruit-visual ..
  cd ..
fi
#==============================================================================
#==============================================================================
#==============================================================================


