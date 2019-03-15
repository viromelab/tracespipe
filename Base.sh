#!/bin/bash
#
GET_TOOLS=1;
#
GET_SMASHPP=1;
GET_FRUIT=1;
GET_FALCON=1;
GET_GULL=1;
GET_GOOSE=1;
GET_MAGNET=1;
GET_GECO=1;
GET_AC=1;
GET_CRYFA=1;
GET_GECO=1;
GET_CHESTER=1;
#
#
#==============================================================================
# GET TOOLS
#
if [[ "$GET_TOOLS" -eq "1" ]]; then
  sudo apt-get install terminator;
  sudo apt-get install git;
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
#
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
  cp falcon/scripts/*.pl .
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
# GET ... chester , ac, geco , cryfa

