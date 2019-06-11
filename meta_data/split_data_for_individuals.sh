#!/bin/bash
#
declare -a INDIVIDUALS=('I1' 'I3' 'I4' 'I6' 'I7' 'I8' 'I12' 'I14' 'I15' 'I16' 'I18' 'I19' 'I24' 'I29' 'I30');
#
#
./TRACES.sh --gen-adapters
./TRACES.sh --get-phix
./TRACES.sh --get-mito
./TRACES.sh --get-y-chromo
#./TRACES.sh --build-viral
#./TRACES.sh --build-unviral
#
for x in "${INDIVIDUALS[@]}"
  do
  echo "Building $x ...";
  mkdir -p ../../$x
  rm -fr ../../$x/*
  cp adapters.fa ../../$x
  cp F_PHIX.fa ../../$x
  cp mtDNA.fa ../../$x
  cp cy.fa ../../$x
  cp VDB.fa ../../$x
  #cp DB.fa ../../$x
  cp TRACE*.sh ../../$x
  cp ../meta_data/$x-meta_info.txt meta_info.txt
  echo "Done!";
  done

