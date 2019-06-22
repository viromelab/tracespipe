#!/bin/bash
#
declare -a INDIVIDUALS=('I1' 'I3' 'I4' 'I6' 'I7' 'I8' 'I12' 'I14' 'I15' 'I16' 'I18' 'I19' 'I24' 'I29' 'I30');
#
#
./TRACESPipe.sh --gen-adapters
./TRACESPipe.sh --get-phix
./TRACESPipe.sh --get-mito
./TRACESPipe.sh --get-y-chromo
#./TRACESPipe.sh --build-viral
#./TRACESPipe.sh --build-unviral
#
for x in "${INDIVIDUALS[@]}"
  do
  echo "Building $x ...";
  mkdir -p ../../$x
  rm -fr ../../$x/*
  cp adapters.fa ../../$x/adapters.fa
  cp F_PHIX.fa ../../$x/F_PHIX.fa
  cp mtDNA.fa ../../$x/mtDNA.fa
  cp cy.fa ../../$x/cy.fa
  #cp VDB.fa ../../$x/VDB.fa
  #cp DB.fa ../../$x/DB.fa
  cp TRACE*.sh ../../$x
  cp ../meta_data/$x-meta_info.txt ../../$x/meta_info.txt
  # 
  scp dpratas@atlas.fimm.fi:/fas/NGS/pipes/fastq/fimm_ngs_mtoppinen/autopsy_NGS_novaseq1/V*.fa.gz .
  #
  mapfile -t READS < ../../$x/meta_info.txt
  #
  for read in "${READS[@]}" # 
    do
    #
    forward=`echo $read | tr ':' '\t' | awk '{ print $2 }'`;
    reverse=`echo $read | tr ':' '\t' | awk '{ print $3 }'`;
    #
    mv $forward ../../$x/$forward
    mv $reverse ../../$x/$reverse    
    #
    done
  #
  echo "Done!";
  done
