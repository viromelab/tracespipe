#!/bin/bash
#
declare -a INDIVIDUALS=('I1' 'I3' 'I4' 'I6' 'I7' 'I8' 'I12' 'I14' 'I15' 'I16' 'I18' 'I19' 'I24' 'I29' 'I30');
#
#
./TRACES.sh --gen-adapters
./TRACES.sh --get-phix
./TRACES.sh --get-mito
#./TRACES.sh --get-y-chromo
#./TRACES.sh --build-viral
#./TRACES.sh --build-unviral
#
for x in "${INDIVIDUALS[@]}"
  do
  echo "Building $x ...";
  mkdir -p ../../$x
  rm -fr ../../$x/*
  ln -s adapters.fa ../../$x/adapters.fa
  ln -s F_PHIX.fa ../../$x/F_PHIX.fa
  ln -s mtDNA.fa ../../$x/mtDNA.fa
  ln -s cy.fa ../../$x/cy.fa
  ln -s VDB.fa ../../$x/VDB.fa
  #ln -s DB.fa ../../$x/DB.fa
  cp TRACE*.sh ../../$x
  cp ../meta_data/$x-meta_info.txt ../../$x/meta_info.txt
  # 
  mapfile -t READS < ../../$x/meta_info.txt
  #
  for read in "${READS[@]}" # 
    do
    #
    forward=`echo $read | tr ':' '\t' | awk '{ print $2 }'`;
    reverse=`echo $read | tr ':' '\t' | awk '{ print $3 }'`;
    #
    ln -s /fas/NGS/pipes/fastq/fimm_ngs_mtoppinen/autopsy_NGS_novaseq1/$forward ../../$x/$forward
    ln -s /fas/NGS/pipes/fastq/fimm_ngs_mtoppinen/autopsy_NGS_novaseq1/$reverse ../../$x/$reverse
    #
    # FOR REMOTELY, RUN THIS:
    # scp dpratas@atlas.fimm.fi:/fas/NGS/pipes/fastq/fimm_ngs_mtoppinen/autopsy_NGS_novaseq1/V*.fa.gz .
    # ln -s $forward ../../$x/$forward
    # ln -s $reverse ../../$x/$reverse    
    #
    done
  #
  echo "Done!";
  done
