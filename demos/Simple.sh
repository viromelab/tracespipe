#!/bin/bash
#
# INSTALLATION DEPENDENCIES:
#
# sudo apt-get install libopenblas-base
# conda install -c 'bioconda' art
# conda install -c 'bioconda' mummer4
#
RUN_ALL=1;
RUN_TRACES=1;
#
zcat MINIDB.fa.gz > ../src/VDB.fa
cd ../src/
#
if [[ "$RUN_ALL" -eq "1" ]];
  then
  ./TRACESPipe.sh --get-all-aux
  #
  gto_fasta_extract_read_by_pattern -p "AY386330.1" < VDB.fa > B19.fa
  gto_fasta_extract_read_by_pattern -p "X69198.1"   < VDB.fa > VARV.fa
  gto_fasta_extract_read_by_pattern -p "MG921180.1" < VDB.fa > HPV.fa
  #
  gto_fasta_mutate -s 0 -e 0.005 < mtDNA.fa > mtDNA_2.fa
  gto_fasta_mutate -s 0 -e 0.002 < B19.fa   > B19_2.fa
  gto_fasta_mutate -s 0 -e 0.003 < HPV.fa   > HPV_2.fa
  gto_fasta_mutate -s 0 -e 0.004 < VARV.fa  > VARV_2.fa
  #
  cp mtDNA_2.fa mtDNA_A.fa
  gto_fasta_mutate -s 0 -e 0.001 < mtDNA_2.fa > mtDNA_B.fa
  gto_fasta_mutate -s 0 -e 0.002 < mtDNA_2.fa > mtDNA_C.fa
  #
  cp B19_2.fa B19_A.fa
  gto_fasta_mutate -s 0 -e 0.001 < B19_2.fa > B19_B.fa
  gto_fasta_mutate -s 7 -e 0.001 < B19_2.fa > B19_C.fa
  #
  cp VARV_2.fa VARV_A.fa
  gto_fasta_mutate -s 0 -e 0.005 < VARV_2.fa > VARV_B.fa
  gto_fasta_mutate -s 7 -e 0.005 < VARV_2.fa > VARV_C.fa
  #
  cp HPV_2.fa HPV_A.fa
  gto_fasta_mutate -s 0 -e 0.01 < HPV_2.fa > HPV_B.fa
  gto_fasta_mutate -s 0 -e 0.01 < HPV_2.fa > HPV_C.fa
  #
  cat mtDNA_A.fa B19_A.fa VARV_A.fa HPV_A.fa > blood.fa
  cat mtDNA_B.fa B19_B.fa VARV_B.fa HPV_B.fa > brain.fa
  cat mtDNA_C.fa B19_C.fa VARV_C.fa HPV_C.fa > bone.fa
  #
  art_illumina -rs 0 -ss HS25 -sam -i blood.fa -p -l 150 -f 40 -m 200 -s 10 -o blood
  art_illumina -rs 0 -ss HS25 -sam -i brain.fa -p -l 150 -f 10 -m 200 -s 10 -o brain
  art_illumina -rs 0 -ss HS25 -sam -i bone.fa  -p -l 150 -f 30 -m 200 -s 10 -o bone
  #
  cp blood*.fq ../input_data
  cp brain*.fq ../input_data
  cp bone*.fq ../input_data
  #
  cd ../input_data
  rm -f blood*.fq.gz brain*.fq.gz bone*.fq.gz skin*.fq.gz teeth*.fq.gz kidney*.fq.gz lung*.fq.gz liver*.fq.gz heart*.fq.gz hair*.fq.gz
  gzip -f blood*.fq
  gzip -f brain*.fq
  gzip -f bone*.fq
  cd ../src
  #
  echo "blood:blood1.fq.gz:blood2.fq.gz"  > ../meta_data/meta_info.txt
  echo "brain:brain1.fq.gz:brain2.fq.gz" >> ../meta_data/meta_info.txt
  echo "bone:bone1.fq.gz:bone2.fq.gz"    >> ../meta_data/meta_info.txt
  #
fi
# ============================================================================
if [[ "$RUN_TRACES" -eq "1" ]];
  then
  ./TRACESPipe.sh --run-meta --inter-sim-size 5 --run-all-v-alig --run-mito --remove-dup --run-de-novo --run-hybrid --min-similarity 2 --view-top 5 --best-of-bests --very-sensitive --run-multiorgan-consensus
fi
#
# ============================================================================
