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
#  ./TRACESPipe.sh --get-all-aux
  #
  gto_fasta_extract_read_by_pattern -p "AY386330.1" < VDB.fa > B19.fa
  gto_fasta_extract_read_by_pattern -p "JN561323.2" < VDB.fa > HV2.fa
  gto_fasta_extract_read_by_pattern -p "X04370.1" < VDB.fa > HV3.fa
  gto_fasta_extract_read_by_pattern -p "DQ279927.1" < VDB.fa > HV4.fa
  gto_fasta_extract_read_by_pattern -p "AF148805.2" < VDB.fa > HV8.fa
  gto_fasta_extract_read_by_pattern -p "X69198.1" < VDB.fa > VARV.fa
  gto_fasta_extract_read_by_pattern -p "AB041963.1" < VDB.fa > TTV.fa
  gto_fasta_extract_read_by_pattern -p "MG921180.1" < VDB.fa > HPV.fa
  #
  gto_fasta_mutate -s 0 -e 0.01 < B19.fa > B19_2.fa
  gto_fasta_mutate -s 0 -e 0.20 < B19.fa > B19_20.fa
  gto_fasta_mutate -s 0 -e 0.01 < HV4.fa > HV4_2.fa
  gto_fasta_mutate -s 0 -e 0.01 < mtDNA.fa > mtDNA1.fa
  gto_fasta_mutate -s 0 -e 0.02 < mtDNA.fa > mtDNA2.fa
  gto_fasta_mutate -s 0 -e 0.03 < mtDNA.fa > mtDNA3.fa
  gto_fasta_mutate -s 0 -e 0.05 < mtDNA.fa > mtDNA5.fa
  gto_fasta_mutate -s 0 -e 0.05 < B19.fa > B19_5.fa
  gto_fasta_mutate -s 0 -e 0.10 < TTV.fa > TTV_6.fa
  gto_fasta_mutate -s 0 -e 0.15 < TTV.fa > TTV_7.fa
  gto_fasta_mutate -s 0 -e 0.10 < HPV.fa > HPV_2.fa
  gto_fasta_mutate -s 0 -e 0.05 < VARV.fa > VARV_5.fa
  #
  cat mtDNA.fa VARV.fa HV2.fa HV3.fa HV8.fa B19.fa > blood.fa
  #
  cp mtDNA.fa B-mtDNA-blood.fa;
  cp VARV.fa  B-VARV-blood.fa;
  cp HV2.fa   B-HV2-blood.fa;
  cp HV3.fa   B-HV3-blood.fa;
  cp HV8.fa   B-HV8-blood.fa;
  cp B19.fa   B-B19-blood.fa;
  #
  cat mtDNA1.fa VARV.fa HV4.fa B19.fa > brain.fa
  cp mtDNA1.fa B-mtDNA-brain.fa;
  cp VARV.fa   B-VARV-brain.fa;
  cp HV4.fa    B-HV4-brain.fa;
  cp B19.fa    B-B19-brain.fa;
  #
  cat mtDNA.fa B19_2.fa TTV_6.fa > bone.fa
  cp mtDNA.fa B-mtDNA-bone.fa;
  cp B19_2.fa B-B19-bone.fa;
  cp TTV_6.fa B-TTV-bone.fa;
  #
  cat mtDNA.fa B19.fa HV3.fa TTV.fa > skin.fa
  cp mtDNA.fa B-mtDNA-skin.fa;
  cp B19.fa B-B19-skin.fa;
  cp HV3.fa B-HV3-skin.fa;
  cp TTV.fa B-TTV-skin.fa;
  #
  cat mtDNA5.fa B19_5.fa TTV.fa HV2.fa > teeth.fa
  cp mtDNA5.fa B-mtDNA-teeth.fa;
  cp B19_5.fa B-B19-teeth.fa;
  cp TTV.fa B-TTV-teeth.fa;
  cp HV2.fa B-HV2-teeth.fa;
  #
  cat mtDNA1.fa HV2.fa TTV_7.fa B19.fa > kidney.fa
  cp mtDNA1.fa B-mtDNA-kidney.fa;
  cp HV2.fa B-HV2-kidney.fa;
  cp TTV_7.fa B-TTV-kidney.fa;
  cp B19.fa B-B19-kidney.fa;
  #
  cat mtDNA.fa HV4_2.fa VARV.fa > lung.fa
  cp mtDNA.fa B-mtDNA-lung.fa;
  cp HV4_2.fa B-HV4-lung.fa;
  cp VARV.fa B-VARV-lung.fa;
  #
  cat mtDNA2.fa HPV.fa VARV.fa HV4_2.fa > liver.fa
  cp mtDNA2.fa B-mtDNA-liver.fa;
  cp HPV.fa B-HPV-liver.fa;
  cp VARV.fa B-VARV-liver.fa;
  cp HV4_2.fa B-HV4-liver.fa;
  #
  cat mtDNA.fa HPV_2.fa B19_20.fa VARV_5.fa > heart.fa
  cp mtDNA.fa B-mtDNA-heart.fa;
  cp HPV_2.fa B-HPV-heart.fa;
  cp B19_20.fa B-B19-heart.fa;
  cp VARV_5.fa B-VARV-heart.fa;
  #
  cat mtDNA.fa HPV_2.fa HV4.fa > hair.fa
  cp mtDNA.fa B-mtDNA-hair.fa;
  cp HPV_2.fa B-HPV-hair.fa;
  cp HV4.fa B-HV4-hair.fa;
  #
  art_illumina -rs 0 -ss HS25 -sam -i blood.fa -p -l 150 -f 40 -m 200 -s 10 -o blood
  art_illumina -rs 0 -ss HS25 -sam -i brain.fa -p -l 150 -f 10 -m 200 -s 10 -o brain
  art_illumina -rs 0 -ss HS25 -sam -i bone.fa -p -l 150 -f 30 -m 200 -s 10 -o bone
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
  echo "blood:blood1.fq.gz:blood2.fq.gz" > ../meta_data/meta_info.txt
  echo "brain:brain1.fq.gz:brain2.fq.gz" >> ../meta_data/meta_info.txt
  echo "bone:bone1.fq.gz:bone2.fq.gz" >> ../meta_data/meta_info.txt
  #
fi
# ============================================================================
if [[ "$RUN_TRACES" -eq "1" ]];
  then
  ./TRACESPipe.sh --run-meta --inter-sim-size 5 --run-all-v-alig --run-mito --remove-dup --run-de-novo --run-hybrid --min-similarity 2 --view-top 5 --best-of-bests --very-sensitive --run-multiorgan-consensus 
# --run-de-novo-specific AF148805.2
fi
#
# ============================================================================
#
