#!/bin/bash
#
# sudo apt-get install libopenblas-base
# conda install -c 'bioconda' art
#
RUN_ALL=1;
zcat MINIDB.fa.gz > ../src/VDB.fa
cd ../src/
#
if [[ "$RUN_ALL" -eq "1" ]];
then
  ./TRACESPipe.sh --get-all-aux
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
  cp mtDNA.fa mtDNA-blood.fa;
  cp VARV.fa  VARV-blood.fa;
  cp HV2.fa   HV2-blood.fa;
  cp HV3.fa   HV3-blood.fa;
  cp HV8.fa   HV8-blood.fa;
  cp B19.fa   B19-blood.fa;
  #
  cat mtDNA1.fa VARV.fa HV4.fa B19.fa > brain.fa
  cp mtDNA1.fa mtDNA-brain.fa;
  cp VARV.fa   VARV-brain.fa;
  cp HV4.fa    HV4-brain.fa;
  cp B19.fa    B19-brain.fa;
  #
  cat mtDNA.fa B19_2.fa TTV_6.fa > bone.fa
  cp mtDNA.fa mtDNA-bone.fa;
  cp B19_2.fa B19-bone.fa;
  cp TTV_6.fa TTV-bone.fa;
  #
  cat mtDNA.fa B19.fa HV3.fa TTV.fa > skin.fa
  cp mtDNA.fa mtDNA-skin.fa;
  cp B19.fa B19-skin.fa;
  cp HV3.fa HV3-skin.fa;
  cp TTV.fa TTV-skin.fa;
  #
  cat mtDNA5.fa B19_5.fa TTV.fa HV2.fa > teeth.fa
  cp mtDNA5.fa mtDNA-teeth.fa;
  cp B19_5.fa B19-teeth.fa;
  cp TTV.fa TTV-teeth.fa;
  cp HV2.fa HV2-teeth.fa;
  #
  cat mtDNA1.fa HV2.fa TTV_7.fa B19.fa > kidney.fa
  cp mtDNA1.fa mtDNA-kidney.fa;
  cp HV2.fa HV2-kidney.fa;
  cp TTV_7.fa TTV-kidney.fa;
  cp B19.fa B19-kidney.fa;
  #
  cat mtDNA.fa HV4_2.fa VARV.fa > lung.fa
  cp mtDNA.fa mtDNA-lung.fa;
  cp HV4_2.fa HV4-lung.fa;
  cp VARV.fa VARV-lung.fa;
  #
  cat mtDNA2.fa HPV.fa VARV.fa HV4_2.fa > liver.fa
  cp mtDNA2.fa mtDNA-liver.fa;
  cp HPV.fa HPV-liver.fa;
  cp VARV.fa VARV-liver.fa;
  cp HV4_2.fa HV4-liver.fa;
  #
  cat mtDNA.fa HPV_2.fa B19_20.fa VARV_5.fa > heart.fa
  cp mtDNA.fa mtDNA-heart.fa;
  cp HPV_2.fa HPV-heart.fa;
  cp B19_20.fa B19-heart.fa;
  cp VARV_5.fa VARV-heart.fa;
  #
  cat mtDNA.fa HPV_2.fa HV4.fa > hair.fa
  cp mtDNA.fa mtDNA-hair.fa;
  cp HPV_2.fa HPV-hair.fa;
  cp HV4.fa HV4-hair.fa;
  #
  art_illumina -rs 0 -ss HS25 -sam -i blood.fa -p -l 150 -f 40 -m 200 -s 10 -o blood
  art_illumina -rs 0 -ss HS25 -sam -i brain.fa -p -l 150 -f 10 -m 200 -s 10 -o brain
  art_illumina -rs 0 -ss HS25 -sam -i bone.fa -p -l 150 -f 30 -m 200 -s 10 -o bone
  art_illumina -rs 0 -ss HS25 -sam -i skin.fa -p -l 150 -f 25 -m 200 -s 10 -o skin
  art_illumina -rs 0 -ss HS25 -sam -i teeth.fa -p -l 150 -f 30 -m 200 -s 10 -o teeth
  art_illumina -rs 0 -ss HS25 -sam -i kidney.fa -p -l 150 -f 20 -m 200 -s 10 -o kidney
  art_illumina -rs 0 -ss HS25 -sam -i lung.fa -p -l 150 -f 10 -m 200 -s 10 -o lung
  art_illumina -rs 0 -ss HS25 -sam -i liver.fa -p -l 150 -f 20 -m 200 -s 10 -o liver
  art_illumina -rs 0 -ss HS25 -sam -i heart.fa -p -l 150 -f 20 -m 200 -s 10 -o heart
  art_illumina -rs 0 -ss HS25 -sam -i hair.fa -p -l 150 -f 5 -m 200 -s 10 -o hair
  #
  cp blood*.fq ../input_data
  cp brain*.fq ../input_data
  cp bone*.fq ../input_data
  cp skin*.fq ../input_data
  cp teeth*.fq ../input_data
  cp kidney*.fq ../input_data
  cp lung*.fq ../input_data
  cp liver*.fq ../input_data
  cp heart*.fq ../input_data
  cp hair*.fq ../input_data
  #
  cd ../input_data
  rm -f blood*.fq.gz brain*.fq.gz bone*.fq.gz skin*.fq.gz teeth*.fq.gz kidney*.fq.gz lung*.fq.gz liver*.fq.gz heart*.fq.gz hair*.fq.gz
  gzip -f blood*.fq
  gzip -f brain*.fq
  gzip -f bone*.fq
  gzip -f skin*.fq
  gzip -f teeth*.fq
  gzip -f kidney*.fq
  gzip -f lung*.fq
  gzip -f liver*.fq
  gzip -f heart*.fq
  gzip -f hair*.fq
  cd ../src
  #
  echo "blood:blood1.fq.gz:blood2.fq.gz" > ../meta_data/meta_info.txt
  echo "brain:brain1.fq.gz:brain2.fq.gz" >> ../meta_data/meta_info.txt
  echo "bone:bone1.fq.gz:bone2.fq.gz" >> ../meta_data/meta_info.txt
  echo "skin:skin1.fq.gz:skin2.fq.gz" >> ../meta_data/meta_info.txt
  echo "teeth:teeth1.fq.gz:teeth2.fq.gz" >> ../meta_data/meta_info.txt
  echo "kidney:kidney1.fq.gz:kidney2.fq.gz" >> ../meta_data/meta_info.txt
  echo "lung:lung1.fq.gz:lung2.fq.gz" >> ../meta_data/meta_info.txt
  echo "liver:liver1.fq.gz:liver2.fq.gz" >> ../meta_data/meta_info.txt
  echo "heart:heart1.fq.gz:heart2.fq.gz" >> ../meta_data/meta_info.txt
  echo "hair:hair1.fq.gz:hair2.fq.gz" >> ../meta_data/meta_info.txt
  #
  # ============================================================================
  #
  ./TRACESPipe.sh --run-meta --inter-sim-size 5 --run-all-v-alig --run-mito --remove-dup --run-de-novo --run-hybrid --min-similarity 1 --view-top 10 --very-sensitive
  #
  # ============================================================================
  #
fi
#
# EVALUATION:
#
declare -a ORGANS=("blood" "bone" "brain" "hair" "heart" "kidney" "liver" "lung" "skin" "teeth");
declare -a VIRUSES=("B19" "HV2" "HV3" "HV4" "HV8" "HPV" "TTV" "VARV");
#
#D_PATH="../output_data/TRACES_hybrid_consensus";
D_PATH="../output_data/TRACES_hybrid_R2_consensus";
for organ in "${ORGANS[@]}"
  do
  printf "$organ\n";	  
  for virus in "${VIRUSES[@]}"
    do	  
    if [ -f $virus-$organ.fa ];
      then
      #echo "Diff -> $virus | $organ :";
      printf "$virus\t";	  
      cp $D_PATH/$virus-consensus-$organ.fa G_A.fa;
      cp $virus-$organ.fa G_B.fa
      dnadiff G_A.fa G_B.fa ;
      IDEN=`cat out.report | grep "AvgIdentity " | head -n 1 | awk '{ print $2;}'`; 
      SNPS=`cat out.report | grep TotalSNPs | awk '{ print $2;}'`;
      printf "$IDEN\t$SNPS\n";
      rm -f G_A.fa G_B.fa ;
      fi
    done
  #echo "Diff -> mtDNA | $organ :";
  printf "mtDNA\t";	  
  cp ../output_data/TRACES_mtdna_consensus/mt-consensus-$organ.fa G_A.fa;
  cp mtDNA-$organ.fa G_B.fa
  dnadiff G_A.fa G_B.fa ;
  IDEN=`cat out.report | grep "AvgIdentity " | head -n 1 | awk '{ print $2;}'`;
  SNPS=`cat out.report | grep TotalSNPs | awk '{ print $2;}'`;
  printf "$IDEN\t$SNPS\n\n";
  rm -f G_A.fa G_B.fa ;
  done
#
# ==============================================================================
#
