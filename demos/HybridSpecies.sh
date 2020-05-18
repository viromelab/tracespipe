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
  gto_fasta_extract_read_by_pattern -p "X69198.1" < VDB.fa > VARV.fa
  #
  gto_fasta_mutate -s 100 -e 0.01 < B19.fa > B19_m1.fa
  gto_fasta_mutate -s 100 -e 0.05 < B19.fa > B19_m5.fa
  gto_fasta_mutate -s 100 -e 0.01 < VARV.fa > VARV_m1.fa
  #
  gto_fasta_extract_by_read -i 1 -e 3000 < VARV_m1.fa | gto_fasta_to_seq > SAMPLE_VARV_m1.seq
  gto_fasta_extract_by_read -i 1 -e 4000 < VARV_m1.fa | gto_fasta_to_seq > SAMPLE_VARV_m2.seq
  gto_fasta_extract_by_read -i 1 -e 3000 < B19_m1.fa | gto_fasta_to_seq > SAMPLE_B19_m1.seq
  gto_fasta_extract_by_read -i 1 -e 1000 < B19_m5.fa | gto_fasta_to_seq > SAMPLE_B19_m5.seq
  #
  cat SAMPLE_B19_m1.seq SAMPLE_VARV_m1.seq | gto_fasta_from_seq -n "Mutant1" > Mutant1.fa
  cat SAMPLE_B19_m5.seq SAMPLE_VARV_m2.seq | gto_fasta_from_seq -n "Mutant2" > Mutant2.fa
  cat SAMPLE_B19_m1.seq SAMPLE_VARV_m2.seq | gto_fasta_from_seq -n "Mutant3" > Mutant3.fa
  #
  cat Mutant1.fa > blood.fa
  cat Mutant2.fa > bone.fa
  cat Mutant3.fa > brain.fa
  #
  art_illumina -rs 0 -ss HS25 -sam -i blood.fa -p -l 150 -f 30 -m 200 -s 10 -o blood
  art_illumina -rs 0 -ss HS25 -sam -i brain.fa -p -l 150 -f 30 -m 200 -s 10 -o brain
  art_illumina -rs 0 -ss HS25 -sam -i bone.fa  -p -l 150 -f 20 -m 200 -s 10 -o bone
  #
  cp blood*.fq ../input_data
  cp brain*.fq ../input_data
  cp bone*.fq  ../input_data
  #
  cd ../input_data
  rm -f blood*.fq.gz brain*.fq.gz bone*.fq.gz
  gzip blood*.fq
  gzip brain*.fq
  gzip bone*.fq
  cd ../src
  #
  echo "blood:blood1.fq.gz:blood2.fq.gz" >  ../meta_data/meta_info.txt
  echo "brain:brain1.fq.gz:brain2.fq.gz" >> ../meta_data/meta_info.txt
  echo "bone:bone1.fq.gz:bone2.fq.gz"    >> ../meta_data/meta_info.txt
  #
  # ============================================================================
  #
  ./TRACESPipe.sh --flush-output --flush-logs --run-meta --inter-sim-size 2 --run-all-v-alig --run-mito --remove-dup --run-de-novo --run-hybrid --min-similarity 1.5 --best-of-bests --very-sensitive
  #
  # ============================================================================
  #
  fi
# EVALUATION:
#
declare -a ORGANS=("blood" "bone" "brain");
declare -a VIRUSES=("B19" "VARV");
#
D_PATH="../output_data/TRACES_hybrid_R5_consensus";
#
printf "\nBood\n";
printf "B19\t";
cp $D_PATH/B19-consensus-blood.fa G_A.fa;
dnadiff G_A.fa Mutant1.fa ;
IDEN=`cat out.report | grep "AvgIdentity " | head -n 1 | awk '{ print $2;}'`;
ALBA=`cat out.report | grep "AlignedBases " | head -n 1 | awk '{ print $2;}'`;
SNPS=`cat out.report | grep TotalSNPs | awk '{ print $2;}'`;
printf "%s\t%s\t%s\n" "$ALBA" "$IDEN" "$SNPS";
rm -f G_A.fa
#
printf "VARV\t";
cp $D_PATH/VARV-consensus-blood.fa G_A.fa;
dnadiff G_A.fa Mutant1.fa ;
IDEN=`cat out.report | grep "AvgIdentity " | head -n 1 | awk '{ print $2;}'`;
ALBA=`cat out.report | grep "AlignedBases " | head -n 1 | awk '{ print $2;}'`;
SNPS=`cat out.report | grep TotalSNPs | awk '{ print $2;}'`;
printf "%s\t%s\t%s\n" "$ALBA" "$IDEN" "$SNPS";
rm -f G_A.fa
#
printf "\nBone\n";
printf "B19\t";
cp $D_PATH/B19-consensus-bone.fa G_A.fa;
dnadiff G_A.fa Mutant2.fa ;
IDEN=`cat out.report | grep "AvgIdentity " | head -n 1 | awk '{ print $2;}'`;
ALBA=`cat out.report | grep "AlignedBases " | head -n 1 | awk '{ print $2;}'`;
SNPS=`cat out.report | grep TotalSNPs | awk '{ print $2;}'`;
printf "%s\t%s\t%s\n" "$ALBA" "$IDEN" "$SNPS";
rm -f G_A.fa
#
printf "VARV\t";
cp $D_PATH/VARV-consensus-bone.fa G_A.fa;
dnadiff G_A.fa Mutant2.fa ;
IDEN=`cat out.report | grep "AvgIdentity " | head -n 1 | awk '{ print $2;}'`;
ALBA=`cat out.report | grep "AlignedBases " | head -n 1 | awk '{ print $2;}'`;
SNPS=`cat out.report | grep TotalSNPs | awk '{ print $2;}'`;
printf "%s\t%s\t%s\n" "$ALBA" "$IDEN" "$SNPS";
rm -f G_A.fa
#
printf "\nBrain\n";
printf "B19\t";
cp $D_PATH/B19-consensus-brain.fa G_A.fa;
dnadiff G_A.fa Mutant3.fa ;
IDEN=`cat out.report | grep "AvgIdentity " | head -n 1 | awk '{ print $2;}'`;
ALBA=`cat out.report | grep "AlignedBases " | head -n 1 | awk '{ print $2;}'`;
SNPS=`cat out.report | grep TotalSNPs | awk '{ print $2;}'`;
printf "%s\t%s\t%s\n" "$ALBA" "$IDEN" "$SNPS";
rm -f G_A.fa
#
printf "VARV\t";
cp $D_PATH/VARV-consensus-brain.fa G_A.fa;
dnadiff G_A.fa Mutant3.fa ;
IDEN=`cat out.report | grep "AvgIdentity " | head -n 1 | awk '{ print $2;}'`;
ALBA=`cat out.report | grep "AlignedBases " | head -n 1 | awk '{ print $2;}'`;
SNPS=`cat out.report | grep TotalSNPs | awk '{ print $2;}'`;
printf "%s\t%s\t%s\n" "$ALBA" "$IDEN" "$SNPS";
rm -f G_A.fa
#
# ==============================================================================
#
