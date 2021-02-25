#!/bin/bash
#
echo -e "\e[34m[TRACESPipe]\e[32m Simulation of human ref VDB.fa and experimental samples (organ=sample_blood) ... \e[0m";
#
efetch -db nucleotide -format fasta -id "MT682520" > sample_B19.fa 2>> ../logs/Log-stderr-sample_blood.txt
efetch -db nucleotide -format fasta -id "NC_007605" > sample_EBV.fa 2>> ../logs/Log-stderr-sample_blood.txt
efetch -db nucleotide -format fasta -id "MT682522" > sample_MT.fa 2>> ../logs/Log-stderr-sample_blood.txt
cat sample_B19.fa sample_EBV.fa sample_MT.fa > sample_blood.fa 2>> ../logs/Log-stderr-sample_blood.txt
art_illumina -rs 0 -ss HS25 -i sample_blood.fa -p -l 150 -f 20 -m 200 -s 10 -o sample_blood 1>> ../logs/Log-stdout-sample_blood.txt 2>> ../logs/Log-stderr-sample_blood.txt
exit;
gzip sample_blood1.fq
gzip sample_blood2.fq
mv sample_blood1.fq.gz ../input_data/
mv sample_blood2.fq.gz ../input_data/
echo "sample_blood:sample_blood1.fq.gz:sample_blood2.fq.gz" > ../meta_data/meta_info.txt
rm -f sample_blood.fa sample_B19.fa sample_EBV.fa sample_blood1.aln sample_blood2.aln
# 
if test -f "VDB.fa"; then
  echo -e "\e[33mWARNING: VDB already found. Creating a copy to VDB.fa.tmp!\e[0m";
  mv VDB.fa VDB.fa.tmp;
fi
cp ../db_human_viral_ref/HRVDB.fasta VDB.fa
#
echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m";
