#!/bin/bash
gto_fasta_extract_read_by_pattern -p "AY386330.1" < VDB.fa > B19.fa
gto_fasta_extract_read_by_pattern -p "X69198.1" < VDB.fa > VARV.fa
#
gto_fasta_mutate -e 0.01 < B19.fa > B19_m1.fa
gto_fasta_mutate -e 0.05 < B19.fa > B19_m5.fa
gto_fasta_mutate -e 0.05 < VARV.fa > VARV_m1.fa
#
echo ">mutant1" > HEADER1;
echo ">mutant2" > HEADER2;
gto_fasta_extract_by_read -i 100 -e 200 < VARV_m1.fa | gto_fasta_to_seq > SAMPLE_VARV_m1.seq
gto_fasta_extract_by_read -i 100 -e 500 < VARV_m1.fa | gto_fasta_to_seq > SAMPLE_VARV_m1_2.seq
gto_fasta_extract_by_read -i 1 -e 3000 < B19_m1.fa | gto_fasta_to_seq > SAMPLE_B19_m1.seq
gto_fasta_extract_by_read -i 1 -e 1000 < B19_m5.fa | gto_fasta_to_seq > SAMPLE_B19_m5.seq
#
cat HEADER1 SAMPLE_B19_m1.seq SAMPLE_VARV_m1.seq > blood.fa
cat HEADER1 SAMPLE_B19_m5.seq SAMPLE_VARV_m1_2.seq > brain.fa
#
art_illumina -ss HS25 -sam -i blood.fa -p -l 150 -f 30 -m 200 -s 10 -o blood
art_illumina -ss HS25 -sam -i brain.fa -p -l 150 -f 30 -m 200 -s 10 -o brain
#
cp blood*.fq ../input_data
cp brain*.fq ../input_data
#
cd ../input_data
rm -f blood*.fq.gz brain*.fq.gz 
gzip blood*.fq
gzip brain*.fq
cd ../src
#
echo "blood:blood1.fq.gz:blood2.fq.gz" > ../meta_data/meta_info.txt
echo "brain:brain1.fq.gz:brain2.fq.gz" >> ../meta_data/meta_info.txt
#
