#!/bin/bash
gto_fasta_extract_read_by_pattern -p "AY386330.1" < VDB.fa > B19.fa
gto_fasta_extract_read_by_pattern -p "JN561323.2" < VDB.fa > HV2.fa
gto_fasta_extract_read_by_pattern -p "X04370.1" < VDB.fa > HV3.fa
gto_fasta_extract_read_by_pattern -p "DQ279927.1" < VDB.fa > HV4.fa
gto_fasta_extract_read_by_pattern -p "AF148805.2" < VDB.fa > HV8.fa
gto_fasta_extract_read_by_pattern -p "X69198.1" < VDB.fa > VARV.fa
gto_fasta_extract_read_by_pattern -p "AB041963.1" < VDB.fa > TTV.fa
gto_fasta_extract_read_by_pattern -p "MG921180.1" < VDB.fa > HPV.fa
#
gto_fasta_mutate -e 0.01 < B19.fa > B19_2.fa
gto_fasta_mutate -e 0.01 < HV4.fa > HV4_2.fa
gto_fasta_mutate -e 0.01 < mtDNA.fa > mtDNA2.fa
gto_fasta_mutate -e 0.05 < mtDNA.fa > mtDNA5.fa
gto_fasta_mutate -e 0.05 < B19.fa > B19_5.fa
gto_fasta_mutate -e 0.10 < TTV.fa > TTV_6.fa
gto_fasta_mutate -e 0.20 < TTV.fa > TTV_7.fa
gto_fasta_mutate -e 0.10 < HPV.fa > HPV_2.fa
#
cat mtDNA.fa VARV.fa HV2.fa HV3.fa HV8.fa B19.fa > blood.fa
cat mtDNA2.fa VARV.fa HV4.fa B19.fa cy.fa > brain.fa
cat mtDNA.fa B19_2.fa TTV_6.fa > bone.fa
cat mtDNA.fa B19.fa HV3.fa TTV.fa cy.fa > skin.fa
cat mtDNA5.fa B19_5.fa TTV.fa HV2.fa > teeth.fa
cat mtDNA.fa HV2.fa TTV_7.fa B19.fa > kidney.fa
cat mtDNA.fa HV4_2.fa VARV.fa > lung.fa
cat mtDNA.fa HPV.fa VARV.fa HV4_2.fa > liver.fa
cat mtDNA.fa HPV_2.fa > heart.fa
cat mtDNA.fa HPV_2.fa HV4.fa > hair.fa
#
art_illumina -ss HS25 -sam -i blood.fa -p -l 150 -f 40 -m 200 -s 10 -o blood
art_illumina -ss HS25 -sam -i brain.fa -p -l 150 -f 10 -m 200 -s 10 -o brain
art_illumina -ss HS25 -sam -i bone.fa -p -l 150 -f 30 -m 200 -s 10 -o bone
art_illumina -ss HS25 -sam -i skin.fa -p -l 150 -f 25 -m 200 -s 10 -o skin
art_illumina -ss HS25 -sam -i teeth.fa -p -l 150 -f 30 -m 200 -s 10 -o teeth
art_illumina -ss HS25 -sam -i kidney.fa -p -l 150 -f 20 -m 200 -s 10 -o kidney
art_illumina -ss HS25 -sam -i lung.fa -p -l 150 -f 10 -m 200 -s 10 -o lung
art_illumina -ss HS25 -sam -i liver.fa -p -l 150 -f 20 -m 200 -s 10 -o liver
art_illumina -ss HS25 -sam -i heart.fa -p -l 150 -f 20 -m 200 -s 10 -o heart
art_illumina -ss HS25 -sam -i hair.fa -p -l 150 -f 5 -m 200 -s 10 -o hair
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
gzip blood*.fq
gzip brain*.fq
gzip bone*.fq
gzip skin*.fq
gzip teeth*.fq
gzip kidney*.fq
gzip lung*.fq
gzip liver*.fq
gzip heart*.fq
gzip hair*.fq
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
