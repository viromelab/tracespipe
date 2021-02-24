#!/bin/bash
#
efetch -db nucleotide -format fasta -id "MT682520" > sample_B19.fa
efetch -db nucleotide -format fasta -id "NC_007605" > sample_EBV.fa
cat sample_B19.fa sample_EBV.fa > sample_blood.fa
art_illumina -rs 0 -ss HS25 -sam -i sample_blood.fa -p -l 150 -f 20 -m 200 -s 10 -o blood
mv blood1.fq ../input_data/
mv blood2.fq ../input_data/
echo "sample_blood:blood1.fq:blood2.fq" > ../meta_data/input_data.txt
rm -f sample_blood.fa sample_B19.fa sample_EBV.fa blood1.aln blood2.aln blood.sam
#

