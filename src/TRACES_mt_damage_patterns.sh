#!/bin/bash
#
# ./TRACES_mt_damage_patterns.sh AB550331.1 
#
GID="$1"
ORGAN="$2";
THREADS=$3;
#
CHECK_ADAPTERS_AR () {
  if [ ! -f adapters_ar.fa ];
    then
    echo -e "\e[31mERROR: adapter sequences (adapters_ar.fa) not found!\e[0m"
    echo "TIP: before this, include valid adapters for AdapterRemoval"
    echo "For addition information, see the instructions at the web page."
    exit 1;
    fi
  }
#  
CHECK_ADAPTERS_AR;
#
echo -e "\e[34m[TRACESPipe]\e[32m Trimming, filtering, and collapsing with AdapterRemoval ...\e[0m";
AdapterRemoval --threads $THREADS --file1 FW_READS.fq.gz --file2 RV_READS.fq.gz --outputcollapsed reads.fq --trimns --trimqualities --minlength 30 --collapse --adapter-list adapters_ar.fa 1>> ../logs/Log-stdout-$ORGAN.txt 2>> ../logs/Log-stderr-$ORGAN.txt;
echo -e "\e[34m[TRACESPipe]\e[32m Aligning data using bwa ...\e[0m";
bwa index $GID.fa 1>> ../logs/Log-stdout-$ORGAN.txt 2>> ../logs/Log-stderr-$ORGAN.txt;
bwa mem -t $THREADS -I 0 -O 2 -N 0.02 -L 1024 -E 7 $GID.fa reads.fq > $GID-$ORGAN.sam 2>> ../logs/Log-stderr-$ORGAN.txt;
echo -e "\e[34m[TRACESPipe]\e[32m Adapting data with samtools ...\e[0m";
samtools view -bSh $GID-$ORGAN.sam > $GID-$ORGAN.bam
samtools view -bh -F4 $GID-$ORGAN.bam > FIL-$GID-$ORGAN.bam;
echo -e "\e[34m[TRACESPipe]\e[32m Estimating the damage of mtDNA using mapDamage2 ...\e[0m";
rm -fr ../output_data/TRACES_damage_$GID-$ORGAN
#mapDamage --rescale -d ../output_data/TRACES_mtdna_damage_$ORGAN_T -i FIL-$ORGAN_T.bam -r mtDNA.fa;
mapDamage -d ../output_data/TRACES_damage_$GID-$ORGAN -i FIL-$GID-$ORGAN.bam -r $GID.fa 1>> ../logs/Log-stdout-$ORGAN.txt 2>> ../logs/Log-stderr-$ORGAN.txt;
echo -e "\e[34m[TRACESPipe]\e[32m Done!\e[0m"
#
