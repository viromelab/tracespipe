#!/bin/bash
#
VIRUS="$3";
ORGAN="$4";
THREADS="$5";
#
if [ ! -f "$1" ] || [ ! -f "$2" ];
  then
  echo -e "\e[38;5;208mWARNING: $1 or $2 file not found!\e[0m"
  echo -e "\e[38;5;198mDon't panic! This warning may be cause by absence of previous results\e[0m."
  echo "For addition information, see the instructions at the web page."
  exit 1;
  fi
#
rm -f $1.* 
bwa index $1
bwa mem -t $THREADS -I 0 -O 2 -N 0.03 -L 1024 -E 7 $1 $2 > specific_aligned_$VIRUS-$ORGAN.sam 
samtools sort specific_aligned_$VIRUS-$ORGAN.sam > specific_aligned_sorted_$VIRUS-$ORGAN.bam
samtools index specific_aligned_sorted_$VIRUS-$ORGAN.bam specific_aligned_sorted_$VIRUS-$ORGAN.bam.bai
rm -f specific_aligned_$VIRUS-$ORGAN.sam;
#
bedtools genomecov -ibam specific_aligned_sorted_$VIRUS-$ORGAN.bam -bga > $VIRUS-coverage-$ORGAN.bed
awk '$4 < 1' $VIRUS-coverage-$ORGAN.bed > $VIRUS-zero-coverage-$ORGAN.bed
samtools faidx $1 # -P 9.9e-1                                      # here!
bcftools mpileup -Ou -f $1 specific_aligned_sorted_$VIRUS-$ORGAN.bam \
| bcftools call --ploidy 1 -mv -Oz -o $VIRUS-$ORGAN-calls.vcf.gz
bcftools index $VIRUS-$ORGAN-calls.vcf.gz
bcftools norm -f $1 $VIRUS-$ORGAN-calls.vcf.gz -Oz -o $VIRUS-$ORGAN-calls.norm.vcf.gz
zcat $VIRUS-$ORGAN-calls.norm.vcf.gz |vcf2bed --snvs > $VIRUS-calls-$ORGAN.bed
tabix -f $VIRUS-$ORGAN-calls.norm.vcf.gz
bcftools consensus -m $VIRUS-zero-coverage-$ORGAN.bed -f $1 $VIRUS-$ORGAN-calls.norm.vcf.gz > $VIRUS-consensus-$ORGAN.fa
#
# Give new header name for the consensus sequence
tail -n +2 $VIRUS-consensus-$ORGAN.fa > $VIRUS-$ORGAN-TMP2_FILE.xki
echo ">Denovo specific consensus $VIRUS (organ=$ORGAN)" > $VIRUS-consensus-$ORGAN.fa
cat $VIRUS-$ORGAN-TMP2_FILE.xki >> $VIRUS-consensus-$ORGAN.fa
rm -f $VIRUS-$ORGAN-TMP2_FILE.xki;
#
rm -f $VIRUS-$ORGAN-calls.vcf.gz $VIRUS-$ORGAN-calls.vcf.gz.csi $VIRUS-$ORGAN-calls.norm.bcf $VIRUS-$ORGAN-calls.norm.vcf.gz.csi $VIRUS-$ORGAN-calls.norm.vcf.gz;
#
