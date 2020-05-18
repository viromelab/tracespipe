#!/bin/bash
#
# EXAMPLE: ./TRACES_multiorgan_consensus.sh B19 B19.fa B19-multiorgans.fa 8
#
# ==============================================================================
#
FIELD="$1";
REFERENCE="$2";
SCAFFOLDS="$3";
THREADS="$4";
#
# ==============================================================================
#
if [ ! -f "$2" ] || [ ! -f "$3" ];
  then
  echo -e "\e[38;5;208mWARNING: $2 or $3 files not found (multiorgan consensus)!\e[0m"
  echo "TIP: before this, run: ./TRACESPipe.sh --run-meta --run-all-v-alig --run-de-novo --run-hybrid"
  echo -e "\e[38;5;198mDon't panic! This warning may be cause by absence of previous results\e[0m."
  echo "For addition information, see the instructions at the web page."
  exit 1;
  fi
#
# ==============================================================================
#
rm -f $REFERENCE.* $FIELD-data_aligned.bam $FIELD-data_aligned_sorted.bam $FIELD-data_aligned_sorted.bam.bai
#
bwa index $REFERENCE
#
bwa mem -t $THREADS -I 0 -O 2 -N 0.02 -L 1024 -E 7 $REFERENCE $SCAFFOLDS > $FIELD-data_aligned.bam 
samtools sort $FIELD-data_aligned.bam > $FIELD-data_aligned_sorted.bam
samtools index $FIELD-data_aligned_sorted.bam $FIELD-data_aligned_sorted.bam.bai
#
bedtools genomecov -ibam $FIELD-data_aligned_sorted.bam -bga > $FIELD-coverage.bed
awk '$4 < 1' $FIELD-coverage.bed > $FIELD-zero-coverage.bed
samtools faidx $REFERENCE
echo "here"
bcftools mpileup -Ou -f $REFERENCE $FIELD-data_aligned_sorted.bam \
| bcftools call --ploidy 1 -mv -Oz -o $FIELD-calls.vcf.gz
bcftools index $FIELD-calls.vcf.gz
bcftools norm -f $REFERENCE $FIELD-calls.vcf.gz -Oz -o $FIELD-calls.norm.vcf.gz
zcat $FIELD-calls.norm.vcf.gz |vcf2bed --snvs > $FIELD-calls.bed
tabix -f $FIELD-calls.norm.vcf.gz
bcftools consensus -f $REFERENCE $FIELD-calls.norm.vcf.gz > $FIELD-consensus.fa
#
# Give new header name for the consensus sequence
tail -n +2 $FIELD-consensus.fa > $FIELD-TMP2_FILE.xki
echo ">$FIELD Multiorgan Consensus" > $FIELD-consensus.fa
cat $FIELD-TMP2_FILE.xki >> $FIELD-consensus.fa
rm -f $FIELD-TMP2_FILE.xki;
#
# ==============================================================================
