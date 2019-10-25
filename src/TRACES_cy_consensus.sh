#!/bin/bash
#
# moreFilters: http://samtools.github.io/bcftools/howtos/variant-calling.html
#
Reference=$1;     # EXAMPLE: cy.fa
Alignments=$2;    # EXAMPLE: cy_aligned_sorted-heart.bam
Organ=$3;         # Example: heart
#
echo "Using Reference   : $Reference";
echo "Using Alignments  : $Alignments";
echo "Using Organ       : $Organ";
#
# MASK LOW COVERAGE (<1) TO N
bedtools genomecov -ibam $Alignments -bga > cy-coverage-$Organ.bed
awk '$4 < 1' cy-coverage-$Organ.bed > cy-zero-coverage-$Organ.bed 
#bedtools maskfasta -fi $Reference -bed cy-zero-coverage-$Organ.bed -fo MASKED-cy-consensus-$Organ.fa;
#
# CALLS W CONSENSUS
samtools faidx $Reference # -P 9.9e-1
bcftools mpileup -Ou -f $Reference $Alignments | bcftools call --ploidy 1 -P 9.9e-1 -mv -Oz -o calls.vcf.gz
bcftools index calls.vcf.gz
#
# normalize indels
bcftools norm -f $Reference calls.vcf.gz -Oz -o calls.norm.vcf.gz
#
# filter adjacent indels within 5bp
bcftools filter --IndelGap 5 calls.norm.vcf.gz -Oz -o calls.norm.flt-indels.vcf.gz
#
# create bed file
zcat calls.norm.flt-indels.vcf.gz |vcf2bed --snvs > cy-calls-$Organ.bed
#
# CONSENSUS
tabix calls.norm.flt-indels.vcf.gz
bcftools consensus -m cy-zero-coverage-$Organ.bed -f $Reference calls.norm.flt-indels.vcf.gz > cy-consensus-$Organ.fa
#
# Give new header name for the consensus sequence
tail -n +2 cy-consensus-$Organ.fa > CY_TMP_FILE_$Organ.xki
echo "> $Organ CY consensus" > cy-consensus-$Organ.fa
cat CY_TMP_FILE_$Organ.xki >> cy-consensus-$Organ.fa
rm -f CY_TMP_FILE_$Organ.xki;
#
#
rm -f calls.vcf.gz calls.vcf.gz.csi calls.norm.bcf calls.norm.flt-indels.bcf calls.norm.flt-indels.vcf.gz calls.norm.vcf.gz;
#
#
