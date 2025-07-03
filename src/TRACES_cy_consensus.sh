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
rawCalls="raw-cy-$Organ.vcf.gz"
bcftools mpileup -Ov -f $Reference $Alignments | bcftools call --ploidy 1 -P 9.9e-1 -mv -Oz -o "$rawCalls"
bcftools index "$rawCalls"
#
# normalize indels
bcftools norm -f $Reference "$rawCalls" -Oz -o calls-cy-$Organ.norm.vcf.gz
#
# filter adjacent indels within 5bp
bcftools filter --IndelGap 5 calls-cy-$Organ.norm.vcf.gz -Oz -o calls-cy-$Organ.norm.flt-indels.vcf.gz
# filter incompatible variants
finalVCF="calls-cy-$Organ.vcf.gz"
./TRACES_filter_incompatible_variants.sh calls-cy-$Organ.norm.flt-indels.vcf.gz >| $finalVCF;
#
# create bed file
zcat "$finalVCF" | vcf2bed --snvs > cy-calls-$Organ.bed
#
# CONSENSUS
tabix -f "$finalVCF"
bcftools consensus -m cy-zero-coverage-$Organ.bed -f $Reference $finalVCF > cy-consensus-$Organ.fa
#
# Give new header name for the consensus sequence
tail -n +2 cy-consensus-$Organ.fa > CY_TMP_FILE_$Organ.xki
echo ">CY consensus $Organ" > cy-consensus-$Organ.fa
cat CY_TMP_FILE_$Organ.xki >> cy-consensus-$Organ.fa
rm -f CY_TMP_FILE_$Organ.xki;
#
#
rm -f "$rawCalls" "$rawCalls.csi" calls-cy-$Organ.norm.bcf calls-cy-$Organ.norm.flt-indels.bcf calls-cy-$Organ.norm.flt-indels.vcf.gz calls-cy-$Organ.norm.vcf.gz "$finalVCF" ${finalVCF/.vcf.gz/.bcf}; "$finalVCF.tbi"
#
#
