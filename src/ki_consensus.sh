#!/bin/bash
#
# moreFilters: http://samtools.github.io/bcftools/howtos/variant-calling.html
#
Reference=$1;     # mtDNA.fa
Alignments=$2;    # aligned_sorted-heart.bam
Organ=$3;         # organ name
#
rm -f calls.vcf.gz calls.vcf.gz.csi calls.norm.bcf calls.norm.flt-indels.bcf 
#
# call variants
bcftools mpileup -Ou -f $Reference $Alignments | bcftools call -mv -Oz -o calls.vcf.gz
bcftools index calls.vcf.gz
#
# normalize indels
bcftools norm -f $Reference calls.vcf.gz -Ob -o calls.norm.bcf
#
# filter adjacent indels within 5bp
bcftools filter --IndelGap 5 calls.norm.bcf -Ob -o calls.norm.flt-indels.bcf
#
bcftools index calls.norm.flt-indels.bcf
bcftools consensus -f mtDNA.fa calls.norm.flt-indels.bcf > consensus-$Organ.fa
#
