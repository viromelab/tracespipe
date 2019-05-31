#!/bin/bash
#
# moreFilters: http://samtools.github.io/bcftools/howtos/variant-calling.html
#
Reference=$1;     # B19.fa
Alignments=$2;    # b19_aligned_sorted-heart.bam
Organ=$3;         # organ name
#
rm -f calls.vcf.gz calls.vcf.gz.csi calls.norm.bcf calls.norm.flt-indels.bcf calls.norm.flt-indels.vcf.gz calls.norm.vcf.gz 
#
# call variants
bcftools mpileup -Ou -f $Reference $Alignments | bcftools call --ploidy 1 -mv -Oz -o calls.vcf.gz
bcftools index calls.vcf.gz
#
# normalize indels
bcftools norm -f $Reference calls.vcf.gz -Oz -o calls.norm.vcf.gz
#
# filter adjacent indels within 5bp
bcftools filter --IndelGap 5 calls.norm.vcf.gz -Oz -o calls.norm.flt-indels.vcf.gz
#
# create consensus sequence
bcftools index calls.norm.flt-indels.vcf.gz
bcftools consensus -f $Reference calls.norm.flt-indels.vcf.gz > B19-consensus-$Organ.fa
tail -n +2 B19-consensus-$Organ.fa > TMP_FILE_X_KI.xki
echo "> $Organ B19 consensus" > B19-consensus-$Organ.fa
cat TMP_FILE_X_KI.xki >> B19-consensus-$Organ.fa
#
# create bed file
zcat calls.norm.flt-indels.vcf.gz |vcf2bed --snvs > B19-calls-$Organ.bed
#
