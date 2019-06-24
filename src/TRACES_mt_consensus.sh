#!/bin/bash
#
# moreFilters: http://samtools.github.io/bcftools/howtos/variant-calling.html
# MITOCHONDRIAL SNIPS DATABASE: http://www.mtdb.igp.uu.se/
#
Reference=$1;     # mtDNA.fa
Alignments=$2;    # mt_aligned_sorted-heart.bam
Organ=$3;         # organ name
#
rm -f calls.vcf.gz calls.vcf.gz.csi calls.norm.bcf calls.norm.flt-indels.bcf 
#
# call variants
bcftools mpileup -Ou -f $Reference $Alignments | bcftools call --ploidy 1 -mv -Oz -o calls.vcf.gz
bcftools index calls.vcf.gz
#
# normalize indels
bcftools norm -f $Reference calls.vcf.gz -Oz -o calls.norm.vcf.gz
rm -f calls.vcf.gz;
#
# filter adjacent indels within 5bp
bcftools filter --IndelGap 5 calls.norm.vcf.gz -Oz -o calls.norm.flt-indels.vcf.gz
rm -f calls.norm.vcf.gz;
#
# create consensus sequence
bcftools index calls.norm.flt-indels.vcf.gz
#
bcftools consensus -f $Reference calls.norm.flt-indels.vcf.gz > mt-consensus-$Organ.fa
tail -n +2 mt-consensus-$Organ.fa > TMP_FILE_X_KI.xki
echo "> $Organ Mitochondrial consensus" > mt-consensus-$Organ.fa
cat TMP_FILE_X_KI.xki >> mt-consensus-$Organ.fa
rm -f TMP_FILE_X_KI.xki;
#
# create bed file
zcat calls.norm.flt-indels.vcf.gz |vcf2bed --snvs > mt-calls-$Organ.bed
rm -f calls.norm.flt-indels.vcf.gz;
#
