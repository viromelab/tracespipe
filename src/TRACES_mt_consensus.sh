#!/bin/bash
#
# moreFilters: http://samtools.github.io/bcftools/howtos/variant-calling.html
# MITOCHONDRIAL SNIPS DATABASE: http://www.mtdb.igp.uu.se/
#
Reference=$1;     # EXAMPLE: mtDNA.fa
Alignments=$2;    # EXAMPLE: mt_aligned_sorted-heart.bam
Organ=$3;         # Example: heart
#
echo "Using Reference   : $Reference";
echo "Using Alignments  : $Alignments";
echo "Using Organ       : $Organ";
#
# MASK LOW COVERAGE (<1) TO N
bedtools genomecov -ibam $Alignments -bga > mt-coverage-$Organ.bed
awk '$4 < 1' mt-coverage-$Organ.bed > mt-zero-coverage-$Organ.bed # CHANGE VALUE TO CHANGE MINIMUM OF DEPTH COVERAGE
#bedtools maskfasta -fi $Reference -bed zero_coverage.bed -fo MASKED-$LABEL-consensus-$Organ.fa;
#
# CALLS 
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
zcat calls.norm.flt-indels.vcf.gz |vcf2bed --snvs > mt-calls-$Organ.bed
#
# CONSENSUS
tabix -f calls.norm.flt-indels.vcf.gz
bcftools consensus -m mt-zero-coverage-$Organ.bed -f $Reference calls.norm.flt-indels.vcf.gz > mt-consensus-$Organ.fa
#
# Give new header name for the consensus sequence
tail -n +2 mt-consensus-$Organ.fa > MT_TMP_FILE_$Organ.xki
echo "> $Organ Mitochondrial DNA consensus" > mt-consensus-$Organ.fa
cat MT_TMP_FILE_$Organ.xki >> mt-consensus-$Organ.fa
rm -f MT_TMP_$Organ.xki;
#
#
rm -f calls.vcf.gz calls.vcf.gz.csi calls.norm.bcf calls.norm.flt-indels.bcf calls.norm.flt-indels.vcf.gz calls.norm.vcf.gz;
#
#
