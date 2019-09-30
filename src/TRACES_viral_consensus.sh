#!/bin/bash
#
# moreFilters: http://samtools.github.io/bcftools/howtos/variant-calling.html
#
Reference=$1;     # EXAMPLE: TTV.fa
Alignments=$2;    # EXAMPLE: ttv_aligned_sorted-heart.bam
Organ=$3;         # Example: heart
LABEL=$4;         # Example: TTV
#

rm -f calls.vcf.gz calls.vcf.gz.csi calls.norm.bcf calls.norm.flt-indels.bcf calls.norm.flt-indels.vcf.gz calls.norm.vcf.gz calls.norm.flt-indels.vcf.gz.csi

# MASK LOW COVERAGE (<1) TO N
bedtools genomecov -ibam $Alignments -bga > coverage.bed
awk '$4 < 1' coverage.bed > zero_coverage.bed # CHANGE VALUE TO CHANGE MINIMUM OF DEPTH COVERAGE
#bedtools maskfasta -fi $Reference -bed zero_coverage.bed -fo MASKED-$LABEL-consensus-$Organ.fa;

# CALLS W CONSENSUS
samtools faidx $Reference # -P 9.9e-1
samtools mpileup -Ou -f $Reference $Alignments | bcftools call --ploidy 1 -P 9e-1 -mv -Oz -o calls.vcf.gz
bcftools index calls.vcf.gz

# normalize indels
bcftools norm -f $Reference calls.vcf.gz -Oz -o calls.norm.vcf.gz

# filter adjacent indels within 5bp
bcftools filter --IndelGap 5 calls.norm.vcf.gz -Oz -o calls.norm.flt-indels.vcf.gz

# CONSENSUS
tabix calls.norm.flt-indels.vcf.gz
bcftools consensus -m zero_coverage.bed -f $Reference calls.norm.flt-indels.vcf.gz > $LABEL-consensus-$Organ.fa

# Give new header name for the consensus sequence
tail -n +2 $LABEL-consensus-$Organ.fa > TMP_FILE_X_KI.xki
echo "> $Organ $LABEL consensus" > $LABEL-consensus-$Organ.fa
cat TMP_FILE_X_KI.xki >> $LABEL-consensus-$Organ.fa
rm -f TMP_FILE_X_KI.xki;

#TODO: falta criar o bed ...

exit 1;

# call variants
bcftools mpileup -Ou -f $Reference $Alignments | bcftools call --ploidy 1 -mv -Oz -o calls.vcf.gz
bcftools index calls.vcf.gz

# normalize indels
bcftools norm -f $Reference calls.vcf.gz -Oz -o calls.norm.vcf.gz

# filter adjacent indels within 5bp
bcftools filter --IndelGap 5 calls.norm.vcf.gz -Oz -o calls.norm.flt-indels.vcf.gz

# create consensus index sequence
bcftools index calls.norm.flt-indels.vcf.gz

# create consensus sequence Using masker for N's (from bed)
bcftools consensus -f $Reference calls.norm.flt-indels.vcf.gz > $LABEL-consensus-$Organ.fa

# MASK NNN
bedtools genomecov -ibam $Alignments -bga > coverage.bed
awk '$4 < 1' coverage.bed > zero_coverage.bed # CHANGE VALUE TO CHANGE MINIMUM OF DEPTH COVERAGE
bedtools maskfasta -fi $LABEL-consensus-$Organ.fa -bed zero_coverage.bed -fo MASKED-$LABEL-consensus-$Organ.fa;

# Give new header name for the consensus sequence
tail -n +2 MASKED-$LABEL-consensus-$Organ.fa > TMP_FILE_X_KI.xki
echo "> $Organ $LABEL consensus" > $LABEL-consensus-$Organ.fa
cat TMP_FILE_X_KI.xki >> $LABEL-consensus-$Organ.fa
rm -f TMP_FILE_X_KI.xki;
#
rm -f calls.vcf.gz calls.vcf.gz.csi calls.norm.bcf calls.norm.flt-indels.bcf calls.norm.flt-indels.vcf.gz calls.norm.vcf.gz 
