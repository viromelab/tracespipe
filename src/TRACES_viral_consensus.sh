#!/bin/bash
#
# moreFilters: http://samtools.github.io/bcftools/howtos/variant-calling.html
#
Reference=$1;     # EXAMPLE: TTV.fa
Alignments=$2;    # EXAMPLE: ttv_aligned_sorted-heart.bam
Organ=$3;         # Example: heart
Label=$4;         # Example: TTV
#
echo "Using Reference   : $Reference";
echo "Using Alignments  : $Alignments";
echo "Using Organ       : $Organ";
echo "Using Viral Label : $Label";
#
rm -f $Label-$Organ-calls.vcf.gz $Label-$Organ-calls.vcf.gz.csi $Label-$Organ-calls.norm.bcf $Label-$Organ-calls.norm.flt-indels.bcf $Label-$Organ-calls.norm.flt-indels.vcf.gz $Label-$Organ-calls.norm.flt-indels.vcf.gz.csi $Label-$Organ-calls.norm.vcf.gz;
#
# https://wikis.utexas.edu/display/bioiteam/Removing+duplicates+from+alignment+output
# Carefull with this filtering... [ambiguity]
#
# MASK LOW COVERAGE (<1) TO N
bedtools genomecov -ibam $Alignments -bga > $Label-coverage-$Organ.bed
awk '$4 < 1' $Label-coverage-$Organ.bed > $Label-zero-coverage-$Organ.bed 
#bedtools maskfasta -fi $Reference -bed $Label-$Organ-zero_coverage.bed -fo MASKED-$Label-consensus-$Organ.fa;
#
# CALLS 
samtools faidx $Reference # -P 9.9e-1                                         #here!
rawCalls="$Label-$Organ-raw.vcf.gz"
bcftools mpileup -Ov -f $Reference $Alignments | bcftools call --ploidy 1 -P 9.9e-1 -mv -Oz -o "$rawCalls"
bcftools index "$rawCalls"
#
# normalize indels
bcftools norm -f $Reference "$rawCalls" -Oz -o $Label-$Organ-calls.norm.vcf.gz
#
# filter adjacent indels within 5bp
bcftools filter --IndelGap 5 $Label-$Organ-calls.norm.vcf.gz -Oz -o $Label-$Organ-calls.norm.flt-indels.vcf.gz
# filter incompatible variants
finalVCF="$Label-$Organ-calls.vcf.gz"
./TRACES_filter_incompatible_variants.sh $Label-$Organ-calls.norm.filt-incompat.vcf.gz >| "$finalVCF";

#
# create bed file
zcat "$finalVCF" | vcf2bed --snvs > $Label-calls-$Organ.bed
#
# CONSENSUS
tabix -f "$finalVCF"
bcftools consensus -m $Label-zero-coverage-$Organ.bed -f $Reference "$finalVCF" > $Label-consensus-$Organ.fa
#
# Give new header name for the consensus sequence
tail -n +2 $Label-consensus-$Organ.fa > $Label-$Organ-TMP_FILE.xki
echo ">$Label consensus ($Organ)" > $Label-consensus-$Organ.fa
cat $Label-$Organ-TMP_FILE.xki >> $Label-consensus-$Organ.fa
rm -f $Label-$Organ-TMP_FILE.xki;
#
#
rm -f "$rawCalls" "$rawCalls.csi" $Label-$Organ-calls.norm.bcf $Label-$Organ-calls.norm.vcf.gz $Label-$Organ-calls.norm.flt-indels.bcf $Label-$Organ-calls.norm.flt-indels.vcf.gz "$finalVCL.csi";
#
#
