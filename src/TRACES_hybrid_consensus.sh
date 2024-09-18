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
rm -f $Label-$Organ-calls.vcf.gz $Label-$Organ-calls.vcf.gz.csi $Label-$Organ-calls.norm.bcf $Label-$Organ-calls.norm.vcf.gz $Label-$Organ-calls.norm.vcf.gz.csi;
#
# MASK LOW COVERAGE (<1) TO N
bedtools genomecov -ibam $Alignments -bga > $Label-coverage-$Organ.bed
awk '$4 < 1' $Label-coverage-$Organ.bed > $Label-zero-coverage-$Organ.bed 
#
# CALLS 
samtools faidx $Reference # -P 9.9e-1                                      # uses the new variations for consensus
bcftools mpileup -Ou -f $Reference $Alignments | bcftools call --ploidy 1 -M -mv -Oz -o $Label-$Organ-calls.vcf.gz
bcftools index $Label-$Organ-calls.vcf.gz
#
# normalize indels
bcftools norm -f $Reference $Label-$Organ-calls.vcf.gz -Oz -o $Label-$Organ-calls.norm.vcf.gz
#
# filter adjacent indels within 5bp
#bcftools filter --IndelGap 5 $Label-$Organ-calls.norm.vcf.gz -Oz -o $Label-$Organ-calls.norm.flt-indels.vcf.gz
# filter incompatible variants
finalVCF="$Label-$Organ-calls.norm.filt-incompat.vcf.gz"
./TRACES_filter_incompatible_variants.sh $Label-$Organ-calls.norm.vcf.gz >| "$finalVCF";
#
# create bed file
#zcat $Label-$Organ-calls.norm.flt-indels.vcf.gz |vcf2bed --snvs > $Label-calls-$Organ.bed
zcat "$finalVCF" | vcf2bed --snvs > $Label-calls-$Organ.bed
#
# CONSENSUS
tabix -f "$finalVCF"
bcftools consensus -m $Label-zero-coverage-$Organ.bed -f $Reference "$finalVCF" > $Label-consensus-$Organ.fa
#
# Give new header name for the consensus sequence
tail -n +2 $Label-consensus-$Organ.fa > $Label-$Organ-TMP_FILE.xki
echo ">Hybrid $Organ $Label consensus" > $Label-consensus-$Organ.fa
cat $Label-$Organ-TMP_FILE.xki >> $Label-consensus-$Organ.fa
rm -f $Label-$Organ-TMP_FILE.xki;
#
#
rm -f $Label-$Organ-calls.vcf.gz $Label-$Organ-calls.vcf.gz.csi $Label-$Organ-calls.norm.bcf $Label-$Organ-calls.norm.vcf.gz.csi $Label-$Organ-calls.norm.vcf.gz "$finalVCF" "$finalVCF.csi";
#
#
