#!/bin/bash
#
# moreFilters: http://samtools.github.io/bcftools/howtos/variant-calling.html
# MITOCHONDRIAL SNIPS DATABASE: http://www.mtdb.igp.uu.se/
#
Reference=$1;     # EXAMPLE: mtDNA.fa
Alignments=$2;    # EXAMPLE: mt_aligned_sorted-heart.bam
Organ=$3;         # Example: heart
Label="mt";         # MT
#
echo "Using Reference   : $Reference";
echo "Using Alignments  : $Alignments";
echo "Using Organ       : $Organ";
#
rm -f $Label-$Organ-calls.vcf.g* $Label-$Organ-calls.vcf.gz.csi $Label-$Organ-calls.norm.bcf $Label-$Organ-calls.norm.flt-indels.bcf $Label-$Organ-calls.norm.flt-indels.vcf.gz $Label-$Organ-calls.norm.flt-indels.vcf.gz.csi $Label-$Organ-calls.norm.vcf.gz;
#
# MASK LOW COVERAGE (<1) TO N
#
bedtools genomecov -ibam $Alignments -bga > $Label-coverage-$Organ.bed
awk '$4 < 1' $Label-coverage-$Organ.bed > $Label-zero-coverage-$Organ.bed # CHANGE VALUE TO CHANGE MINIMUM OF DEPTH COVERAGE
#bedtools maskfasta -fi $Reference -bed zero_coverage.bed -fo MASKED-$LABEL-consensus-$Organ.fa;
#
# CALLS 
samtools faidx $Reference # -P 9.9e-1
bcftools mpileup -Ou -f $Reference $Alignments | bcftools call --ploidy 1 -P 9.9e-1 -mv -Oz -o $Label-$Organ-calls.vcf.gz
bcftools index $Label-$Organ-calls.vcf.gz
#
# normalize indels
bcftools norm -f $Reference $Label-$Organ-calls.vcf.gz -Oz -o $Label-$Organ-calls.norm.vcf.gz
#
# filter adjacent indels within 5bp
bcftools filter --IndelGap 5 $Label-$Organ-calls.norm.vcf.gz -Oz -o $Label-$Organ-calls.norm.flt-indels.vcf.gz
# filter incompatible variants
finalVCF="$Label-$Organ-calls.norm.flt-indels-incompat.vcf.gz"
./TRACES_filter_incompatible_variants.sh $Label-$Organ-calls.norm.flt-indels.vcf.gz >| "$finalVCF";
#
# create bed file
zcat "$finalVCF" |vcf2bed --snvs > $Label-calls-$Organ.bed
#
# CONSENSUS
tabix -f "$finalVCF"
bcftools consensus -m $Label-zero-coverage-$Organ.bed -f $Reference "$finalVCF" > $Label-consensus-$Organ.fa
#
# Give new header name for the consensus sequence
tail -n +2 $Label-consensus-$Organ.fa > $Label-TMP_FILE_$Organ.xki
echo ">Mitochondrial DNA consensus ($Organ)" > $Label-consensus-$Organ.fa
cat $Label-TMP_FILE_$Organ.xki >> $Label-consensus-$Organ.fa
rm -f $Label-TMP_FILE_$Organ.xki;
#
#
rm -f $Label-$Organ-calls.vcf.gz $Label-$Organ-calls.vcf.gz.csi $Label-$Organ-calls.norm.bcf $Label-$Organ-calls.norm.flt-indels.bcf $Label-$Organ-calls.norm.flt-indels.vcf.gz $Label-$Organ-calls.norm.flt-indels.vcf.gz.csi $Label-$Organ-calls.norm.vcf.gz "$finalVCF" "$finalVCF.csi" "${finalVCF/.vcf.gz/.bcf}";
#
#
