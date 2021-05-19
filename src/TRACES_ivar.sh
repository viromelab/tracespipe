#!/bin/bash
#
# ./TRACES_ivar.sh ../output_data/TRACES_viral_alignments/H_N-HV6B.fa ../output_data/TRACES_viral_alignments/viral_aligned_sorted-H_N-HV6B.bam H_N HV6B 0.2
#
Reference=$1;     # EXAMPLE: TTV.fa
Alignments=$2;    # EXAMPLE: ttv_aligned_sorted-heart.bam
Organ=$3;         # Example: heart
Label=$4;         # Example: TTV
Threshold=$5;     # Example: 0.10
#
echo "Using Reference   : $Reference";
echo "Using Alignments  : $Alignments";
echo "Using Organ       : $Organ";
echo "Using Viral Label : $Label";
echo "Threshold         : $Threshold";
#
rm -f $SID.gff3;
SID=`head $Reference | grep ">" | tr -d ">" | awk '{ print $1;}'`;
wget -O $SID.gff3 "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nucleotide&report=gff3&id=$SID" 2>> .tmp
samtools mpileup -aa -A -d 600000 -B -Q 0 $Alignments | ivar variants -p ../output_data/TRACES_results/Variants-$Organ-$Label -q 20 -t $Threshold -r $Reference -g $SID.gff3
cp $SID.gff3 ../output_data/TRACES_results/Gff-$Organ-$SID.gff3
#

