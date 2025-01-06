#!/bin/bash

#This script is intended to pass over a vcf file and remove conflicting variants
#The choice is made first by higher quality, then by VCF order
#Confilicting alleles are snps at the same position as an insertion
#	and other variants which appear in the region covered by a deletion
#If there are multiple variants overlapping a position, the best one is retained
#	variants are processed in an order such that come conflicts may be resolved
#	Example two cariant sites covered by a deletion, where the deletion is lower
#	quality than either variant, in this way the maximum number of variants are
#	retained
#Written by Zachery Dickson, September 2024
#	Contact:	dicksoz@mcmaster.ca
#				zachery.dickson@helski.fi

if [ "$#" -lt 1 ]; then
	>&2 echo -e "Usage: $(basename $0) In.vcf[.gz] > Out.vcf.gz\n" \
		"\tFinds conflicting variants, selects the maximum conflict free set\n" \
		"\tAssumes variants are sorted by position";
	exit 1;
fi

file=$1;shift;
if file $file | grep -q gzip; then
	zcat $file;
else
	cat $file;
fi | 
	awk -F '\\t' -v call="$0 $file" '
		#j > i
		function testOverlap(i,j) {
			if(chrom[j] != chrom[i]){return 0;}
			if(left[j] > right[i]){return 0;}
			return 1;
		}
		BEGIN{OFS="\t"}
		/^##/{print;next;}
		/^#CHROM/{
			print "##incompatFilterCommand="call;
			print;
			next
		}
		{
			counter++;
			fullLine[counter] = $0;
			qual[counter] = $6;
			chrom[counter] = $1;
			left[counter] = $2;
			len = length($4)
#			if(length($5) > len){len = length($5)}
			right[counter] = $2+len-1;
			sortableLength[counter] = sprintf("%d%06d",len,counter);
		}
		END {
			#Process Variants in order of length
			n=asort(sortableLength);
			for(i = n; i >= 1; i--){
				#Iterate from longest to shortest 
				idx = substr(sortableLength[i],length(sortableLength[i])-5)+0
				if(Filtered[idx]){continue;}
				jdx = idx;
				bestQual=qual[idx];
				#Check all following entries as long as this variant is best,
				# we have not checked all variants, and the next one overlaps
				while(bestQual==qual[idx] && jdx < counter &&
						testOverlap(idx,jdx+1)){
					jdx++;
					if(qual[jdx] > qual[idx]){
						Filtered[idx] = 1;
						bestQual = qual[jdx];
					} else {
						Filtered[jdx] = 1;
					}
				}
			}
			for(idx=1;idx<=counter; idx++){
				if(!Filtered[idx]){
					print fullLine[idx];
				}
			}
		}
	' |
	bgzip;

