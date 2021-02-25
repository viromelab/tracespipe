#!/bin/bash
#
# ==============================================================================
#
echo "Trimmomatic    : `trimmomatic -version`";
echo "Cryfa          : `cryfa 2> VX.tmp ; grep "Cryfa v" VX.tmp | awk '{ print $2;}'`";
echo "MAGNET         : `MAGNET --version 2> VX.tmp | grep MAGNET VX.tmp | awk '{ print $3;}'`";
echo "FALCON         : `FALCON -V 2> VX.tmp | grep VERSION VX.tmp | awk '{ print $2;}'`";
echo "GTO            : `gto 2> VX.tmp | grep "GTO v" VX.tmp | awk '{ print $2;}'`";
echo "SPAdes         : `spades.py --version | awk '{ print $2;}'`";
echo "IGV            : [Graphical version]";
echo "Bowtie2        : `bowtie2 --version | head -n 1 | grep version | awk '{ print $3; }'`";
echo "samtools       : `samtools --version | head -n 1 | awk '{ print $2;}'`";
echo "bcftools       : `bcftools --version | head -n 1 | awk '{ print $2;}'`";
echo "Bedops         : `bedops --version | grep version | awk '{ print $2; }'`";
echo "Bedtools       : `bedtools --version | awk '{ print $2;}'`";
echo "fastq_pair     : [version not included]";
echo "efetch         : `efetch -version`";
echo "mapDamage      : `mapDamage --version`";
echo "Tabix          : `tabix 2> VX.tmp; grep "Version:" VX.tmp | awk '{ print $2; }'`";
echo "AdapterRemoval : `AdapterRemoval --version 2> VX.tmp; cat VX.tmp | awk '{ print $3;}'`";
echo "bwa            : `bwa 2> VX.tmp; grep Version VX.tmp | awk '{ print $2;}'`";
echo "art_illumina   : `art_illumina | grep Version | awk '{ print $3;}'`";
echo "blastn         : `blastn -version | head -n 1 | awk '{print $2;}'`";
echo "dnadiff        : `dnadiff --version 2> VX.tmp; grep version VX.tmp | awk '{ print $3; }'`";
#
rm -f VX.tmp;
#
# ==============================================================================

