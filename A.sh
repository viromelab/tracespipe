# conda install -c bioconda trimmomatic

# INSTALL FALCON
git clone https://github.com/pratas/falcon.git
cd falcon/src/
cmake .
make
cp FALCON ../../
cp FALCON-FILTER ../../
cp FALCON-EYE ../../
cd ../../

# ==============================================================================
# ==============================================================================
#
# GET NCBI DATABASE FILE
rm -f assembly_summary.txt;
#
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/assembly_summary.txt
awk -F '\t' '{if($12=="Complete Genome") print $20}' assembly_summary.txt > ASCG.txt
#
# '{if($2=="Complete Genome"||$2=="Complete sequence") # <= IF YOU NEED ADDITIONAL FIELDS
#
# (awk -F '\t' '{print $12}' assembly_summary.txt | sort -d | uniq)
# KNOWN FIELDS:
# * assembly_level
# * Chromosome
# * Complete Genome
# * Scaffold
#
mkdir -p GB_DB_VIRAL/
rm -f GB_DB_VIRAL/*.fa.gz
#
cat ASCG.txt | xargs -I{} -n1 -P8 wget -P GB_DB_VIRAL {}/*_genomic.fna.gz
#
# MOVE CDS AND RNA SEQUENCES TO RESPECTIVE DIRS
mkdir -p GB_DB_VIRAL_CDS/
mkdir -p GB_DB_VIRAL_RNA/
rm -f GB_DB_VIRAL_CDS/*.fa.gz
rm -f GB_DB_VIRAL_RNA/*.fa.gz
mv GB_DB_VIRAL/*_cds_from_genomic.fna.gz GB_DB_VIRAL_CDS/
mv GB_DB_VIRAL/*_rna_from_genomic.fna.gz GB_DB_VIRAL_RNA/
# 
rm -f VDB.fa.gz;
zcat GB_DB_VIRAL/*.fna.gz | gzip -9 > VDB.fa.gz
#
# CHECK FOR ">" SYMBOLS INSIDE THE HEADER
zcat VDB.fa.gz | grep ">" | cut -c2- | grep ">"
#
# ==============================================================================
# ==============================================================================

#
# Needs a primer.fa
# $1 -> forward
# $2 -> reverse
#
trimmomatic PE -phred33 $1 $2 out1.fq out2.fq ILLUMINACLIP:primer.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35 CROP:75
# CHECK QUALITY WITH FASTQC
#
cat out1.fq out2.fq > reads.fq
#
time ./FALCON -F -t 1000 -v -l 47 -c 20 -n 8 -x top.txt merged.fq ../VDB.fa
