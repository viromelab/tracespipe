<div align="center">
  
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)
[![Speed](https://img.shields.io/static/v1.svg?label=Ultra-Fast&message=High%20speed%20performance&color=green)](#)
[![Release](https://img.shields.io/static/v1.svg?label=Release&message=v1.1.1&color=orange)](#)
[![Conda](https://img.shields.io/static/v1.svg?label=Conda&message=Bioconda&color=green)](#)
[![Conda](https://img.shields.io/static/v1.svg?label=Conda&message=Cobilab&color=green)](#)
[![TinyURL](https://img.shields.io/static/v1.svg?label=TinyURL&message=traces-pipe&color=blue)](https://tinyurl.com/traces-pipe)
<!--[![Build Status](https://travis-ci.org/pratas/traces.svg?branch=master)](https://travis-ci.org/pratas/traces)-->

</div>
<br>
<p align="center">
<img src="imgs/logo.png" alt="TRACES Pipeline" height="200" border="0" /><br>
<i>A next-generation sequencing pipeline for identification, assembly,<br> and analysis of viral and human-host genomes at multi-organ level.</i>
<br><br>

## 1. About ##


TRACESPipe is a next-generation sequencing pipeline for identification, assembly, and analysis of viral and human-host genomes at multi-organ level. The identification and assembly of viral genomes rely on cooperation between three modalities:
<ul>
<li>compression-based predictors;</li>
<li>sequence alignments;</li>
<li><i>de-novo</i> assembly.</li>
</ul>
The compression-based prediction applies FALCON technology with ultra-fast comparative quantification to find the best reference genome (from a large viral database) containing the highest similarity relative to the sequenced reads. After identification, the reads are aligned according to the best reference by Bowtie2. A consensus sequence is produced with specific filters using Bcftools. Then, <i>de-novo</i> assembly (SPAdes) is involved in building scaffolds. The high coverage scaffolds that overlap totally or partially the consensus sequence are used to validate or either augment the new genome. The final analysis of the assembly is interactively supervised with the IGV with the goal of drafting the final sequence.

For the human-host variant call identification, the same procedure is followed although directly starting within the second point, given the use of the same reference (Cambridge Reference) to all the cases.

<br>
<p align="center">
<img src="imgs/pipeline.png" alt="TRACESPipe architecture" height="500" border="0" />
</p>
<br>

The previous image shows the architecture of TRACESPipe, where the green line stands for the mitochondrial human line. This pipeline has been tested in Illumina HiSeq and NovaSeq platforms. The operating system required to run it is Linux. In windows use cygwin (https://www.cygwin.com/) and make sure that it is included in the installation: cmake, make, zcat, unzip, wget, tr, grep (and any dependencies). If you install the complete cygwin packet then all these will be installed. After, all steps will be the same as in Linux.


## 2. Installation and Configuration ##

### 2.1 Installation ###

To install TRACESPipe, run the following commands in a Linux OS:
```
git clone https://github.com/pratas/traces.git
cd traces/src/
chmod +x TRACES*.sh
./TRACESPipe.sh --install
```
<img src="imgs/conda_logo.png" alt="CONDA" height="14" border="0" /> is needed for installation.<br>
To install Conda use the following steps:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
Additional instructions can be found here:
```
https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html
```

### 2.2 Configuration ###

To configure TRACESPipe add your FASTQ files gziped at folder
```
input_data/
```
Then, add a file exclusively with name <b>meta\_info.txt</b> at folder
```
meta_data/
```
This file needs to specify the organ type (with a single word name) and the filenames for the paired end reads. An example of the content of meta\_info.txt is the following:
```
skin:V1_S44_R1_001.fastq.gz:V1_S44_R2_001.fastq.gz
brain:V2_S29_R1_001.fastq.gz:V2_S29_R2_001.fastq.gz
colon:V3_S45_R1_001.fastq.gz:V3_S45_R2_001.fastq.gz
```

## 3. Running ##

To run TRACES Pipeline, use the following command:
```
./TRACESPipe.sh <parameters>
```
There are many parameters and configurations that can be used.
See the next section for more information about the usage.

## 4. Usage ##

```
./TRACESPipe.sh -h
```

```                                                     
                                                                
         ████████╗ ██████╗   █████╗   ██████╗ ███████╗ ███████╗   
         ╚══██╔══╝ ██╔══██╗ ██╔══██╗ ██╔════╝ ██╔════╝ ██╔════╝   
            ██║    ██████╔╝ ███████║ ██║      █████╗   ███████╗   
            ██║    ██╔══██╗ ██╔══██║ ██║      ██╔══╝   ╚════██║   
            ██║    ██║  ██║ ██║  ██║ ╚██████╗ ███████╗ ███████║   
            ╚═╝    ╚═╝  ╚═╝ ╚═╝  ╚═╝  ╚═════╝ ╚══════╝ ╚══════╝   
                                                                  
                             P I P E L I N E                            
                                                                
    A Next-generation sequencing pipeline for identification, assembly,
    and analysis of viral and human-host genomes at a multi-organ level.
                                                                
    Usage: ./TRACESPipe.sh [options]                             
                                                                
    -h,    --help            Show this help message and exit,     
    -v,    --version         Show the version and some information,  
                                                                  
    -i,    --install         Installation of all the tools,       
                                                                  
    -vdb,  --build-viral     Build viral database (all sequences), 
    -vdbr, --build-viral-r   Build viral database (references only),  
    -udb,  --build-unviral   Build non viral database (control),  
                                                                  
    -gx,   --get-extra-vir   Downloads/appends (VDB) extra viral seq, 
    -gad,  --gen-adapters    Generate FASTA file with adapters,   
    -gp,   --get-phix        Extracts PhiX genomes (Needs viral DB),  
    -gm,   --get-mito        Downloads human Mitochondrial genome,
    -gy,   --get-y-chromo    Downloads human Y-chromosome,        
                                                                  
    -rm,   --run-meta        Run viral metagenomic identification,    
    -ro,   --run-meta-nv     Run NON-viral metagenomic identification,   
                                                                  
    -rava, --run-all-v-alig  Run all viral align/sort/consensus seqs,    
                                                                 
    -rb,   --run-b19         Run B19   align and consensus seq,    
    -rh1,  --run-hv1         Run HV1   align and consensus seq,    
    -rh2,  --run-hv2         Run HV2   align and consensus seq,    
    -rh3,  --run-hv3         Run HV3   align and consensus seq,    
    -rh4,  --run-hv4         Run HV4   align and consensus seq,    
    -rh5,  --run-hv5         Run HV5   align and consensus seq,    
    -rh6,  --run-hv6         Run HV6   align and consensus seq,    
    -rh6a, --run-hv6a        Run HV6A  align and consensus seq,    
    -rh6b, --run-hv6b        Run HV6B  align and consensus seq,    
    -rh7,  --run-hv7         Run HV7   align and consensus seq,    
    -rh8,  --run-hv8         Run HV8   align and consensus seq,    
    -rh8,  --run-hv8         Run HV8   align and consensus seq,    
    -rp1,  --run-poly1       Run Polyoma 1  align and consensus seq,  
    -rp2,  --run-poly2       Run Polyoma 2  align and consensus seq,  
    -rp3,  --run-poly3       Run Polyoma 3  align and consensus seq,  
    -rp4,  --run-poly4       Run Polyoma 4  align and consensus seq,  
    -rp5,  --run-poly5       Run Polyoma 5  align and consensus seq,  
    -rp6,  --run-poly6       Run Polyoma 6  align and consensus seq,  
    -rp7,  --run-poly7       Run Polyoma 7  align and consensus seq,  
    -rp8,  --run-poly8       Run Polyoma 8  align and consensus seq,  
    -rp9,  --run-poly9       Run Polyoma 9  align and consensus seq,  
    -rp10, --run-poly10      Run Polyoma 10 align and consensus seq,  
    -rp11, --run-poly11      Run Polyoma 11 align and consensus seq,  
    -rp12, --run-poly12      Run Polyoma 12 align and consensus seq,  
    -rp13, --run-poly13      Run Polyoma 13 align and consensus seq,  
    -rp14, --run-poly14      Run Polyoma 14 align and consensus seq,  
    -rtt,  --run-ttv         Run TTV   align and consensus seq,    
    -rbv1, --run-hbov1       Run HBoV1 align and consensus seq,    
    -rbv0, --run-hbovnot1    Run HBoV (2,3,...) align/consensus seq, 
    -rhbv, --run-hbv         Run HBV   align and consensus seq,    
    -rhpv, --run-hpv         Run HPV   align and consensus seq,    
    -rvar, --run-varv        Run VARV  align and consensus seq,    
                                                                 
    -rsr,  --run-specific    Run specific REF align/consensus seq, 
                                                                 
    -rmt,  --run-mito        Run Mito  align and consensus seq,   
    -rcy,  --run-y-chromo    Run CY    align and consensus seq,    
                                                                  
    -rda,  --run-de-novo     Run de-novo assembly,               
                                                                 
    -ra,   --run-analysis    Run data analysis,                   
    -all,  --run-all         Run all the options.                 
                                                                
    Example: ./TRACESPipe.sh --run-meta --run-b19 --run-mito 
                                                                
    meta_info.txt -> 'organ:reads_forward.fa.gz:reads_reverse.fa.gz'  
    The reads and meta_info.txt must be in the src/ folder.  
 
```
## 5. Examples ##

### 5.1 Building a Parvovirus consensus sequence (if exists in the FASTQ samples): ###
```
./TRACESPipe.sh --run-meta --run-b19
```
The output consensus sequence is included at 
```
output_data/TRACES_viral_consensus
```
while the alignments at
```
output_data/TRACES_viral_alignments
```
and the BED files at
```
output_data/TRACES_viral_bed
```

### 5.2 Building a mitochondrial consensus sequence (if exists in the FASTQ samples): ###

```
./TRACESPipe.sh --run-mito
```
The output consensus sequence is included at
```
output_data/TRACES_mtdna_consensus
```
while the alignments at
```
output_data/TRACES_mtdna_alignments
```
and the BED files at
```
output_data/TRACES_mtdna_bed
```

## 6. Programs ##

TRACES Pipeline uses a combination of the following tools:

| Tool | URL | Article |
| --- | --- | --- |
| &#x1F49A;&nbsp; Cryfa | [[https://github.com/cobilab/cryfa]](https://github.com/cobilab/cryfa) | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://doi.org/10.1093/bioinformatics/bty645) |
| &#x1F49A;&nbsp; GTO | [[https://github.com/cobilab/gto]](https://github.com/cobilab/gto) | 
| &#x1F49A;&nbsp; Trimmomatic | [[http://www.usadellab.org/cms/?page=trimmomatic]](http://www.usadellab.org/cms/?page=trimmomatic) | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://academic.oup.com/bioinformatics/article/30/15/2114/2390096) |
| &#x1F49A;&nbsp; MAGNET | [[https://github.com/cobilab/magnet]](https://github.com/cobilab/magnet) | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://www.eurasip.org/Proceedings/Eusipco/Eusipco2018/papers/1570439333.pdf) |
| &#x1F49A;&nbsp; FALCON-meta | [[https://github.com/cobilab/falcon]](https://github.com/cobilab/falcon) | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://www.mdpi.com/2073-4425/9/9/445) |
| &#x1F49A;&nbsp; Bowtie2 | [[http://bowtie-bio.sourceforge.net/bowtie2]](http://bowtie-bio.sourceforge.net/bowtie2) | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)]( https://www.nature.com/articles/nmeth.1923) |
| &#x1F49A;&nbsp; SPAdes | [[http://cab.spbu.ru/software/spades/]](http://cab.spbu.ru/software/spades/) | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://www.liebertpub.com/doi/full/10.1089/cmb.2012.0021) | 
| &#x1F49A;&nbsp; Samtools | [[http://samtools.sourceforge.net/]](http://samtools.sourceforge.net/) | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://academic.oup.com/bioinformatics/article/25/16/2078/204688) | 
| &#x1F49A;&nbsp; Bcftools | [[http://www.htslib.org/doc/bcftools.html]](http://www.htslib.org/doc/bcftools.html) | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://academic.oup.com/bioinformatics/article/27/21/2987/217423) |
| &#x1F49A;&nbsp; IGV | [[https://software.broadinstitute.org/software/igv/]](https://software.broadinstitute.org/software/igv/) | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](http://cancerres.aacrjournals.org/content/77/21/e31.long) |


## 7. Citation ##

If you use this pipeline, please cite:
```
[Manuscript in preparation]
```

## 8. Issues ##

For any issue let us know at [issues link](https://github.com/pratas/traces/issues).

## 9. License ##

GPL v3.

For more information see LICENSE file or visit
<pre>http://www.gnu.org/licenses/gpl-3.0.html</pre>

