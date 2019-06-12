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
<b>A next-generation sequencing pipeline for identification, assembly, and analysis of viral and host genomes at multi-organ level<b/>.
<br>

## 1. Installation ##

To install TRACES Pipeline, run the following commands in a Linux OS:
```
git clone https://github.com/pratas/traces.git
cd traces/src/
chmod +x TRACE*.sh
./TRACES.sh --install
```
Conda is needed for installation. <br>
To install Conda use the following steps:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
Additional instructions can be found here:
```
https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html
```


## 2. Running ##

To run TRACES Pipeline, use the following commands in a Linux OS:
```
./TRACES.sh <parameters>
```

## 3. Usage ##

```
./TRACES.sh -h
```

```                                                     
                                                              
     ████████╗ ██████╗   █████╗   ██████╗ ███████╗ ███████╗         
     ╚══██╔══╝ ██╔══██╗ ██╔══██╗ ██╔════╝ ██╔════╝ ██╔════╝         
        ██║    ██████╔╝ ███████║ ██║      █████╗   ███████╗         
        ██║    ██╔══██╗ ██╔══██║ ██║      ██╔══╝   ╚════██║         
        ██║    ██║  ██║ ██║  ██║ ╚██████╗ ███████╗ ███████║         
        ╚═╝    ╚═╝  ╚═╝ ╚═╝  ╚═╝  ╚═════╝ ╚══════╝ ╚══════╝         
                                                                
    An automatic pipeline for viral and human genome analysis
    in the contexts of clinical virology and forensics.         
                                                                
    Usage: ./TRACES.sh [options]                                
                                                                
    -h,    --help            Show this help message and exit,     
                                                                  
    -i,    --install         Installation of all the tools,       
                                                                  
    -vdb,  --build-viral     Build viral database (all sequences), 
    -vdbr, --build-viral-r   Build viral database (references only),  
    -udb,  --build-unviral   Build non viral database (control),  
                                                                  
    -gx,   --get-extra-vir   Downloads/appends (VDB) extra viral seq, 
    -gad,  --gen-adapters    Generate FASTA file with adapters,   
    -gp,   --get-phix        Extracts PhiX genomes (Needs viral DB),  
    -gm,   --get-mito        Downloads human Mitochondrial genome,
    -gy,   --get-y-chromo    Downloads human Y-chromosome,        
                                                                  
    -ra,   --run-analysis    Run data analysis,                   
                                                                  
    -rm,   --run-meta        Run viral metagenomic identification,    
    -ro,   --run-meta-nv     Run NON-viral metagenomic identification,   
                                                                  
    -rava, --run-all-v-alig  Run all viral align/sort/consensus seqs,    
                                                                 
    -rb,   --run-b19         Run B19 align, sort and consensus seq,    
    -rh1,  --run-hv1         Run HV1 align, sort and consensus seq,    
    -rh2,  --run-hv2         Run HV2 align, sort and consensus seq,    
    -rh3,  --run-hv3         Run HV3 align, sort and consensus seq,    
    -rh4,  --run-hv4         Run HV4 align, sort and consensus seq,    
    -rh5,  --run-hv5         Run HV5 align, sort and consensus seq,    
    -rh6,  --run-hv6         Run HV6 align, sort and consensus seq,    
    -rh6a, --run-hv6a        Run HV6A align, sort and consensus seq,    
    -rh6b, --run-hv6b        Run HV6B align, sort and consensus seq,    
    -rh7,  --run-hv7         Run HV7 align, sort and consensus seq,    
    -rh8,  --run-hv8         Run HV8 align, sort and consensus seq,    
    -rtt,  --run-ttv         Run TTV align, sort and consensus seq,    
    -rjc,  --run-jcv         Run JCV align, sort and consensus seq,    
    -rmc,  --run-mcv         Run MCV align, sort and consensus seq,    
    -rbk,  --run-bk          Run BK align, sort and consensus seq,    
    -rbv1, --run-hbov1       Run HBoV1 align, sort and consensus seq,    
    -rbv0, --run-hbovnot1    Run HBoV 2,3,... align, sort and consensus seq,    
    -rhbv, --run-hbv         Run HBV align, sort and consensus seq,    
    -rhpv, --run-hpv         Run HPV align, sort and consensus seq,    
    -rvar, --run-varv        Run VARV align, sort and consensus seq,    
                                                                 
    -rsr,  --run-specific    Run specific REF align/consensus seq, 
                                                                 
    -rmt,  --run-mito        Run Mito align, sort and consensus seq,   
    -rcy,  --run-y-chromo    Run CY align, sort and consensus seq,    
                                                                  
    -rda,  --run-de-novo     Run de-novo assembly,               
                                                                 
    -all,  --run-all         Run all the options.                 
                                                                
    Example: ./TRACES.sh --run-meta --run-mito           
                                                                
    meta_info.txt -> 'name:reads_forward.fa.gz:reads_reverse.fa.gz'  
    The reads and meta_info.txt must be in the src/ folder.     
                                                                
```
## 4. Programs ##

TRACES Pipeline uses a combination of the following tools:
```
Cryfa
GTO
Trimmomatic
MAGNET
FALCON-meta
Bowtie2
SPAdes
Samtools
Bcftools
IGV
```

## 5. Citation ##

If you use this pipeline, please cite:
```
null
```

## 6. Issues ##

For any issue let us know at [issues link](https://github.com/pratas/traces/issues).

## 7. License ##

GPL v3.

For more information see LICENSE file or visit
<pre>http://www.gnu.org/licenses/gpl-3.0.html</pre>

