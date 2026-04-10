<br>
<div align="center">
  
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)
[![Speed](https://img.shields.io/static/v1.svg?label=Speed&message=Ultra-fast&color=green)](#)
<!-- BEGIN AUTO-GENERATED SECTION - VersionShield : DO NOT MODIFY-->
<!-- ================================================== -->
<!-- This snippet generated on Fri Apr 10 16:44:56 EEST 2026, -->
<!-- using src/TRACES_generator_scripts/gen_readme.sh -->
<!-- based on information in system_files/dependencies.yml -->
<!-- ================================================== -->
[![Version](https://img.shields.io/static/v1.svg?label=Release&message=v1.14.4&color=orange)](#)
<!-- ================================================== -->
<!-- END AUTO-GENERATED SECTION : DO NOT MODIFY-->
[![TinyURL](https://img.shields.io/static/v1.svg?label=TinyURL&message=traces-pipe&color=yellow)](https://tinyurl.com/traces-pipe)
[![DOI](https://zenodo.org/badge/175676601.svg)](https://zenodo.org/badge/latestdoi/175676601)

</div>
<br>
<p align="center">
<img src="imgs/logo.png" alt="TRACES Pipeline" height="200" border="0" /><br>
<i>A hybrid pipeline for reconstruction & analysis<br> 
of viral and host genomes at multi-organ level.</i>
<br><br>

## 1. About ##

TRACESPipe is a next-generation sequencing pipeline for identification, assembly, and analysis of viral and human-host genomes at multi-organ level. The identification and assembly of viral genomes rely on cooperation between three modalities:
<ul>
<li>compression-based predictors;</li>
<li>sequence alignments;</li>
<li><i>de-novo</i> assembly.</li>
</ul>
The compression-based prediction applies FALCON-meta technology with ultra-fast comparative quantification to find the best reference genome (from a large viral database) containing the highest similarity relative to the sequenced reads. After identification, the reads are aligned according to the best reference by Bowtie2. A consensus sequence is produced with specific filters using Bcftools. Then, <i>de-novo</i> assembly (metaSPAdes) is involved in building scaffolds. The high coverage scaffolds that overlap totally or partially the consensus sequence (aligned by bwa) are used to validate or either augment the new genome. The final analysis of the assembly is interactively supervised with the IGV with the goal of drafting the final sequence.

For the human-host variant call identification, the same procedure is followed although directly starting within the second point, given the use of the same reference (revised Cambridge Reference) to all the cases.

<br>
<p align="center">
<img src="imgs/pipeline.png" alt="TRACESPipe architecture" height="500" border="0" />
</p>
<br>

The previous image shows the architecture of TRACESPipe, where the green line stands for the mitochondrial human line. This pipeline has been tested in Illumina HiSeq and NovaSeq platforms. The operating system required to run it is Linux. In windows use cygwin (https://www.cygwin.com/) and make sure that it is included in the installation: cmake, make, zcat, unzip, wget, tr, grep (and any dependencies). If you install the complete cygwin packet then all these will be installed. After, all steps will be the same as in Linux.

The TRACESPipe includes methods for ancient DNA authentication, namely using the quantification of damage (in the tips of the reads) relative to a reference. Other feature is the quantification of y-chromosome presence through compression-based predictors.

Additionally, the TRACESPipe includes read trimming and filtering, PhiX removal, and redundancy controls (at the Database level and for each candidate reference genomes) to improve the consistency and quality of the data.

## 2. Installation, Structure and Configuration ##

### 2.1 Installation ###

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
To install TRACESPipe, run the following commands in a Linux OS:
```
git clone https://github.com/viromelab/tracespipe.git
cd tracespipe/src/
chmod +x TRACES*.sh
./TRACESPipe.sh --install
./TRACESPipe.sh --get-all-aux
```

#### Development note 

Install, Update, Version, and Check scripts, as well as this README have sections which are
automatically generated based on the dependencies using:

```
make all
```

As user you should not need to run this command, but if you have yq (a CLI YAML parser), you may.
As a developer this should be run whenever there are changes to the `system_files/dependencies.yml` file,
the generator scripts, or the relevant files. A suggestion is to add this to a `pre_commit` git hook.


### 2.2 Structure ###

In the tracespipe/ folder the following structure exists:
```bash
tracespipe/
│   
├── meta_data/         # information about the filenames in input_data/ and organ names
│   └── meta_info.txt  # see Configuration section for this file.
│   
├── input_data/        # where the NGS reads must be placed (and compressed with gzip)
│   
├── output_data/       # where the results will appear using the following subfolders: 
│   │
│   │
│   ├── TRACES_preprocessed_reads/     # trimmed and adapter removed fastq files
│   ├── TRACES_results/                # where the files regarding the metagenomic 
│   │                                  # analysis, redundancy (complexity) and control will appear
│   ├── TRACES_results/profiles/       # where the redundancy (complexity) profiles appear 
│   │
│   ├── TRACES_viral_alignments/       # where viral alignments and index will appear
│   ├── TRACES_viral_consensus/        # where viral consensus (FASTA) will appear
│   ├── TRACES_viral_bed/              # where viral BED files will appear (SNPs and Coverage)
│   ├── TRACES_viral_statistics/       # where viral statistics appear (depth/wide coverage)
│   │
│   ├── TRACES_mtdna_alignments/       # where mtdna alignments and index will appear
│   ├── TRACES_mtdna_consensus/        # where mtdna consensus (FASTA) will appear
│   ├── TRACES_mtdna_bed/              # where mtdna BED files will appear (SNPs and Coverage)
│   ├── TRACES_mtdna_statistics/       # where mtdna statistics appear (depth/wide coverage)
│   ├── TRACES_mtdna_authentication/   # where mtdna species and population authentication appears
│   │
│   ├── TRACES_cy_alignments/          # where cy alignments and index will appear
│   ├── TRACES_cy_consensus/           # where cy consensus (FASTA) will appear
│   ├── TRACES_cy_bed/                 # where cy BED files will appear (SNPs and Coverage)
│   ├── TRACES_cy_statistics/          # where cy statistics appear (depth/wide coverage)
│   │
│   ├── TRACES_specific_alignments/    # where specific alignments and index will appear
│   ├── TRACES_specific_consensus/     # where specific consensus (FASTA) will appear
│   ├── TRACES_specific_bed/           # where specific BED files will appear
│   ├── TRACES_specific_statistics/    # where specific statistics appear (depth/wide coverage)
│   │
│   ├── TRACES_mtdna_damage_<ORGAN>/   # where the mtdna damage estimation files will appear
│   │
│   ├── TRACES_denovo_<ORGAN>/         # where the output of de-novo assembly appears
│   │
│   ├── TRACES_hybrid_alignments/      # where the hybrid data appears
│   ├── TRACES_hybrid_consensus/       # where the hybrid data appears
│   ├── TRACES_hybrid_bed/             # where the hybrid data appears
│   │
│   ├── TRACES_hybrid_R2_alignments/   # where the second round hybrid data appears
│   ├── TRACES_hybrid_R2_consensus/    # where the second round hybrid data appears
│   ├── TRACES_hybrid_R2_bed/          # where the second round hybrid data appears
│   │
│   ├── TRACES_hybrid_R3_alignments/   # where the third round hybrid data appears
│   ├── TRACES_hybrid_R3_consensus/    # where the third round hybrid data appears
│   ├── TRACES_hybrid_R3_bed/          # where the third round hybrid data appears
│   │
│   ├── TRACES_hybrid_R4_alignments/   # where the fourth round hybrid data appears
│   ├── TRACES_hybrid_R4_consensus/    # where the fourth round hybrid data appears
│   ├── TRACES_hybrid_R4_bed/          # where the fourth round hybrid data appears
│   │
│   ├── TRACES_hybrid_R5_consensus/    # where the automatic choosen hybrid consensus 
│   │                                  # appears (diff will be made using this data)
│   │
│   ├── TRACES_multiorgan_alignments/  # where the multi-organ alignments data appears
│   ├── TRACES_multiorgan_consensus/   # where the multi-organ consensus data appears
│   │
│   ├── TRACES_diff/                   # where the dnadiff results appear (identity & SNPs)
│   ├── TRACES_specific_diff/          # where the dnadiff results appear for specific
│   │
│   └── TRACES_blasts/                 # where the specific blasted results appears
│   
├── to_encrypt_data/    # where the NGS files to encrypt must be before encryption
├── encrypted_data/     # where the encrypted data will appear
├── decrypted_data/     # where the decrypted data will appear
│   
├── logs/               # where the logs (stdout, stderr, and system) will appear
│   
├── src/                # where the bash code is and where the commands must be call
│   
└── imgs/               # images related with the pipeline
```

### 2.3 Configuration ###

To configure TRACESPipe add your <b>FASTQ files gziped</b> at the folder
```
input_data/
```
Then, add a file exclusively with name <b>meta\_info.txt</b> at the folder
```
meta_data/
```
This file needs to specify the organ type (with a single word name) and the filenames for the paired end reads. An example of the content of meta\_info.txt is the following:
```
skin:V1_S44_R1_001.fastq.gz:V1_S44_R2_001.fastq.gz
brain:V2_S29_R1_001.fastq.gz:V2_S29_R2_001.fastq.gz
colon:V3_S45_R1_001.fastq.gz:V3_S45_R2_001.fastq.gz
```
Then, at the <b>src/</b> folder run:
```
./TRACESPipe.sh --get-all-aux
```

## 3. Running ##

To run TRACES Pipeline, use the following command:
```
./TRACESPipe.sh <parameters>
```
There are many parameters and configurations that can be used.<br>
See the next section for more information about the usage.

## 4. Usage ##

```
./TRACESPipe.sh -h
```

<!-- BEGIN AUTO-GENERATED SECTION - TracesHelp : DO NOT MODIFY -->
<!-- ================================================== -->
<!-- This snippet generated on Fri Apr 10 16:44:56 EEST 2026, -->
<!-- using src/TRACES_generator_scripts/gen_readme.sh -->
<!-- based on information in system_files/dependencies.yml -->
<!-- ================================================== -->
 ```
                                                                       
                                                              
         ████████╗ ██████╗   █████╗   ██████╗ ███████╗ ███████╗        
         ╚══██╔══╝ ██╔══██╗ ██╔══██╗ ██╔════╝ ██╔════╝ ██╔════╝        
            ██║    ██████╔╝ ███████║ ██║      █████╗   ███████╗        
            ██║    ██╔══██╗ ██╔══██║ ██║      ██╔══╝   ╚════██║        
            ██║    ██║  ██║ ██║  ██║ ╚██████╗ ███████╗ ███████║        
            ╚═╝    ╚═╝  ╚═╝ ╚═╝  ╚═╝  ╚═════╝ ╚══════╝ ╚══════╝        
                                                                       
                             P I P E L I N E                           
                                                                       
          |  A hybrid pipeline for reconstruction & analysis  | 
          |  of viral and host genomes at multi-organ level.  | 
                                                                
    Usage: ./TRACESPipe.sh [options]                     
                                                                       
 ===========            GENERAL OPTIONS              ==========
                                                                       
    -h,     --help            Show this help message and exit,         
    -v,     --version         Show the version and some information,   
    -flog,  --flush-logs      Flush logs (delete logs),                
    -fout,  --flush-output    Flush output data (delete all output_data), 
    -t <THREADS>, --threads <THREADS> Number of threads to use,        
                                                                       
 ===========            SETUP COMMANDS               ==========
                                                                       
    -i,     --install         Installation of all the tools,           
    -up,    --update          Update all the tools in TRACESPipe,      
    -spv,   --show-prog-ver   Show included programs versions,         
                                                                       
    -st,    --sample          Creates human ref. VDB and sample organ, 
                                                                       
    -gmt,   --get-max-threads Get the number of maximum machine threads,
                                                                       
    -dec,   --decrypt         Decrypt (all files in ../encrypted_data), 
    -enc,   --encrypt         Encrypt (all files in ../to_encrypt_data),
                                                                       
    -vdb,   --build-viral     Build viral database (all) [Recommended], 
    -vdbr,  --build-viral-r   Build viral database (references only),  
    -udb,   --build-unviral   Build non viral database (control),      
    -lcm,   --lcr-mask-vdb    Construct an LCR masked viral database   
                              Uses alt-vdb if specified                
                                                                       
    -afs <FASTA>, --add-fasta <FASTA>                                  
                              Add a FASTA sequence to the viral database,
    -aes <ID>, --add-extra-seq <ID>                                    
                              Add extra sequence to the viral database,
    -gx,    --get-extra-vir   Downloads/appends (VDB) extra viral seq, 
                                                                       
    -gad,   --gen-adapters    Generate FASTA file with adapters,       
    -gp,    --get-phix        Extracts PhiX genomes (Needs viral DB),  
    -gm,    --get-mito        Downloads human Mitochondrial genome,    
                                                                       
    -dwms,  --download-mito-species                                    
                              Downloads the complete NCBI mitogenomes  
                              database containing the existing species,
                                                                       
    -dwmp,  --download-mito-population                                 
                              Downloads two complete mitogenome databases
                              with healthy and pathogenic sequences,   
                                                                       
    -aums,  --auth-mito-species                                        
                              Autheticate the mitogenome species,      
                                                                       
    -aump,  --auth-mito-population                                     
                              Authenticate closest population,         
                                                                       
    -cmt <ID>, --change-mito <ID>                                      
                              Set any Mitochondrial genome by ID,      
                                                                       
    -gy,    --get-y-chromo    Downloads human Y-chromosome,            
    -gax,   --get-all-aux     Runs -gad -gp -gm -gy,                   
                                                                       
    -cbn,   --create-blast-db It creates a nucleotide blast database,  
    -ubn,   --update-blast-db It updates a nucleotide blast database,  
                                                                       
 ===========           ANALYSIS COMMANDS             ==========
                                                                       
    --- Some commands can only be run with or after others.            
        At the bottom of this section is a dependency tree             
                                                                       
    -ra,    --run-analysis    Run data analysis (core),                      
    -all,   --run-all         Run all the options (excluding the specific).  
    -proc,  --run-preprocess  Run adapter removal, quality trimming,   
                              length filtering, base correction,       
                              and poly-g tail removal with fastp       
                                                                       
    -sfs <FASTA>, --search-blast-db <FASTA>                            
                              It blasts the nucleotide (nt) blast DB,  
    -sfrs <FASTA>, --search-blast-remote-db <FASTA>                    
                              It blasts remotly thenucleotide (nt) blast 
                              database (it requires internet connection), 
                                                                       
    -gbb,   --best-of-bests   Identifies the best of bests references  
                              between multiple organs [similar reference], 
                                                                       
    -rm,    --run-meta        Run viral metagenomic identification,    
    -ro,    --run-meta-nv     Run NON-viral metagenomic identification,
    -rpro,  --run-profiles    Run complexity and relative profiles (control), 
                                                                       
    -rpgi <ID>,  --run-gid-complexity-profile <ID>                     
                              Run complexity profiles by GID,          
                                                                       
    -rava,  --run-all-v-alig  Run all viral align/sort/consensus seqs  
                              from a specific list,                    
                                                                       
    -rsr <ID>, --run-specific <ID/PATTERN>                             
                              Run specific reference align/consensus,  
                                                                       
    -rsx <ID>, --run-extreme <ID/PATTERN>                              
                              Run specific reference align/consensus   
                              using extreme sensitivity;               
                              Retained for backwards compatibility;    
                              Now an alias for -vhs -rsr <ID/PATTERN>, 
                                                                       
    -rmt,   --run-mito        Run Mito align and consensus seq,        
    -rmtd,  --run-mito-dam    Run Mito damage only,                    
                                                                       
    -rgid <ID>, --run-gid-damage <ID>                                  
                              Run damage pattern analysis by GID,      
                                                                       
    -rya,   --run-cy-align    Run CY align and consensus seq,          
    -ryq,   --run-cy-quant    Estimate the quantity of CY DNA,         
                                                                       
    -rda,   --run-de-novo     Run de-novo assembly,                    
                                                                       
    -rhyb,  --run-hybrid      Run hybrid assembly (align/de-novo),     
    -rsd <ID>, --run-de-novo-specific <ID/PATTERN>                     
                              Run specific alignments of the de-novo   
                              to the reference genome,                 
                                                                       
    -cast,  --compile-aln-stats <TYPE>                                 
                              Combine breadth, depth, similarity,      
                              and selected mapping statistics into a   
                              report. Valid types are: viral, mtdna    
                              cy, specific, and all. Auto-enabled      
                              with any alignment command.              
    -rmhc,  --run-multiorgan-consensus                                 
                              Run alignments/consensus between all the 
                              reconstructed organ sequences,           
                                                                       
    -vis,   --visual-align    Run Visualization tool for alignments,   
    -covl,  --coverage-latex  Run coverage table in Latex format,      
    -covc,  --coverage-csv    Run coverage table in CSV format,        
                                                                       
    -covp <NAME>,  --coverage-profile <BED_NAME_FILE>                   
                              Run coverage profile for specific BED file, 
                                                                       
    -diff,  --run-diff        Run diff -> reference and hybrid (ident/SNPs), 
                                                                       
    -sdiff <V_NAME> <ID/PATTERN>, --run-specific-diff <V_NAME> <ID/PATTERN>  
                              Run specific diff of reconstructed to a virus  
                              pattern of ID. Example: -sdiff B19 AY386330.1, 
                                                                       
    -brec,  --blast-reconstructed                                      
                              Run local blast over reconstructed genomes, 
                                                                             
                                                                       
    --- Dependency Tree ---                                            
                                                                       
    run-preprocess                                                     
    ├───run-meta                                                       
    │   ├───run-profiles                                               
    │   └───run-all-v-alig                                             
    │       ├───compile-aln-stats viral                                
    │       ├───coverage-latex                                         
    │       ├───coverage-csv                                           
    │       └───run-diff                                               
    ├───run-specific                                                   
    │   └───compile-aln-stats specific                                 
    ├───run-mito                                                       
    │   └───compile-aln-stats mtdna                                    
    ├───run-cy-align                                                   
    │   └───compile-aln-stats cy                                       
    ├───run-cy-quant                                                   
    └───run-de-novo                                                    
        ├───run-hybrid                                                 
        │   ├───run-multiorgan-consensus                               
        │   ├───run-diff                                               
        │   ├───run-specific-diff                                      
        │   └───blast-reconstructed                                    
        └───run-de-novo-specific                                       
    search-blast-db                                                    
    search-blast-remote-db                                             
    best-of-bests                                                      
    run-gid-complexity-profile                                         
    run-mito-dam                                                       
    run-gid-damage                                                     
    visual-align                                                       
    coverage-profile                                                   
                                                                       
 ===========           ANALYSIS OPTIONS              ==========
                                                                       
    -avdb <FASTA>, --alt-viral-db <FASTA>                              
                              Specify a path to fasta file containing  
                               viral sequences                         
                              Sequence names must include the accession
                               as the first field, either whitespace or
                               underscore (_) delimited (NC_ handled)  
    -vdbm <PATH>, --viral-db-metadata <PATH>                           
                              Specify a path to a tab del file which   
                               has sequence GID/ACC in the first column
                               and a Name representing a virus label   
                               in the second                           
                              This changes the meta-analysis behaviour:
                               rather than using internal virus names  
                               and inclusion criteria, the relationships
                               defined in the provided file are used.  
                              This allows the user to define groupings 
                               of interest in the default, or user     
                               provided viral database.                
                                                                       
    -rdup,  --remove-dup      Remove duplications (e.g. PCR dup),      
    -vhs,   --very-sensitive  Aligns with very high sensitivity (slower),  
                                                                       
    -adr,   --attempt-denovo-restart                                   
                              TRACESPipe will attempt to run previous  
                              de-novo assemblies from their last       
                              checkpoint, useful if an external fault  
                              halted assembly                          
                              Warning: If resuming initiates but later 
                              fails, assembly will restart from scratch
    -dmef <VALUE>,  --denovo-mem-estimate-factor <VALUE>     Default:50
                              A value multiplied by the file size of   
                              compressed raw foward reads to estimate  
                              the spades memory usage. Controls trigger
                              for bbnorm digital normalization.        
                              A value of 0 disables bbnorm.            
    -mdm,   --max-denovo-mem                                Default:350
                              Maximum memory in GB to allocate for     
                              de novo assemly. A value of 0 removes the
                              limit and disables bbnorm.               
                                                                       
    -ulcm, --use-lcr-masked-vdb                                                                   
                              Use the pre-constructed lcr masked vdb   
                              for input to FALCON; disabled unless -lcm
                              has been run                             
    -pdep <FILE>, --pattern-depletion <FILE>                           
                              A path to a file containing regular      
                              expressions. Matches are filtered out    
                              of FALCON input. Primary use is to filter
                              patterns which exist in both a virus and 
                              background sequences (eq. telomeres)     
    -iss <SIZE>, --inter-sim-size <SIZE>                               
                              Inter-genome similarity top size (control), 
                                                                       
    -cpwi <VALUE>, --complexity-profile-window <VALUE>                 
                              Complexity profile window size,          
    -cple <VALUE>, --complexity-profile-level <VALUE>                  
                              Complexity profile compression level [1;10], 
                                                                       
    -mis <VALUE>, --min-similarity <VALUE>                             
                              Minimum similarity value to consider the 
                              sequence for alignment-consensus (filter), 
                                                                       
    -misl <VALUE>, --min-similarity-len <VALUE>                        
                              Minimum product of similarity value and  
                              best hit sequence length for             
                              alignment-consensus (filter),            
    -misv <PATH>, --min-similarity-virus <PATH>                        
                              Path to a tab sep file with tow columns  
                              containing the virus and min sim values  
                              Any values lower than --misl are ignored 
                                                                       
    -top <VALUE>, --view-top <VALUE>                                   
                              Display the top <VALUE> with the highest 
                              similarity (by descending order),        
    -amax <VALUE>, --max-alignments <VALUE>                            
                              The maximum number of alignments to      
                              report for any one read; 0 = no limit,   
                                                                       
    -c <VALUE>,   --cache <VALUE>                                      
                              Cache to be used by FALCON-meta,         
    -tsv <VALUE>,   --top-size-virus <VALUE>                           
                              Top size to be used by FALCON-meta when  
                              using TRACES_metagenomic_viral.sh;       
                              default:0 -> seq count in viral db       
    -ts <VALUE>,   --top-size <VALUE>                                  
                              Top size to be used by FALCON-meta when  
                              using TRACES_metagenomic.sh;             
                              default:0 -> seq count in  non-viral db  
                                                                       
    -cmax <MAX>,   --max-coverage <MAX_COVERAGE>                       
                              Maximum depth coverage (depth normalization), 
    -clog <VALUE>, --coverage-log-scale <VALUE>                        
                              Coverage profile logarithmic scale VALUE=Base, 
    -cwis <VALUE>, --coverage-window-size <VALUE>                      
                              Coverage window size for low-pass filter, 
    -cdro <VALUE>, --coverage-drop <VALUE>                             
                              Coverage drop size (sampling),           
    -covm <VALUE>, --coverage-min-x <VALUE>                             
                              Coverage minimum value for x-axis        
                                                                       
 ===========                EXAMPLES                 ==========
                                                                       
    Ex: ./TRACESPipe.sh --flush-output --flush-logs --run-mito --run-meta 
    --remove-dup --run-de-novo --run-hybrid --min-similarity 1 --run-diff 
    --very-sensitive --best-of-bests --run-multiorgan-consensus 
                                                                       
    Add the file meta_info.txt at ../meta_data/ folder. Example:       
    meta_info.txt -> 'organ:reads_forward.fa.gz:reads_reverse.fa.gz'   
    The reads must be GZIPed in the ../input_data/ folder.             
    The output results are at ../output_data/ folder.                  
                                                                       
    Contact: tracespipe@gmail.com                        
                                                                       
 ```
<!-- ================================================== -->
<!-- END AUTO-GENERATED SECTION : DO NOT MODIFY -->

## 5. Examples ##

The common use of TRACESPipe as command is:
```
./TRACESPipe.sh \
--flush-logs \
--run-preprocess \
--run-meta \
--inter-sim-size 2 \
--run-all-v-alig \
--run-mito \
--remove-dup \
--run-de-novo \
--run-hybrid \
--min-similarity 1.5 \
--view-top 5 \
--best-of-bests \
--very-sensitive \
--run-multiorgan-consensus \
--run-diff
```
From the run all the output is provided at folder output\_data and it can be human inspected using IGV.

Nevertheless, for specific runs, below some examples are described.

### 5.1 Building viral consensus sequences with fixed reference sequence in all organs (if exists in the FASTQ samples): ###
```
./TRACESPipe.sh --run-meta --run-all-v-alig --remove-dup --min-similarity 3 --best-of-bests
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
./TRACESPipe.sh --run-mito --remove-dup
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

### 5.3 Encrypt and Decrypt NGS data: ###

TRACESPipe supports secure encryption of genomic data. This allows outsourcing of the sequencing service while maintaining secure transmission and storage of the files.

#### 5.3.1 Encrypt #### 

Place the files from sequencing (e.g. FASTQ gziped files) in the folder <b>to_encrypt_data</b> and, then, run:
```
./TRACESPipe.sh --encrypt
```
Insert a strong password.<br>
The encrypted files are in the <b>encrypted_data</b> folder.

#### 5.3.2 Decrypt #### 

Place the encrypted files in the folder <b>encrypted_data</b> and, then, run:
```
./TRACESPipe.sh --decrypt
```
Insert the password that has been used in encryption.<br>
The decrypted files are in the <b>decrypted_data</b> folder.

### 5.4 Run all viral genome alignments, variation, and consensus sequences: ###

```
./TRACESPipe.sh --run-meta --run-all-v-alig
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

### 5.5 Quantify the presence of y-chromosome: ###

```
./TRACESPipe.sh --run-cy-quant
```
The output quantify is included at
```
output_data/TRACES_results/REP_CY_<organ_name>.txt
```

### 5.6 Full viral metagenomic composition for all the organs: ###

```
./TRACESPipe.sh --run-meta
```
The output is included at
```
../output_data/TRACES_results/REPORT_META_VIRAL_ALL.txt
```

### 5.7 Run NON viral metagenomic composition for all the organs (fungi, archaea, etc): ###

```
./TRACESPipe.sh --run-meta-nv
```
The output is included at
```
../output_data/TRACES_results/REPORT_META_NON_VIRAL_<organ_name>.txt
```

### 5.8 Run de-novo assembly (all data): ###

```
./TRACESPipe.sh --run-de-novo
```
The outputs are included at
```
../output_data/TRACES_denovo_alignments
../output_data/TRACES_denovo_consensus
../output_data/TRACES_denovo_bed
```

### 5.9 Run specific viral alignment (AF037218.1) for all organs using extreme sensitivity without duplications: ###

```
./TRACESPipe.sh --remove-dup --run-extreme AF037218.1
```
The output is included at
```
../output_data/TRACES_specific_alignments
```
and the depth and breadth coverage values at
```
cat ../output_data/TRACES_specific_statistics
```

### 5.10 Evaluate damage of mitochondrial DNA ###

```
./TRACESPipe.sh --run-mito-dam
```
The output is included at
```
../output_data/TRACES_mtdna_damage_<organ_name>
```

### 5.11 Remote blastn search over nucleotide NCBI database ###

```
./TRACESPipe.sh --search-blast-remote-db AF037218.1
```
The output is included at 
```
../output_data/TRACES_blastn
```

### 5.12 Calculate depth coverage with normalized value ###

This approach assumes that the reconstruction has already been processed:
```
./TRACES_normalized_depth.sh ../output_data/TRACES_viral_bed/B19-coverage-blood.bed 200
```
The output is provided to the stdout.

## 6. Programs ##

TRACES Pipeline uses a combination of the following tools:

<!-- BEGIN AUTO-GENERATED SECTION - DepInfoTable : DO NOT MODIFY -->
<!-- ================================================== -->
<!-- This snippet generated on Fri Apr 10 16:44:57 EEST 2026, -->
<!-- using src/TRACES_generator_scripts/gen_readme.sh -->
<!-- based on information in system_files/dependencies.yml -->
<!-- ================================================== -->
| Tool | TestedVersion | Article |
| --- | --- | --- |
| &#x1F49A;&nbsp; [AdapterRemoval](https://adapterremoval.readthedocs.io) | 2.3.4 | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://doi.org/10.1186/s13104-016-1900-2) |
| &#x1F49A;&nbsp; [ART\_illumina](https://www.niehs.nih.gov/research/resources/software/biostatistics/art) | 2016.06.05 | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://doi.org/10.1093/bioinformatics/btr708) |
| &#x1F49A;&nbsp; [BBNorm](https://bbmap.org) | 39.79 | NA |
| &#x1F49A;&nbsp; [BFCtools](http://www.htslib.org/doc/bcftools.html) | 1.21 | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://doi.org/10.1093/gigascience/giab008) |
| &#x1F49A;&nbsp; [BEDOPS](https://bedops.readthedocs.io) | 2.4.42 | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://doi.org/10.1093/bioinformatics/bts277) |
| &#x1F49A;&nbsp; [BEDTools](https://bedtools.readthedocs.io/en/latest) | v2.31.1 | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://doi.org/10.1002/0471250953.bi1112s47) |
| &#x1F49A;&nbsp; [BLASTn](https://blast.ncbi.nlm.nih.gov) | 2.16.0+ | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://doi.org/10.1089/10665270050081478) |
| &#x1F49A;&nbsp; [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2) | 2.5.4 | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://doi.org/10.1038/nmeth.1923) |
| &#x1F49A;&nbsp; [BWA](http://bio-bwa.sourceforge.net) | 0.7.18-r1243-dirty | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://doi.org/10.1093/bioinformatics/btp324) |
| &#x1F49A;&nbsp; [Cryfa](https://github.com/cobilab/cryfa) | v18.06 | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://doi.org/10.1093/bioinformatics/bty645) |
| &#x1F49A;&nbsp; [dnadiff](https://mummer4.github.io) | 1.3 | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://doi.org/10.1371/journal.pcbi.1005944) |
| &#x1F49A;&nbsp; [efetch](https://www.ncbi.nlm.nih.gov/books/NBK179288) | 24.4 | NA |
| &#x1F49A;&nbsp; [FALCON](https://github.com/cobilab/falcon) | 2.3 | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://doi.org/10.3390/genes9090445) |
| &#x1F49A;&nbsp; [fastp](https://github.com/OpenGene/fastp) | 0.23.4 | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://doi.org/10.1002/imt2.70078) |
| &#x1F49A;&nbsp; [grepq](https://github.com/Rbfinch/grepq) | 1.5.4 | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://doi.org/10.21105/joss.08048) |
| &#x1F49A;&nbsp; [GTO](https://github.com/cobilab/gto) | v1.5.9 | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://doi.org/10.1101/2020.01.07.882845) |
| &#x1F49A;&nbsp; [IGV](https://software.broadinstitute.org/software/igv) | 2.19.3 | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://doi.org/10.1038/nbt.1754) |
| &#x1F49A;&nbsp; [iVar](https://andersen-lab.github.io/ivar) | 1.4.4 | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://doi.org/10.1186/s13059-018-1618-7) |
| &#x1F49A;&nbsp; [MAGNET](https://github.com/cobilab/magnet) | 19.4 | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://doi.org/10.23919/eusipco63174.2024.10715165) |
| &#x1F49A;&nbsp; [mapDamage2](https://ginolhac.github.io/mapDamage) | 2.2.3 | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://doi.org/10.1093/bioinformatics/btt193) |
| &#x1F49A;&nbsp; [SAMtools](https://www.htslib.org/) | 1.21 | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://doi.org/10.1093/gigascience/giab008) |
| &#x1F49A;&nbsp; [Sdust](https://github.com/lh3/sdust) | 0.1 | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://doi.org/10.1089/cmb.2006.13.1028) |
| &#x1F49A;&nbsp; [SPAdes](http://cab.spbu.ru/software/meta-spades) | 4.2.0 | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://doi.org/10.1089/cmb.2012.0021) |
| &#x1F49A;&nbsp; [Tabix](http://htslib.org/doc/tabix.html) | 1.21 | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://doi.org/10.1093/gigascience/giab007) |
| &#x1F49A;&nbsp; [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) | 0.39 | [![Article](https://img.shields.io/static/v1.svg?label=View&message=Article&color=green)](https://doi.org/10.1093/bioinformatics/btu170) |
<!-- ================================================== -->
<!-- END AUTO-GENERATED SECTION : DO NOT MODIFY -->

## 7. Citation ##

If you use this pipeline, please cite:
```
Pratas, D., Toppinen, M., Pyöriä, L., Hedman, K., Sajantila, A. and Perdomo, M.F., 2020. 
A hybrid pipeline for reconstruction and analysis of viral genomes at multi-organ level. 
GigaScience, 9(8), p.giaa086.
```
[PDF Link](https://doi.org/10.1093/gigascience/giaa086)

## 8. Issues ##

For any issue let us know at [issues link](https://github.com/viromelab/tracespipe/issues).

## 9. License ##

GPL v3.

For more information see LICENSE file or visit
<pre>http://www.gnu.org/licenses/gpl-3.0.html</pre>

