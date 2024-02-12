# Human viral detection in scRNA-seq data 
## Human known viruses identification from scRNA-seq data using Bowtie2 and blastn
There's a different purpose of using these tools, Blastn and Bowtie2 

|             |    Blastn     |     Bowtie2   |
|-------------| ------------- | ------------- |
| Purpose     | Searching our scRNA-seq reads in known human viruses database to find viral nucleotide sequence similarity | Mapping our scRNA-seq reads to specific viral reference genomes |

## Workflow of our study

```mermaid
graph TD;
    A[scRNA-seq]-->|Filter read quality by fast[| B[trimmed.fast.gz];
    B-->|Quality check by fastqc| C[Pass];
    C--> |Map with human reference genome using Bowtie2| E[.sam];
    B-->|Quality check by fastqc| D[No];
    D-->|Filter read quality by fast[| B[trimmed.fast.gz];
    E-->|Convert to binary file| F[.bam];
    F-->|Generate BAM index| G[BAM indexes];
```

## Dependencies
- Docker
- Anaconda
- Entrez Direct
```
sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
echo "export PATH=\$HOME/edirect:\$PATH" >> $HOME/.bash_profile
```
- Blast
```
wget  https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.15.0+-x64-linux.tar.gz
tar -xvzf ncbi-blast-2.15.0+-x64-linux.tar.gz
mv ncbi-blast-2.15.0+ blast
export PATH=/path/to/directory/blast/bin:$PATH
```
- fastp
```
docker pull quay.io/biocontainers/fastp:0.23.4--hadf994f_2
```
- fastqc
```
docker pull quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0
```
- bowtie2
```
docker pull quay.io/biocontainers/bowtie2:2.5.3--py310ha0a81b8_0
```
- samtools
```
docker pull quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0
```
- seqtk
```
docker pull quay.io/biocontainers/seqtk:1.4--he4a0461_1
```

## [[For Bowtie2]] Retrieving reference genome and generating index files 
Before mapping our scRNA-seq data with reference sequences, we need to prepare index files of reference genomes.
We need to prepare index files of `1) Human reference genome`  remove human genome background from our scRNA seq data and `2) Viral reference genome` we interested, in order to identify whether interested viral sequences are in our scRNA-seq data or not

### Downloading human reference genome from NCBI
```
curl https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001405.40/download?include_annotation_type=GENOME_FASTA,GENOME_GFF --output hg38.zip
unzip hg38.zip
cd /ncbi_dataset/data/GCF_000001405.40 
mv GCF_000001405.40_GRCh38.p14_genomic.fna hg38.fna
```
### or download human reference genome by using `efetch ` command of Entrez ID tool *We need to know accession number of the reference genome*
```
efetch -db nuccore -id NC_000001.11 -format fasta > hg38.fasta
```

### Downloading viral reference genome from NCBI using `efetch ` command of Entrez ID tool , *Note:We need to know accession number of the reference genome*
e.g.
```
#Cytomegalovirus complete genome (Human betaherpesvirus 5)
efetch -db nuccore -id NC_006273.2 -format fasta > cmv.fasta

# Epstein-barr virus complete genome (Human gammaherpesvirus 4)
efetch -db nuccore -id NC_007605.1 -format fasta > ebv.fasta

#Vericell-zoster virus complete genome (Human alphaherpesvirus 3)
efetch -db nuccore -id NC_001348.1 -format fasta > vzv.fasta
```

### After we downloaded reference genomes (contains in same directory), we will generate index files in order to map with our scRNA-seq with Bowtie2 `bowtie2-build` command
```
#Generating index file of Human reference genome
docker run --rm -v `pwd`:`pwd` -w `pwd` quay.io/biocontainers/bowtie2:2.5.2--py39h6fed5c7_0 bowtie2-build -f hg38.fasta hg38

#Generating index file of interested viral reference genome (e.g. we interested in Cytomegalovirus, Epstein-barr virus, Vericello-zoster virus)
docker run --rm -v `pwd`:`pwd` -w `pwd` quay.io/biocontainers/bowtie2:2.5.2--py39h6fed5c7_0 bowtie2-build -f cmv.fasta cmv       # → Cytomegalovirus(CMV)
docker run --rm -v `pwd`:`pwd` -w `pwd` quay.io/biocontainers/bowtie2:2.5.2--py39h6fed5c7_0 bowtie2-build -f ebv.fasta ebv       # → Epstein-barr virus (EBV)
docker run --rm -v `pwd`:`pwd` -w `pwd` quay.io/biocontainers/bowtie2:2.5.2--py39h6fed5c7_0 bowtie2-build -f vzv.fasta vzv       # → Vericello-zoster virus (VZV)
```


## 1. Download all known human viruses data from `ViralZone` 
```
wget -O Table_human_viruses.txt https://viralzone.expasy.org/resources/Table_human_viruses.txt?
```
## 2. Generate list of `efetch` commands into bash script
```
conda create -n r_env r-essentials r-base
conda activate r_env
Rscript create_fetch.R
```
## 3. Run bash script to fetch all known human viruses .fasta files 
```
./only_blastn.sh
```
## 4. Obtain file of matched human viruses `exp. human_viruses_detected_scRNA.txt`

Reference: https://github.com/caleblareau/serratus-reactivation-screen/tree/main/serratus_data_setup
