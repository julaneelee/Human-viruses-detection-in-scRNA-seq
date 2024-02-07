# Uveitis-scRNA-seq

Human known viruses identification from scRNA-seq data using Bowtie2 and blastn

Dependencies
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
Reference: https://github.com/caleblareau/serratus-reactivation-screen/tree/main/serratus_data_setup
