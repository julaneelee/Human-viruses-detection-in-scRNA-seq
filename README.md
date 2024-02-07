# Uveitis-scRNA-seq

Human known viruses identification from scRNA-seq data using Bowtie2 and blastn

Dependencies
- Blast
```
wget  https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.15.0+-x64-linux.tar.gz
tar -xvzf ncbi-blast-2.15.0+-x64-linux.tar.gz
mv ncbi-blast-2.15.0+ blast
export PATH=/path/to/directory/blast/bin:$PATH
```

```
wget -O Table_human_viruses.txt https://viralzone.expasy.org/resources/Table_human_viruses.txt?
```

Reference: https://github.com/caleblareau/serratus-reactivation-screen/tree/main/serratus_data_setup
