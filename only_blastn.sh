#!/bin/bash

#Activate conda environment (r_env environment)
#run Rscript step2_gen_efetch_fasta.R
#Obtain step3_run_efetch.sh (which contains efetch commands to download fasta of all human viruses)
#run ./step3_run_efetch.sh
#After finished downloading all human viruses fasta file then ... remove .fasta file that has no data (due to step2 trying to write efetch command t oretrieve all possible version of complete genome in NCBI NC_XXXXX.01-.09) 
#After removing files that have not data then concatenate all human viruses files  (1700+ file fasta) into single file (all_human_viruses.fasta)
#Make database for blast using makeblastdb of file name all_human_viruses.fasta
#blastn each scRNA file with the database 
#Obtain .txt file that contain human virus matched with our scRNA data 

cd /mnt/data/julanee/sc_rna
#Blastn
export PATH=/mnt/data/julanee/sc_rna/blast/bin:$PATH
#Merge(cat) all human viruses .fasta file together before using makeblastdb (ex. herpesvirales.fna) 
human_viruses="all_human_viruses"

#Find all file in directory that has no second line then remove !
find . -type f -exec awk -v x=2 'NR==x{exit 1}' {} \; -exec rm -f {} \;

cat NC*.fasta > all_human_viruses.fasta
mv NC*.fasta nc_pulls
cp all_human_viruses.fasta nc_pulls


#BLASTn
makeblastdb -in ${human_viruses}.fasta -dbtype nucl -out $human_viruses

for j in *unmapped_reads_default.fasta; do
	blastn -query $j -db $human_viruses -out  ${j%.fasta}_${human_viruses}.txt -perc_identity 80 -outfmt "6 qseqid stitle pident qseq length mismatch qstart qend sstart send evalue" -num_threads 8
done

#list of unmapped_reads_default.fasta (4 files)
#GEX_AAU1_S3_L006_R2_001_unmapped_reads_default.fasta
#GEX_AAU1_S3_L007_R2_001_unmapped_reads_default.fasta
#GEX_AAU1_S3_L008_R2_001_unmapped_reads_default.fasta
#GEX_AAU2_S1_L006_R2_001_unmapped_reads_default.fasta
