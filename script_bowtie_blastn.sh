#!/bin/bash
cd /PATH/TO/DIRECTORY
viral_reference_genome_array=("cmv" "vzv" "ebv")   # Note: This is array of indexes of our interested viral reference genomes 

#In my study, there are 4 scRNA-seq files (.fastq.gz) in directory and all ended with suffix `2_001.fastq.gz` so it would be selected in each loop
for i in *2_001.fastq.gz;
do
	
	echo $i
	##Filtering reads and Removing background
	#1 Filtering by Fastp
	docker run --user $(id -u):$(id -g) --rm -v `pwd`:`pwd` -w `pwd` quay.io/biocontainers/fastp:0.23.4--hadf994f_2 fastp -i $i -o ${i%.fastq.gz}_trimmed.fastq.gz && 
	#2 Quality checking of trimmed reads using fastqc
 	docker run --rm --user $(id -u):$(id -g) -v `pwd`:`pwd` -w `pwd` quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0 fastqc ${i%.fastq.gz}_trimmed.fastq.gz &&
	#3 Mapping with human reference genome in order to remove background of samples
	docker run --rm --user $(id -u):$(id -g) -v `pwd`:`pwd` -w `pwd` quay.io/biocontainers/bowtie2:2.5.2--py39h6fed5c7_0 bowtie2 -x hg38 -U ${i%.fastq.gz}_trimmed.fastq.gz -S ${i%.fastq.gz}.sam &&
	#4 Converting SAM to BAM
	docker run --rm --user $(id -u):$(id -g) -v `pwd`:`pwd` -w `pwd` quay.io/biocontainers/samtools:1.19--h50ea8bc_0 samtools view -bS ${i%.fastq.gz}.sam -o ${i%.fastq.gz}.bam &&
	#5 Sorting BAM 
	docker run --rm --user $(id -u):$(id -g) -v `pwd`:`pwd` -w `pwd` quay.io/biocontainers/samtools:1.19--h50ea8bc_0 samtools sort ${i%.fastq.gz}.bam -o ${i%.fastq.gz}_sorted.bam &&
	#6 Generating BAM index
	docker run --rm --user $(id -u):$(id -g) -v `pwd`:`pwd` -w `pwd` quay.io/biocontainers/samtools:1.19--h50ea8bc_0 samtools index ${i%.fastq.gz}_removeBG_sorted.bam &&
	#7 Extracting only unmapped reads 
	docker run --rm --user $(id -u):$(id -g) -v `pwd`:`pwd` -w `pwd` quay.io/biocontainers/samtools:1.19--h50ea8bc_0 samtools view -b -f 4 ${i%.fastq.gz}_removeBG_sorted.bam > ${i%.fastq.gz}_unmapped_reads.bam &&
	#8 Generating BAM index of unmapped reads 
	docker run --rm --user $(id -u):$(id -g) -v `pwd`:`pwd` -w `pwd` quay.io/biocontainers/samtools:1.19--h50ea8bc_0 samtools index ${i%.fastq.gz}_unmapped_reads.bam &&
	#9 Converting Bam to Fastq
	docker run --rm --user $(id -u):$(id -g) -v `pwd`:`pwd` -w `pwd` quay.io/biocontainers/samtools:1.19--h50ea8bc_0 samtools fastq -f 4 ${i%.fastq.gz}_unmapped_reads.bam > ${i%.fastq.gz}_unmapped_reads.fastq &&
	#10 Generate Fastq file for blastn
	docker run --rm --user $(id -u):$(id -g) -v `pwd`:`pwd` -w `pwd` quay.io/biocontainers/seqtk:1.4--he4a0461_1 seqtk seq -A ${i%.fastq.gz}_unmapped_reads.fastq > ${i%.fastq.gz}_unmapped_reads.fasta
	echo 'Removing background finished'
done

### FOR BLAST NCBI ### GENERATE BLASTDB
#docker run --rm -v `pwd`:`pwd` -w `pwd` quay.io/biocontainers/seqtk:1.4--he4a0461_1 seqtk seq -A herpesvirales.fna.gz > herpesvirales.fna
#export PATH=/mnt/data/julanee/sc_rna/blast/bin:$PATH
#makeblastdb -in herpesvirales.fna.gz -dbtype nucl -out herpesvirales

for j in *unmapped_reads.fastq;
do
	echo $j

	#Bowtie-2 mapping
	for m in "${viral_reference_genome_array[@]}";
	do 
		echo $m
		#0. Generate index for reference (already generated) 
		#docker run --rm --user 1012:1012 -v `pwd`:`pwd` -w `pwd` quay.io/biocontainers/bowtie2:2.5.2--py39h6fed5c7_0 bowtie2-build -f ${m}.fasta $m &&
		#1. Aligning with viral reference 
		docker run --rm --user $(id -u):$(id -g) -v `pwd`:`pwd` -w `pwd` quay.io/biocontainers/bowtie2:2.5.2--py39h6fed5c7_0 bowtie2 -x $m -U $j -S ${j%.fastq}_${m}.sam &&
		#2. Converting .sam to .bam
		docker run --rm --user $(id -u):$(id -g) -v `pwd`:`pwd` -w `pwd` quay.io/biocontainers/samtools:1.19--h50ea8bc_0 samtools view -bS ${j%.fastq}_${m}.sam -o ${j%.fastq}_${m}.bam &&
		#3. Sort Bam
		docker run --rm --user $(id -u):$(id -g) -v `pwd`:`pwd` -w `pwd` quay.io/biocontainers/samtools:1.19--h50ea8bc_0 samtools sort ${j%.fastq}_${m}.bam -o ${j%.fastq}_sorted_${m}.bam &&
		#4. Generating BAM index
		docker run --rm --user $(id -u):$(id -g) -v `pwd`:`pwd` -w `pwd` quay.io/biocontainers/samtools:1.19--h50ea8bc_0 samtools index ${j%.fastq}_sorted_${m}.bam &&
		#5. Read depth checking
		docker run --rm --user $(id -u):$(id -g) -v `pwd`:`pwd` -w `pwd` quay.io/biocontainers/samtools:1.19--h50ea8bc_0 samtools depth ${j%.fastq}_sorted_${m}.bam > ${j%.fastq}_depth_${m}.txt &&
		echo 'Mapping with Bowtie2 finished'
	done
	
	#Blastn
	export PATH=/path/to/blast/bin:$PATH
	makeblastdb -in herpesvirales.fna -dbtype nucl -out herpesvirales
	blastn -query ${j%.fastq}.fasta -db herpesvirales -out  ${j%.fastq}_herpesvirales.txt -perc_identity 80 -outfmt "6 qseqid stitle pident qseq length mismatch qstart qend sstart send evalue"
	blastn -query ${j%.fastq}.fasta -db 3ref1ct -out  ${j%.fastq}_3ref1ct.txt -perc_identity 80 -outfmt "6 qseqid stitle pident qseq length mismatch qstart qend sstart send evalue" 
done
