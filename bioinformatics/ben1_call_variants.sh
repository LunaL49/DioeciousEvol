#!/bin/bash

# Example calling syntax: bash ben1_call_variants.sh  R2
# Input argument: the pair of reverse demultiplexed file to be used, e.g "R2", "R6"
# Output files: 6 .vcf files containing mpileup result at the ben-1 locus of interest

# Process input
R="$1"

# Step 3: Forward demultiplex - Cutadapt
cutadapt --revcomp -j 4 -e 1 --no-indels \
  -g ^file:ben1F.fasta \
  -o {name}.R1.fastq.gz \
  -p {name}.R2.fastq.gz \
  ben-1${R}.R1.fastq.gz ben-1${R}.R2.fastq.gz

# Step 4: Process forward demultiplexed reads to call variant at locus of interest
for file in ben1F*.R1.fastq.gz
do
	name=$(basename "$file" .R1.fastq.gz)
	R2=$(echo "$file" | sed 's/R1/R2/')
	echo "$R2"
	bwa mem -t 4 ben-1.fasta $file $R2 | samtools view -b - | samtools sort -o ${name}.bam 
	samtools index ${name}.bam
	bcftools mpileup -d 1000000 -r ben-1:48 -f ben-1.fasta ${name}.bam -o ${name}.vcf
done