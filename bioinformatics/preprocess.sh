#!/bin/bash

# Calling syntax: bash pipeline.sh [name of raw R1 file] [name of raw R2 file]
# Input files: 2 raw paired-end sequencing output files 
# Output files: 2*10 reverse demultiplexed fastq.gz files in pairs

# Process inputs
R1="$1"
R2="$2"

# Step 1: QC - Trimmomatic (paired-end)
trimmomatic PE \
  "$R1" "$R2" \
  R1_paired_QCed.fastq.gz R1_unpaired.fastq.gz \
  R2_paired_QCed.fastq.gz R2_unpaired.fastq.gz \
  SLIDINGWINDOW:4:25 MINLEN:50

# Step 2: Reverse demultiplex - Cutadapt
cutadapt --revcomp -j 4 -e 1 --no-indels \
  -g ^file:barcodes_R.fasta \
  -o "${NAME}.R1.fastq.gz" \
  -p "${NAME}.R2.fastq.gz" \
  R1_paired_QCed.fastq.gz R2_paired_QCed.fastq.gz

