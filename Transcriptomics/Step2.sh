#!/bin/bash

#copy the sequencing data to data directory
cp ~/data/Transcriptomics_data/Assignment/Col_leaf_chr4_R*.fastq.gz ~/Assignment6/data/

#QC for raw reads
cd ~/Assignment6/results/1_QC/
fastqc -t 2 -o ./ ~/Assignment6/data/Col_leaf_chr4_R1.fastq.gz 
fastqc -t 2 -o ./ ~/Assignment6/data/Col_leaf_chr4_R2.fastq.gz

#Trim adaptor and low-quality sequences using cutadapt 
cd ~/Assignment6/results/2_clean_data/
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o Col_leaf_chr4_R1.clean.fastq.gz -p Col_leaf_chr4_R2.clean.fastq.gz --minimum-length 25 --quality-cutoff 20 ~/Assignment6/data/Col_leaf_chr4_R1.fastq.gz ~/Assignment6/data/Col_leaf_chr4_R2.fastq.gz 

#Get QC report for clean reads
#QC for clean reads
cd ~/Assignment6/results/2_clean_data/
fastqc -t 2 -o ./ Col_leaf_chr4_R1.clean.fastq.gz
fastqc -t 2 -o ./ Col_leaf_chr4_R2.clean.fastq.gz