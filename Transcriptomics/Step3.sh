#!/bin/bash

#Do de novo transcriptome assembly for clean reads from chromosome 4 using Trinity
cd ~/Assignment6/results/3_denovo_assembly
Trinity --seqType fq --left ~/Assignment6/results/2_clean_data/Col_leaf_chr4_R1.clean.fastq.gz \
--right ~/Assignment6/results/2_clean_data/Col_leaf_chr4_R2.clean.fastq.gz \
--output Col_leaf_chr4_trinity --CPU 2 --max_memory 8G --bypass_java_version_check

#Get the genome coordinates file (GFF3) of de novo assembled transcripts using GMAP
cd ~/Assignment6/results/3_denovo_assembly
gmap -D ~/Assignment6/DB/TAIR10_GMAP -d TAIR10_GMAP -t 2 -f 3 -n 1 ./Col_leaf_chr4_trinity/Trinity.fasta > Trinity.gff3

#Assessing the assembly quality using BUSCO
cd ~/Assignment6/results/3_denovo_assembly
busco -i ./Col_leaf_chr4_trinity/Trinity.fasta -l ~/Assignment6/DB/brassicales_odb10 -o BUSCO_Trinity_brassicales -m transcriptome --cpu 2


