#!/bin/bash

#Do genome guided transcriptome assembly for clean reads from chromosome 4 with STAR and StringTie
cd ~/Assignment6/results/4_genome_guided_assembly
STAR --genomeDir ~/Assignment6/DB/TAIR10_STAR125 --readFilesIn ~/Assignment6/results/2_clean_data/Col_leaf_chr4_R1.clean.fastq.gz \
~/Assignment6/results/2_clean_data/Col_leaf_chr4_R2.clean.fastq.gz --readFilesCommand zcat \
--runThreadN 2 --outSAMstrandField intronMotif --outSAMattributes All \
--outFilterMismatchNoverLmax 0.03 --alignIntronMax 10000 --outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix Col_leaf_chr4. --quantMode GeneCounts

#copy and paste this version
cd ~/Assignment6/results/4_genome_guided_assembly
STAR --genomeDir ~/Assignment6/DB/TAIR10_STAR125 --readFilesIn ~/Assignment6/results/2_clean_data/Col_leaf_chr4_R1.clean.fastq.gz ~/Assignment6/results/2_clean_data/Col_leaf_chr4_R2.clean.fastq.gz --readFilesCommand zcat --runThreadN 2 --outSAMstrandField intronMotif --outSAMattributes All --outFilterMismatchNoverLmax 0.03 --alignIntronMax 10000 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix Col_leaf_chr2. --quantMode GeneCounts

#Assembly using StringTie
cd ~/Assignment6/results/4_genome_guided_assembly
stringtie Col_leaf_chr4.Aligned.sortedByCoord.out.bam -o StringTie.gtf -p 2

#Get a fasta file of genome guided assembled transcripts using gffread
cd ~/Assignment6/results/4_genome_guided_assembly
gffread -w StringTie.fasta -g ~/Assignment6/DB/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa StringTie.gtf

#Assess the genome guided assembly quality using BUSCO
cd ~/Assignment6/results/4_genome_guided_assembly
busco -i StringTie.fasta -l ~/Assignment6/DB/brassicales_odb10 -o BUSCO_StringTie_brassicales -m transcriptome --cpu 2

##Organise the final results
#Copy QC reports for raw reads (two html files) and clean reads (two html files) into 5_final_assembly
cd ~/Assignment6/results/5_final_assembly
cp ~/Assignment6/results/2_clean_data/Col_leaf_chr4_R1.clean_fastqc.html ./
cp ~/Assignment6/results/2_clean_data/Col_leaf_chr4_R2.clean_fastqc.html ./

#Copy de novo assembly results (including one fasta file and one GFF3 file) into 5_final_assembly
cp ~/Assignment6/results/3_denovo_assembly/Trinity.gff3 ./
cp ~/Assignment6/results/3_denovo_assembly/Col_leaf_chr4_trinity/Trinity.fasta ./

#Copy genome guided assembly results (including one GTF file and one fasta file) into 5_final_assembly
cp ~/Assignment6/results/4_genome_guided_assembly/StringTie.fasta ./
cp ~/Assignment6/results/4_genome_guided_assembly/StringTie.gtf ./