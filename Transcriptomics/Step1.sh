#!/bin/bash

#activate the conda environment
conda activate transcriptomics_analysis

#make the directory ~/Assignment6 as the parent directory for all downloads and analysis
mkdir Assignment6


#Assignment6/
#├── data
#├── DB
#├── results
#│   ├── 1_QC
#│   ├── 2_clean_data
#│   ├── 3_denovo_assembly
#│   ├── 4_genome_guided_assembly
#│   └── 5_final_assembly
#└── scripts

#Store files at different processing stages in different folders (refer to the figure above)
cd Assignment6
mkdir data DB results scripts
cd results
mkdir 1_QC 2_clean_data 3_denovo_assembly 4_genome_guided_assembly 5_final_assembly

#Download the genomic sequence & annotation of the Arabidopsis thaliana to the subdirectory DB
cd ~/Assignment6/DB/

#Download genome (fasta) for the model plant Arabidopsis thaliana
#File name: Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-51/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

#Unzip the file
#File name: Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

#Download the annotation (gff3) files for the model plant Arabidopsis thaliana to the subdirectory Refs/Athaliana
#Arabidopsis_thaliana.TAIR10.51.gff3.gz
wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-51/plants/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.51.gff3.gz

#Unzip the file
#File name: Arabidopsis_thaliana.TAIR10.51.gff3
gunzip Arabidopsis_thaliana.TAIR10.51.gff3.gz

#Download the BUSCO “brassicales” lineage dataset
wget https://busco-data.ezlab.org/v5/data/lineages/brassicales_odb10.2020-08-05.tar.gz

#Decompress and untar the BUSCO “brassicales” lineage dataset
tar -zxvf brassicales_odb10.2020-08-05.tar.gz

#Build the genome index for GMAP
cd ~/Assignment6/DB
gmap_build -D ./ -d TAIR10_GMAP Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

#Build the genome index for STAR
cd ~/Assignment6/DB
STAR --runThreadN 2 --runMode genomeGenerate --genomeDir ~/Assignment6/DB/TAIR10_STAR125 --genomeFastaFiles Arabidopsis_thaliana.TAIR10.dna.toplevel.fa --sjdbGTFfile Arabidopsis_thaliana.TAIR10.51.gff3 --sjdbOverhang 124 --genomeSAindexNbases 12
