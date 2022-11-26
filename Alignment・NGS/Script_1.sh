#!/bin/bash

#make subdirectory
mkdir Refs
#Go to the right directory
cd Refs/

#make subdirectory
mkdir Athaliana
#Go to the right directory
cd Refs/
cd Athaliana/

#Download genome (fasta) for the model plant Arabidopsis thaliana to the subdirectory Refs/Athaliana
wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-51/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
#Download the annotation (gff3) files for the model plant Arabidopsis thaliana to the subdirectory Refs/Athaliana
wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-51/plants/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.51.gff3.gz

#Unzip the files 
gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.51.gff3.gz

#Identify the number of chromosomes in the genome, and write this information to standardoutput
echo 'The following value shows the number of chromosomes in the genome'
cut -f3 -s  Arabidopsis_thaliana.TAIR10.51.gff3 | egrep -c 'chromosome'
#Identify the number of unique genes located in the genome, and send this information to stdout
echo 'The following value shows the number of unique genes located in the genome'
cut -f3 -s  Arabidopsis_thaliana.TAIR10.51.gff3 | egrep -c '^gene'


##Create a genome index for use with bwa
mkdir -p ~/Assignment4/genome
cp ~/Assignment4/Refs/Athaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa ~/Assignment4/genome/
#change back to our main project folder
cd ~/Assignment4
#make subdirectory 
cd genome
bwa index Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -p Arabidopsis_thaliana.51
