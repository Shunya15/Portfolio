#!/bin/bash

#Get the sequences
wget https://university-of-adelaide-bx-masters.github.io/BIOTECH-7005-BIOINF-3000/Practicals/evolutionary_prac/bovidae_118_mtDNA.fa

#Reduce the size of the dataset
./subset -id aXXXXXXX -n 50 -in bovidae_118_mtDNA.fa > bovidae_50_mtDNA.fa

#Multiple alignment
#to construct a phylogenetic tree
mafft bovidae_50_mtDNA.fa > bovidae_50_mtDNA-mafft.fa

#Remove non-conserved blocks
sed -e 's/^>[^ ]\+ \([^ ]\+\) \([^ ]\+\).*$/>\1_\2/g' bovidae_50_mtDNA-mafft.fa > bovidae_50_mtDNA-named.fa

#examine the alignments more effectively we will convert the format from FASTA to NEXUS format
seqmagick convert --output-format nexus --alphabet dna bovidae_50_mtDNA-named.fa bovidae_50_mtDNA-named.nex

#Run Mr Bayes
mb
execute
lset nst=<num>
mcmc ngen=50000 relburnin=no burnin=0 filename=model-<num>

#Examine convergence
sump

#Examine trees
sumt

#Bayesean Trees with whole mitochondrial genomes
./fetch -email <youremailaddress> -query "mitochondrion[All Fields] AND \"Metatheria\"[Organism] AND \"complete genome\"[All fields] AND \"RefSeq\"[All 
fields]" -out metatheria_mtDNA.fa


