#!/bin/bash

#Download the gff3 file (Clupea harengus)
wget http://ftp.ensembl.org/pub/release-100/gff3/clupea_harengus/Clupea_harengus.Ch_v2.0.2.100.gff3.gz

#Unzip the file
gunzip Clupea_harengus.Ch_v2.0.2.100.gff3.gz

#Change the name of the file to ch.gff
mv Clupea_harengus.Ch_v2.0.2.100.gff3 ch.gff

#Grep gff to find the genome build info as the header & Export the results to a file 
egrep '^.+build' ch.gff > clupea_harengus_gff_features.txt

#the code used used to generate the summary (counts) data
echo -e '#cut -f3 -s ch.gff | sort  | uniq -c | sort >> clupea_harengus_gff_features.txt' >> clupea_harengus_gff_features.txt

#Count how many of each feature type there is, sorted in numerical order & Overwrite the file
cut -f3 -s ch.gff | sort  | uniq -c | sort >> clupea_harengus_gff_features.txt 


