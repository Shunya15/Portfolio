#!/bin/bash

#make subdirectory
mkdir 01_rawData/
#Go to the right directory
cd 01_rawData/

#make subdirectory
mkdir FastQC
mkdir log
mkdir fastq
#Go to the right directory
cd fastq/

#Download the sequencing data 
wget https://universityofadelaide.box.com/shared/static/egl3n16r0ziaxlvbs9074xqd1liktnuz.gz
wget https://universityofadelaide.box.com/shared/static/g2ly4kzz1blus5juy426i37zl45o38pu.gz

#Change the name of the file to SRR5882797_10M_1.fastq.gz & SRR5882797_10M_2.fastq.gz, respectively
mv egl3n16r0ziaxlvbs9074xqd1liktnuz.gz SRR5882797_10M_1.fastq.gz
mv g2ly4kzz1blus5juy426i37zl45o38pu.gz SRR5882797_10M_2.fastq.gz

#Go to Assignment4
cd ~/Assignment4

#make subdirectory 
mkdir 02_trimmedData
#Go to mkdir 02_trimmedData
cd 02_trimmedData/

#make subdirectory 
mkdir FastQC
mkdir log
mkdir fastq


## Define the key directories
PROJROOT=/home/student/Assignment4
RAWFQ=${PROJROOT}/01_rawData/fastq
RAWQC=${PROJROOT}/01_rawData/FastQC
TRIMFQ=${PROJROOT}/02_trimmedData/fastq
TRIMQC=${PROJROOT}/02_trimmedData/FastQC
TRIMLOG=${PROJROOT}/02_trimmedData/log

## Check the project root exists
if [[ -d ${PROJROOT} ]]; then
  echo -e "Found ${PROJROOT}\n"
else
  echo -e "${PROJROOT} not found.\nExiting with Error code 1"
  exit 1
fi

## Check all directories exist for the raw data
if [[ -d ${RAWFQ} ]] && [[ -d ${RAWQC} ]]; then
  echo -e "Found ${RAWFQ}\nFound ${RAWQC}\n"
else
  echo -e "Raw data directories not found.\nExiting with Error code 2"
  exit 2
fi

## Check all directories exist for the trimmed data
if [[ -d ${TRIMFQ} ]] && [[ -d ${TRIMQC} ]] && [[ -d ${TRIMLOG} ]]; then
  echo -e "Found ${TRIMFQ}\nFound ${TRIMQC}\nFound ${TRIMLOG}\n"
else
  echo -e "Trimmed data directories not found.\nExiting with Error code 3"
  exit 3
fi

## Run FastQC on the raw data
fastqc -t 2 -o ${RAWQC} ${RAWFQ}/*gz

## Trim the data
echo -e "Running cutadapt\n"
cutadapt \
  -m 35 \
  -q 30 \
  -a CTGTCTCTTATACACATCT \
  -A CTGTCTCTTATACACATCT \
  -o ${TRIMFQ}/SRR5882797_10M_1.fastq.gz \
  -p ${TRIMFQ}/SRR5882797_10M_2.fastq.gz \
  ${RAWFQ}/SRR5882797_10M_1.fastq.gz ${RAWFQ}/SRR5882797_10M_2.fastq.gz > \
  ${TRIMLOG}/cutadapt.log
  
## Run FastQC on the trimmed data
fastqc -t 2 -o ${TRIMQC} ${TRIMFQ}/*gz

#make "log" folder and "bam" folder in 03_alignedData
cd ~/Assignment4/
mkdir 03_alignedData
cd 03_alignedData/
mkdir log
mkdir bam
cd ~/Assignment4

#Align paired-end reads to the genome index using bwa mem
bwa mem -t 2 genome/Arabidopsis_thaliana.51 02_trimmedData/fastq/SRR5882797_10M_1.fastq.gz 02_trimmedData/fastq/SRR5882797_10M_2.fastq.gz | samtools view -bhS -F4 - > 03_alignedData/bam/SRR5882797_thaliana.bam


#find out information about the alignments using samtools stats
cd ~/Assignment4
samtools stats 03_alignedData/bam/SRR5882797_thaliana.bam > \
  03_alignedData/log/SRR5882797_thaliana.stats
  
#Sort and index the bam file
mkdir ~/Assignment4/03_alignedData/sorted_bam
cd ~/Assignment4/03_alignedData/
samtools sort bam/SRR5882797_thaliana.bam -o sorted_bam/SRR5882797_thaliana.bam



