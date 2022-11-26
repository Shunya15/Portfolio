#!/bin/bash

#Change directory
cd 03_alignedData/log/

#How many reads were mapped in total?
egrep '^SN.reads.mapped:' SRR5882797_thaliana.stats| cut --fields 2,3 | sed 's/://' >> ~/Assignment4/Practical_questions_Q3.txt 

#How many reads were mapped as a pair?
egrep '^SN.reads.mapped.and.paired:' SRR5882797_thaliana.stats| cut --fields 2,3 | sed 's/://' >> ~/Assignment4/Practical_questions_Q3.txt 

#How many reads were mapped as a “proper” pair?
egrep '^SN.reads.properly.paired:' SRR5882797_thaliana.stats| cut --fields 2,3 | sed 's/://' >> ~/Assignment4/Practical_questions_Q3.txt 

