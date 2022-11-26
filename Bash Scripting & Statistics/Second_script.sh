#!/bin/bash

INFILE=$1

# Check the file has the suffix .fa or .fasta
SUFFIX=$(echo ${INFILE} | sed -r 's/.+(fasta|fa)$/\1/')
if [ ${SUFFIX} == "fa" ] || [ ${SUFFIX} == "fasta" ]; then
  echo File has the suffix ${SUFFIX}
else
  echo File does not have the suffix 'fa' or 'fasta'. Exiting with error.
  exit 1
fi
 

# Define the output file by changing the suffix to .locations
OUTFILE=${INFILE%.${SUFFIX}}.locations
echo Output will be written to $OUTFILE

#Add a header 
echo -e '#ncrna sequences from Drosophila melanogaster in fasta format' > ${OUTFILE} 

#Add column names
echo -e 'gene_id\tchromosome\tstart\tstop\tstrand\tgene_biotype' >> ${OUTFILE}


# Get the header lines which correspond to chromosomes, then collect the
# gene id, chromosome, start and end and write to the output file
egrep '^>.+chromosome' ${INFILE} | \
  sed -r 's/.+BDGP6:(.*):([0-9]+):([0-9]+):(.+)gene:([^ ]+).+gene_biotype:([^ ]+) .+/\5\t\1\t\2\t\3\t\4\t\6/g' \
  >> ${OUTFILE}

echo Done

