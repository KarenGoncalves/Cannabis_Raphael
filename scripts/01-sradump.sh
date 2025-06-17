#!/bin/sh

module load StdEnv/2020  gcc/9.3.0 sra-toolkit/3.0.0
cd $SCRATCH/Cannabis_sativa

#accessions=($(cut -f 1 metadata/metadata.txt | tail -n +125))

acc=($(awk 'NR < 66 || NR > 123 {print $1}' metadata/metadata.txt))

#mkdir RAW_DATA/

for accession in ${acc[@]}; do
 fasterq-dump --seq-defline '@$ac.$sn./$ri'\
 --split-files --outdir RAW_DATA/\
 -x ${accessions} -e 8 & done
