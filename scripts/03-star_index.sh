#!/bin/sh
#SBATCH --cpus-per-task=8
#SBATCH --account=def-laboidp
#SBATCH --time=3:00:00
#SBATCH --mem=30G
#SBATCH --output=/scratch/karencgs/Cannabis_sativa/slurms/STAR_idx-%j.out

# Using STAR version star/2.7.11b
module load StdEnv/2023 star/2.7.11b
NCPUS=8

# Start by downloading genomes
cd /scratch/karencgs/Cannabis_sativa/genome
# Fasta
#wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/cannabis_sativa_female/dna/Cannabis_sativa_female.cs10.dna.toplevel.fa.gz
# GTF
#wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/gtf/cannabis_sativa_female/Cannabis_sativa_female.cs10.61.gtf.gz

genomeIdxDIR=Cannabis
genomeFastaFiles=(Cannabis_sativa_female.cs10.dna.toplevel.fasta)
genomeGTFFile=(Cannabis_sativa_female.cs10.61.gtf)
readLength=(30 70 90 144)

for i in ${readLength[@]}; do
 mkdir ${genomeIdxDIR}_$i

 STAR\
  --runThreadN ${NCPUS}\
  --runMode genomeGenerate\
  --genomeDir ${genomeIdxDIR}_$i\
  --genomeFastaFiles ${genomeFastaFiles}\
  --sjdbGTFfile ${genomeGTFFile}\
  --sjdbOverhang ${i}

done
