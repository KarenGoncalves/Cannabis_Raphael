#!/bin/sh
#SBATCH --account=def-laboidp
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=10G
#SBATCH --output=/scratch/karencgs/Cannabis_sativa/hmm_%j.out

THREADS=6

module load StdEnv/2023 hmmer/3.4

cd $SCRATCH/Cannabis_sativa

db_dir=/project/def-desgagne/karencgs/trinotate_data
THREADS=16


hmmscan\
  --cpu ${THREADS}\
  --domtblout results/Modules_4689_domains.out\
  --tblout results/Modules_4689_seqHits.out\
  --pfamtblout results/Modules_4689_pfam.out\
  --cut_nc\
  $db_dir/Pfam-A.hmm\
  genes_of_interest.pep
