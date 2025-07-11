#!/bin/sh
#SBATCH --account=def-laboidp
#SBATCH --time=0:5:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=1G
#SBATCH --output=/scratch/karencgs/Cannabis_sativa/blastp_%j.out

THREADS=8
        
module load StdEnv/2023 blast+
cd $SCRATCH/Cannabis_sativa


db_dir=/project/def-desgagne/karencgs/trinotate_data

blastp -query genes_of_interest.pep\
  -db $db_dir/uniprot_sprot.pep\
  -evalue 1e-5\
  -task  blastp\
  -num_threads ${THREADS}\
  -max_target_seqs 1\
  -outfmt '6 std ppos stitle'\
  > results/Modules_4689_blastp.out
        
