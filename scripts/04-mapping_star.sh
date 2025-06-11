#!/bin/sh
#SBATCH --cpus-per-task=8
#SBATCH --account=def-desgagne
#SBATCH --time=2:00:00
#SBATCH --mem=15G
#SBATCH --output=/scratch/karencgs/Cannabis_sativa/slurms/STAR_map-%A-%a.out
#SBATCH --array=1-130

# Using STAR version star/2.7.11b
module load StdEnv/2023 star/2.7.11b
NCPUS=8

cd /scratch/karencgs/Cannabis_sativa

inDIR=$PWD/clean_reads
sample_info=($(awk -v r=$(( $SLURM_ARRAY_TASK_ID+1 )) 'NR == r {print $1, $3}' $PWD/metadata/Read_lengths.txt))

if [[ ${sample_info[1]} -lt 70 ]]; then
	indexDIR=$PWD/genome/Cannabis_30
elif [[ ${sample_info[1]} -lt 90 ]]; then
	indexDIR=$PWD/genome/Cannabis_70
elif [[ ${sample_info[1]} -lt 140 ]]; then
	indexDIR=$PWD/genome/Cannabis_90
else
	indexDIR=$PWD/genome/Cannabis_144
fi 

MYREADS=($(ls $PWD/clean_reads/${sample_info[0]}*))
outPrefix=${sample_info[0]}_

STAR --genomeDir $indexDIR/\
 --runThreadN $NCPUS \
 --readFilesIn ${MYREADS[@]}\
 --outFileNamePrefix alignment_STAR/${outPrefix}\
 --outSAMtype BAM SortedByCoordinate \
 --outSAMunmapped Within \
 --outSAMattributes Standard\
 --quantMode GeneCounts
