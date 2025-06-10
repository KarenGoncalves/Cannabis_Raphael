#!/bin/sh
#SBATCH --account=def-laboidp
#SBATCH --mem-per-cpu=8G
#SBATCH --time=0:20:00
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/karencgs/Cannabis_sativa/slurms/fastp_%A-%a.out
#SBATCH --array=124

#go to 130

# This script uses a metadata table ($PWD/metadata/metadata.txt) that has the columns:
# 1) Run ID, 2) sample name (no space), 3) library layout (SE or PE), 
# 4) sequencing platform (Illumina or DNB), read length type (< 50 = short; > 50 normal)

module load StdEnv/2020 scipy-stack fastp
SPECIES=Cannabis_sativa
cd  $SCRATCH/${SPECIES}/
mkdir -p $PWD/clean_reads $PWD/fastpReports/

sample_info=($(awk -v r=$SLURM_ARRAY_TASK_ID 'NR == r {print $1, $2, $3, $4, $5}' $PWD/metadata/metadata.txt))
length_required=$( [[ ${sample_info[4]} == "short" ]]  && echo 30 || echo 50)
 
#This allows us to find reads regardless of their origin (NGDC, SRA or ENA) 
MYREADS=($(ls $PWD/RAW_DATA/${sample_info[0]}*))
# The conditional below will standardize the clean reads names

if [[ ${sample_info[2]} == "PE" ]]; then
        MYREAD1=${MYREADS[0]}
        MYREAD2=${MYREADS[1]}
        MYCLEANREAD1=$PWD/clean_reads/${sample_info[0]}_1.fastq
        MYCLEANREAD2=${MYCLEANREAD1/_1/_2}
else
        MYREAD1=${MYREADS}
        MYCLEANREAD1=$PWD/clean_reads/${sample_info[0]}.fastq
fi


if [[ ${sample_info[2]} == "SE" && ${sample_info[3]} == "Illumina" ]]; then
        analysisType=0
elif [[ ${sample_info[2]} == "PE" && ${sample_info[3]} == "Illumina" ]]; then
        analysisType=1
elif [[ ${sample_info[2]} == "SE" && ${sample_info[3]} == "DNB" ]]; then
        analysisType=2
elif [[ ${sample_info[2]} == "PE" && ${sample_info[3]} == "DNB" ]]; then
        analysisType=3
else
        analysisType=Unknown
fi

# If reads are illumina, we don't give the file with adaptors to fastp
case $analysisType in
	0)
		fastp -i ${MYREAD1}\
		 -o ${MYCLEANREAD1}\
		 --qualified_quality_phred 20\
		 --unqualified_percent_limit 30\
		 --cut_front --cut_front_window_size 5\
		 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 15\
		 --length_required ${length_required}\
		 --json $PWD/fastpReports/${sample_info[0]}.json\
		 --thread ${SLURM_CPUS_PER_TASK}
		;;
	1)
		fastp -i ${MYREAD1}\
		 -I ${MYREAD2}\
		 -o ${MYCLEANREAD1}\
		 -O ${MYCLEANREAD2}\
		 --qualified_quality_phred 20\
		 --unqualified_percent_limit 30\
		 --cut_front --cut_front_window_size 5\
		 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 15\
		 --length_required ${length_required}\
		 --json $PWD/fastpReports/${sample_info[0]}.json\
		 --thread ${SLURM_CPUS_PER_TASK}
		;;
	2)
		fastp -i ${MYREAD1}\
		 -o ${MYCLEANREAD1}\
		 --adapter_fasta $HOME/annotation_scripts/adaptors_DNBSEQ-T7.fasta\
		 --qualified_quality_phred 20\
		 --unqualified_percent_limit 30\
		 --cut_front --cut_front_window_size 5\
		 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 15\
		 --length_required ${length_required}\
		 --json $PWD/fastpReports/${sample_info[0]}.json\
		 --thread ${SLURM_CPUS_PER_TASK}
		;;

	3)
		fastp -i ${MYREAD1}\
		 -I ${MYREAD2}\
		 -o ${MYCLEANREAD1}\
		 -O ${MYCLEANREAD2}\
		 --adapter_fasta $HOME/annotation_scripts/adaptors_DNBSEQ-T7.fasta\
		 --qualified_quality_phred 20\
		 --unqualified_percent_limit 30\
		 --cut_front --cut_front_window_size 5\
		 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 15\
		 --length_required ${length_required}\
		 --json $PWD/fastpReports/${sample_info[0]}.json\
		 --thread ${SLURM_CPUS_PER_TASK}
		;;
	"Unknown" )
		echo "Unknown library layout or sequencing platform"
		echo "Column 3 should be PE|SE, received ${sample_info[2]}"
		echo "Column 4 should be Illumina|DNB, received ${sample_info[3]}"
		echo "Stopping here, check your  $PWD/metadata/metadata.txt file and try again"
		;;
esac
