#!/bin/sh
#SBATCH --account=def-desgagne
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --output=/scratch/karencgs/Cannabis_sativa/slurms/featureCounts_%j.out

# Softwares to use
module load StdEnv/2020 gcc/9.3.0 subread

# Folders where files are and will be stored
DIR=/scratch/karencgs/Cannabis_sativa
GTF=$DIR/genome/Cannabis_sativa_female.cs10.61.gtf
bam_Path=$DIR/alignments
outPutDir=$DIR/counts

mkdir $outPutDir/

# Use the metadata table to separate single and paired-end data
SE=$( awk -F "\t" '$3 == "SE" {print $1}' $DIR/metadata/metadata.txt |\
 while read i; do echo $bam_Path/${i}*bam; done | xargs)
PE=$( awk -F "\t" '$3 == "PE" {print $1}' $DIR/metadata/metadata.txt |\
 while read i; do echo $bam_Path/${i}*bam; done | xargs)


# featureCounts by default considers the data is single-end
# -a is the annotation file, -o is where the output table will be saved
# -t and -g by default are 'exon' and 'gene_id', respectively
# -T is the number of threads. 
# It takes between 5 secs and 5 mins per file


featureCounts\
 -a $GTF\
 -o $outPutDir/Cannabis_SE.txt\
 -t 'exon' -g 'gene_id'\
 -T 12\ 
 ${SE}

# for paired-end, both -p (paired-end) and --countReadPairs arguments are needed
featureCounts\
 -a $GTF\
 -o $outPutDir/Cannabis_PE.txt\
 -t 'exon' -g 'gene_id'\
 -T 12 -p --countReadPairs\ 
 ${PE}

# The output table from featureCounts has a comment in the first line. Remove it with tail -n +2
# It also contains 5 extra columns we don't need:
# Chr, Start, End, Strand, Length . Remove with cut
# Finally, the column names are filenames. Use sed to keep only the read number.
for i in $outPutDir/Cannabis_*.txt; do
	tail -n +2 $i |\
	cut -f 1,7- |\
	sed -E "s@$bam_Path/(SRR[0-9]+)_Aligned.sortedByCoord.out.bam@\1@g" >\
	${i/.txt/_formatted.txt}
done

cd $DIR

module load r/4.3.1

R --no-save --vanilla -e '
suppressMessages(library(tidyverse))

SE_table = read_delim("counts/Cannabis_SE_formatted.txt")
PE_table = read_delim("counts/Cannabis_PE_formatted.txt")

full_join(SE_table, PE_table, by = "Geneid") %>% write_table("counts/Cannabis_fullCounts.txt")
'
