#!/bin/sh

# Folder name
DIR=$SCRATCH/Cannabis_sativa
if [ $# -eq 0 ]
  then
    alignment_folder=$DIR/alignment
  else
    alignment_folder=$1
fi

if [ $# -gt 1 ]
  then
    resultSummary=$2
  else
    resultSummary=$DIR/metadata/alignment_summary.txt
fi

# List of alignment summaries
ALIGNMENT_SUMMARY=($(ls $alignment_folder/*Log.final.out))

# Names of columns of the alignment stats file
echo "Sample N_input_reads N_uniq_mapped pct_uniq_mapped\
 pct_map_multiple_loci pct_map_many_loci\
 pct_unmap_mismatch pct_unmap_short\
 pct_unmap_other pct_chimera" | tr ' ' '\t' > $resultSummary

for i in ${ALIGNMENT_SUMMARY[@]}; do
	# Get stats into a table
	# we need rows 6, 9, 10, 25, 27, 30, 32, 34 and 37 from the alignment summary
	# These contain the number of input reads (6) number of uniquely mapped reads (9) and the % (10), 
	# Then the percentage of reads mapped to multiple loci (25), to too many loci (27)
	# % of reads unmapped due to too many mismatches (30), because they were too short (32) and for other reasons (34)
	# and the % of chimeric reads (37)
	STATS=($(awk -F '\t' 'NR==6 || NR==9 || NR==10 ||\
			NR==25 || NR==27 || NR==30 || NR==32 ||\
			NR==34 || NR==37 { print $2}' ${i} |\
			sed 's/%//g'))
	
	sampleName=$(basename ${i/_Log.final.out})
	
	# Print the sample name and the stats columns to the alignment stats file
	echo ${sampleName} ${STATS[@]} | tr ' ' '\t' >> $resultSummary
	
done
