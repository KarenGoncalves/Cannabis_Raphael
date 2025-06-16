#!/bin/sh
module load StdEnv/2023 gcc scipy-stack r/4.4.0

# install how_are_we_stranded with: pip install how_are_we_stranded
# it will also download RSeQC, which has the function infer_experiment.py

# If you do not have the bed file, use the RSeQC function gtf2bed --gtf $GTF --bed ${GTF/.gft/.bed}

cd $SCRATCH/Cannabis_sativa
rm metadata/Stranded_data.txt


# For each alingment BAM file, use infer_experiment to get the proportion of reads which have:
# the right-most read aligning to the same strand as the transcript is read (forward, +, first)
# the left-most read aligning to the same strand as the transcript is read (forward, +, first)

for bam in alignment_STAR/*_Aligned.sortedByCoord.out.bam; do\
        run=$(basename ${bam/_Aligned.sortedByCoord.out.bam})
        out=alignment_STAR/${run}_strandness.txt
       infer_experiment.py\
        -i $bam\
        -r genome/Cannabis_sativa_female.cs10.61.bed >\
        ${out};

        layout=$(grep "This is " $out | sed -E 's/This is ([PS])[a-z]+(E)nd Data/\1\2/')
        failed=$(grep "reads failed" $out | sed -E 's/.+: //')
        fwd=$(grep "reads failed" $out -A1 | sed -E 's/.+: //' | tail -n 1)
        rev=$(grep "reads failed" $out -A2 | sed -E 's/.+: //' | tail -n 1)

        echo 'layout="'$layout'"
fwd='$fwd'
rev='$rev'
fwd_percent=fwd / (fwd + rev)
rev_percent=rev / (fwd + rev)

# If there is more forward than reverse, the right-most read was sequenced first
if (fwd_percent > 0.8) {
    cat("First_stranded")
} else if (rev_percent > 0.8) {
# On the opposite case, the right-most was the second read sequenced
    cat("Second_stranded")
} else if (max(fwd_percent, rev_percent) < 0.6) {
# If the two have similar proportions, there was no strand specificity
        cat("Unstranded")
}' > tmp.R
        stranded=$(Rscript tmp.R)
        rm tmp.R
        echo $run $layout $stranded | tr ' ' '\t' >> metadata/Stranded_data.txt
done
