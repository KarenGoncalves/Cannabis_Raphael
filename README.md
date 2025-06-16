Read files (list of runs and information about them is presented in [metadata/SRA_info.csv](./metadata/SRA_info.csv) were downloaded from NCBI’s Sequence Read Archive (SRA)\
and filtered and trimmed for quality control using fastp (version 0.23.4) with the following parameters:\
```
--qualified_quality_phred 20 --unqualified_percent_limit 30\
 --cut_front_window_size 5 --cut_right_window_size 4\
 --cut_right_mean_qual 15 --length_requiered 50
```
for reads from bioproject PRJNA80055 (which were single-end reads of 35 bases) and sample SRR351932 (composed of single-end reads of 38 bases), minimum length post-trimming was set to 30 bases.

Reads were aligned to _C. sativa_’s female genome (downloaded from Ensembl Plants release 61: [genome](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/cannabis_sativa_female/dna/Cannabis_sativa_female.cs10.dna.toplevel.fa.gz) and [GTF annotation](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/gtf/cannabis_sativa_female/Cannabis_sativa_female.cs10.61.gtf.gz)) using STAR (version 2.7.11b).

Given the difference in average read lengths between the datasets, the genome indexing step was taking using four different values for the `--sjdbOverhang` parameter: 30, 70, 90 and 144.

The option `--quantMode` was added to STAR to obtain read counts. Strandness was infered using RSeQC "infer_experiment.py" script, alignment BAM files and the genome annotation in BED format.

Subsequently, read counts (for each replicate, according to their strandness) were combined into a count matrix using R 4.4.0 and the tidyverse package.
