This readme file was generated on June 19th 2025 by Karen Cristine Goncalves dos Santos

GENERAL INFORMATION

Title of Dataset: Cannabis sativa tissue-wise gene co-expression

<provide at least two contacts>
Principal Investigator Information
Name: Isabel Desgagne-Penix
ORCID:
Institution: Université du Québec à Trois-Rivieres
Email: isabel.desgagne-penix@uqtr.ca

Co-investigator Information
Name: Hugo Germain
ORCID:
Institution: Université du Québec à Trois-Rivieres
Email: hugo.germain@uqtr.ca

Main author Contact Information
Name: Raphael Boucher
ORCID:
Institution: Université du Québec à Trois-Rivieres
Email: raphael.boucher@uqtr.ca


Co-author (Data analysis) Contact Information
Name: Karen Cristine Goncalves dos Santos
ORCID:
Institution: Université du Québec à Trois-Rivieres
Email: karen.cristine.goncalves.dos.santos@uqtr.ca

Date of data download: June 2025


SHARING/ACCESS INFORMATION

Licenses/restrictions placed on the data: 

Links to other publicly accessible locations of the data: 
Data downloaded online from NCBI SRA (list of SRA IDs provided in [metadata/SRA_info.csv](./metadata/SRA_info.csv)). The file also contains, when available, the DOI of the publication associated with it.


METHODOLOGICAL INFORMATION
	
Data was downloaded from SRA using sra-toolkit 3.0.0 fasterq-dump and filtered and trimmed for quality control using fastp (version 0.23.4) with the following parameters:\
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


DATA-SPECIFIC INFORMATION FOR: [FILENAME]
<repeat this section for each dataset, folder or file, as appropriate>

Number of variables: 

Number of cases/rows: 

Variable List: <list variable name(s), description(s), unit(s) and value labels as appropriate for each>

Missing data codes: <list code/symbol and definition>

Specialized formats or other abbreviations used: 
