This readme file was generated on June 19th 2025 by Karen Cristine Goncalves dos Santos

## GENERAL INFORMATION

Title of Dataset: Cannabis sativa tissue-wise gene co-expression

### Authors

- Principal Investigator Information  
  - Name: Isabel Desgagne-Penix
  - ORCID:
  - Institution: Université du Québec à Trois-Rivieres
  - Email: isabel.<!--.example-->desgagne-penix<!--.example-->@uqtr.ca

- Co-investigator Information
  - Name: Hugo Germain
  - ORCID:
  - Institution: Université du Québec à Trois-Rivieres
  - Email: hugo.<!--.example-->germain<!--.example-->@uqtr.ca

- Main author Contact Information
  - Name: Raphael Boucher
  - ORCID:
  - Institution: Université du Québec à Trois-Rivieres
  - Email: raphael<!--.example-->.boucher@<!--.example-->uqtr.ca

Co-author (Data analysis) Contact Information
- Name: Karen Cristine Goncalves dos Santos
- ORCID:
- Institution: Université du Québec à Trois-Rivieres
- Email: karen.<!--.example-->cristine.goncalves.<!--.example-->dos.santos@<!--.example-->uqtr.ca

## SHARING/ACCESS INFORMATION

Licenses/restrictions placed on the data: 

Data downloaded online from NCBI SRA (list of SRA IDs provided in [metadata/SRA_info.csv](./metadata/SRA_info.csv)). The file also contains, when available, the DOI of the publication associated with it. Date of data download: June 2025.


## METHODOLOGICAL INFORMATION

### Dependencies

Softwares used and versions: sra-toolkit/3.0.0 fastp/0.23.4 star/2.7.11b r/4.4.0

R packages used are listed in [metadata/Rsession.info](./metadata/Rsession.info).

Python modules are listed in [metadata/requirements.txt](./metadata/requirements.txt).

### Methodology

Data was downloaded from SRA using sra-toolkit fasterq-dump and filtered and trimmed for quality control using fastp.

Reads were aligned to _C. sativa_’s female genome (downloaded from Ensembl Plants release 61: [genome](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/cannabis_sativa_female/dna/Cannabis_sativa_female.cs10.dna.toplevel.fa.gz) and [GTF annotation](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/gtf/cannabis_sativa_female/Cannabis_sativa_female.cs10.61.gtf.gz)) using STAR.

Given the difference in average read lengths between the datasets, the genome indexing step was taking using four different values for the `--sjdbOverhang` parameter: 30, 70, 90 and 144.

The option `--quantMode` was added to STAR to obtain read counts. Strandness was infered using RSeQC "infer_experiment.py" script, alignment BAM files and the genome annotation in BED format.

Subsequently, read counts (for each replicate, according to their strandness) were combined into a count matrix using R 4.4.0 and the tidyverse package.

For TPM computation, the GTF file was used to compute gene lengths (through exon overlap) (using R packages rtracklayer and GenomicRanges).

PCA was performed for initial visualization of sample distribution.

Filtering of weakly expressed genes was done using the DAFS method for computation of minimum expression threshold.

Subsequently, gene co-expression was computed with the TidyGeneCoEx method, using TPM relative variance as proxy for variability and keeping the top 20% expressed genes most variable for network construction. 

## DATA-SPECIFIC INFORMATION 

### [FILENAME]
<repeat this section for each dataset, folder or file, as appropriate>

Number of variables: 

Number of cases/rows: 

Variable List: <list variable name(s), description(s), unit(s) and value labels as appropriate for each>

Missing data codes: <list code/symbol and definition>

Specialized formats or other abbreviations used: 
