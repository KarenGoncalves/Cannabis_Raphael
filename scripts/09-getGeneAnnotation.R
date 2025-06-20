#!/usr/bin/env Rscript

# Get annotation for expressed genes
suppressMessages(library(tidyverse))
suppressMessages(library(biomaRt))

## Expressed genes
genes_expressed <- read_delim("results/Filtered_TPM.csv")[["gene_ID"]]

mart <- useMart(host="https://plants.ensembl.org", 
                biomart="plants_mart",
                dataset = "csfemale_eg_gene")
getBM(mart = mart,
        attributes = c("ensembl_gene_id", "description", 
                       "go_id", "name_1006", "namespace_1003"),
        filters = "ensembl_gene_id",
        values = genes_expressed) %>%
  write_delim("results/Annotation_expressed_genes.txt")
