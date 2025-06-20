#!/usr/bin/env Rscript

#### Install and load packages ####
pkgs = c("tidyverse", "ggpubr", "DESeq2", "GenomicFeatures")

for (curPkg in pkgs) suppressMessages(library(curPkg, character.only = T) )

##### VARIABLES #####
DAFS_ran = T
#### CustomSelection ####
metadata <- read_delim("metadata/metadata_pca.txt", 
		       col_names = c("replicateName", "BioProject", "Layout", "Length_type", "tissue", "Part"),
			skip = 1) %>%
	dplyr::select(tissue, replicateName)

replicatesPerTissue <- table(metadata$tissue)

Exp_table <- 
  read_delim("counts/TPM.tsv", col_names = T, delim = " ") %>%
  dplyr::select(all_of(metadata$replicateName), gene_ID) |>
  data.frame()

rownames(Exp_table) = Exp_table$gene_ID
Exp_table = Exp_table %>% dplyr::select(!gene_ID)


if (!DAFS_ran) {
	source("scripts/07-runDAFS.r")
}
cutv_table = read_delim("counts/cutv.tsv", skip=1, col_names = c("Replicates", "thresholds")) %>%
  filter(Replicates %in% names(Exp_table))
cutv = cutv_table$thresholds
names(cutv) = cutv_table$Replicates

# We check, for each replicate, if the gene passed its expression threshold
passCutOff_by_rep = sapply(names(cutv), \(Replicate) {
	Exp_table[[Replicate]] > cutv[Replicate]^2
}) %>%  as.data.frame
rownames(passCutOff_by_rep) = rownames(Exp_table)

# Then, we check, for each sample, if the gene passed the threshold in all replicates
new_Exp_table = 
  sapply(unique(metadata$tissue), simplify = F, USE.NAMES = F,
         \(sampleName) {
           # get the names of the replicates for a given tissue
           replicates = subset(metadata, tissue == sampleName)$replicateName
           print(replicates)
           # Check, for each gene, whether the level of expression is 
           # greater than the threshold for all replicates of that tissue 
           passedCutoff = apply(passCutOff_by_rep[, replicates], 1, all)
           
           # Then for each replicate
           sapply(replicates, \(replicate) {
             # And each gene
             sapply(1:length(passedCutoff), \(gene) {
               # if it passed the cutoff, assign it the expression level from the expression table
               # otherwise, assign it 0
               ifelse(passedCutoff[gene], Exp_table[gene, replicate], 0)
             })
           }) %>% data.frame
         }) %>% list_cbind() %>%
  data.frame(row.names = Exp_table %>% rownames())

# we keep the genes that passed the threshold in at least one sample
Exp_filtered = new_Exp_table[rowSums(new_Exp_table) != 0,] 

Exp_filtered %>% 
	mutate(gene_ID = rownames(Exp_filtered)) %>%
	write_csv("results/Filtered_TPM.csv",
		  quote = "none", append = F,
		  col_names = T)
	
