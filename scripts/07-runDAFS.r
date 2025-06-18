#### Calculate minimum expression threshold
#### Install and load packages ####
#### All packages should be installed already if "prepare_environment.R" was run

pkgs = c("tidyverse", "DESeq2", "devtools", 
         "GenomicFeatures", "ggpubr", "CustomSelection")

for (curPkg in pkgs) suppressMessages(library(curPkg, character.only = T))


#### CustomSelection ####
tpm_matrix <- read_delim("counts/TPM.tsv")
Exp_table <- tpm_matrix %>%
  data.frame

rownames(Exp_table) = Exp_table$Gene_ID
Exp_table = Exp_table %>% dplyr::select(!Gene_ID)

metadata <- read_delim("metadata/metadata_pca.txt", 
                       col_names = T) |>
  dplyr::select(Tissue_group, Run)

cutv = DAFS(tpm = Exp_table)
data.frame(Run = names(cutv),
	   thresholds = unname(cutv)) %>%
	write.table(file = "counts/cutv.tsv",
		    append = F, quote = F, sep = "\t", 
		    row.names = F, col.names = T)
