#!/usr/bin/env Rscript

source("scripts/FUNCTIONS.R") # loads packages too
runPCA=T

## Input files ##
Exp_table <- read.csv("results/Filtered_kallisto_TPM.csv",
                      header = T)
metadata <- read_delim("metadata/metadata_pca.txt",
                       col_names = c("replicateName", "BioProject", "Layout", "Length_type", "tissue", "Part"),
                        skip = 1) %>%
  filter(replicateName %in% names(Exp_table)) %>%
  mutate(SampleName = tissue)


Baits <- read_delim("metadata/Baits_ensembl_ids.txt")
pct_genes=20
division_factor=100/pct_genes # by how much we need to divide the list of genes to get the right percentage

## Long exp_table ##

Exp_table_long <- Exp_table %>% 
  pivot_longer(cols = !gene_ID, 
               names_to = "replicateName", 
               values_to = "tpm") %>% 
  mutate(logTPM = log10(tpm + 1)) 

Exp_table_log_wide <- Exp_table_long %>% 
  dplyr::select(gene_ID, replicateName, logTPM) %>% 
  pivot_wider(names_from = replicateName, 
              values_from = logTPM, 
              id_cols = gene_ID)

## PCA ##
if (runPCA){
  source("scripts/10-geneSelection-PCA.R")
}

#### Gene co-expression analysis ####
# Average up the reps #
Exp_table_long_averaged_z <- Exp_table_long %>% 
  full_join(metadata, 
            by = "replicateName") %>% 
  group_by(gene_ID, SampleName) %>%
  summarise(mean.logTPM = mean(logTPM),
            mean.TPM = mean(tpm))  %>% 
  group_by(gene_ID) %>% 
  mutate(z.score.TPM = zscore(mean.TPM),
         z.score.logTPM = zscore(mean.logTPM)) %>% 
  ungroup() 

head(Exp_table_long_averaged_z)

# Write_Exp_table_long_z #

split_at <- 400000

if(nrow(Exp_table_long_averaged_z) < split_at) {
  Exp_table_long_averaged_z %>% 
    mutate(group = 0) %>%
    write_delim(file = "results/Exp_table_long_averaged_z_0.tsv",
                quote = "none", append = F,
                col_names = T, delim = "\t")
} else {
  Exp_table_long_averaged_z %>% 
    mutate(group = (row_number() - 1) %/% !! split_at) %>%
    group_split(group) %>%
    map(.f = ~{
      fileName = paste0("results/Exp_table_long_averaged_z_",
                        unique(.x$group), ".tsv")
      write_delim(.x, fileName,
                  quote = "none", append = F,
                  col_names = T, delim = "\t")
    }) 
  
}

head(Exp_table_long_averaged_z)

#### Ranking ####
all_coefVar_and_ranks <- Exp_table_long_averaged_z %>% 
  group_by(gene_ID) %>% 
  summarise(coefVar.logTPM = relative.variation(mean.logTPM),
            coefVar.TPM = relative.variation(mean.TPM),
            relVar.TPM = relative.variation(mean.TPM, coef.var = F),
            relVar.logTPM = relative.variation(mean.logTPM, coef.var = F)
  ) %>% 
  ungroup() %>%
  mutate(
    rank.CVlogTPM = rank(coefVar.logTPM, ties.method = "average"),
    rank.CVTPM = rank(coefVar.TPM, ties.method = "average"),
    rank.varTPM = rank(relVar.TPM, ties.method = "average"),
    rank.varlogTPM = rank(relVar.logTPM, ties.method = "average")
    )

write_delim(all_coefVar_and_ranks,
            "results//all_coefVar_and_ranks.tsv",
            col_names = T, quote = "none",
            delim = "\t", na = '')

high_var_genes_pct <- 
  all_coefVar_and_ranks %>% 
  slice_max(order_by = rank.varTPM, 
            n = round(nrow(Exp_table)/division_factor, digits = 0)
  )

bait_var <- high_var_genes_pct %>% 
  filter(gene_ID %in% Baits$gene_ID) %>% 
  group_by(gene_ID)
right_join(all_coefVar_and_ranks, Baits, by="gene_ID") %>% View

#### Coeficiente of variation distribution ####
# Check where the bait genes fall along the variance distribution
(var_plot = all_coefVar_and_ranks %>% 
  ggplot(aes(x = relVar.TPM, y = rank.varTPM))  +
   geom_hline(yintercept = 
                nrow(all_coefVar_and_ranks) - nrow(high_var_genes_pct)
              ) +
  geom_rect( 
    xmax = max(high_var_genes_pct$relVar.TPM), 
    xmin = min(high_var_genes_pct$relVar.TPM),
    ymax = nrow(all_coefVar_and_ranks),
    ymin = nrow(all_coefVar_and_ranks) - nrow(high_var_genes_pct),
    fill = "dodgerblue2", alpha = .01
  )  +
  geom_line(linewidth = 0.8) +
  geom_vline(
    data = bait_var, aes(xintercept = relVar.TPM), 
    color = "tomato1", linewidth = 0.2, alpha = 1
  ) + 
  labs(y = "Gene rank",
       x = "TPM relative variance",
       # caption = 
       #   paste0("Blue box = top ", pct_genes, "% high var genes.\nRed lines = bait genes.")
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(color = "black"),
    plot.caption = element_text(hjust = 0)
  ) 
)

# ggsave(plot = var_plot,
#        filename = "plots/MainAnalysis/gene_var_distribution.svg", 
#        height = 5, width = 5)
ggsave(plot = var_plot, 
       filename = "plots/MainAnalysis/gene_relative_variance_distribution.svg", 
       height = 5, width = 7)

#### Gene selection ####


write_delim(high_var_genes_pct,
            file = paste0("results//high_coefvar_genes_", 
                          pct_genes, "pct.tsv"),
            col_names = T, delim = "\t", quote = "none"
)



Exp_table_long_averaged_z_high_var <- 
  Exp_table_long_averaged_z %>% 
  filter(gene_ID %in% high_var_genes_pct$gene_ID)

save(Exp_table_long_averaged_z_high_var,
     high_var_genes_pct, 
     file = "RDATA/GeneSelection_objects.RData"
)
