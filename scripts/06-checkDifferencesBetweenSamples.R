#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer))
suppressMessages(library(scriptName))

# Load inputs
args <- commandArgs(trailingOnly=T)
if (length(args) != 2) {
	stop("Usage:
	Rscript", current_filename(), "TPM_matrix metadata_pca"
	)
}

tpm_matrix = read_delim(args[1])

metadata = read_delim(args[2])

# Transform tpm to log10(tpm+1)
Exp_table_long <- tpm_matrix %>% 
  pivot_longer(cols = !Gene_ID, names_to = "Run", values_to = "tpm") %>% 
  mutate(logTPM = log10(tpm + 1)) 


Exp_table_log_wide <- Exp_table_long %>% 
  dplyr::select(Gene_ID, Run, logTPM) %>% 
  pivot_wider(names_from = Run, 
              values_from = logTPM, 
              id_cols = Gene_ID)

# Compute PCA
# remove all-zero rows
allZeros = rowSums(Exp_table_log_wide[-1]) == 0
my_pca <- prcomp(t(Exp_table_log_wide[!allZeros, -1]))
pc_importance <- as.data.frame(t(summary(my_pca)$importance))


PCA_coord <- my_pca$x[, 1:10] %>% 
  as.data.frame() %>% 
  mutate(Run = row.names(.)) %>% 
  inner_join(metadata, 
            by = "Run")


(PCA_by_Tissue <- PCA_coord %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Tissue_group, color = Tissue_group,
                 shape = Layout),  
               size = 3) +
  scale_color_manual(values = brewer.pal(length(PCA_coord$Tissue_group %>% unique), "Set1"),) +
    labs(x = paste("PC1 (", pc_importance[1, 2] %>% signif(3)*100, "% of Variance)", sep = ""),
         y = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
         fill = NULL, color = NULL, shape = NULL) +
  facet_wrap(~Project) +
  theme_grey() +
  theme(
#    aspect.ratio = .7,
    text = element_text(size= 14),
    axis.text = element_text(color = "black"),
    # legend.position = "bottom"
  )
)
dir.create("plots/MainAnalysis", recursive = T)
ggsave("plots/MainAnalysis/PCA_project_tissue_layout.svg", 
       height = 10, width=8, dpi=1200)

#### Repeat pca without PRJNA560453 ####
remove_560453 <- 
  metadata$Run[
    metadata$Project != "PRJNA560453"
  ]
my_pca <- 
  prcomp(t(Exp_table_log_wide[!allZeros, remove_560453]))
pc_importance <- as.data.frame(t(summary(my_pca)$importance))


PCA_coord <- my_pca$x[, 1:10] %>% 
  as.data.frame() %>% 
  mutate(Run = row.names(.)) %>% 
  inner_join(metadata, 
             by = "Run")

shoot_sample <- PCA_coord %>%
  filter(Tissue == "Shoot")
(PCA_by_Tissue <- PCA_coord %>%
    ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Tissue_group, color = Tissue_group,
                 shape = Layout),  
               size = 3) +
  scale_color_manual(values = brewer.pal(length(PCA_coord$Tissue_group %>% unique), "Set1"),) +
    labs(x = paste("PC1 (", pc_importance[1, 2] %>% signif(3)*100, "% of Variance)", sep = ""),
         y = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
         fill = NULL, color = NULL, shape = NULL) +
  facet_wrap(~Project) +
    theme_grey() +
    theme(
#      aspect.ratio = .7,
      text = element_text(size= 14),
      axis.text = element_text(color = "black"),
      # legend.position = "bottom"
    ) +
    annotate(geom="text", 
             x = shoot_sample$PC1,
             y = shoot_sample$PC2+5,
             label = "Shoot")
)
ggsave("plots/MainAnalysis/PCA_without560453.svg",
       height = 10, width=8, dpi=1200)
