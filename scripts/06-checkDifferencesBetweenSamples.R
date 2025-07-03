#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer))
suppressMessages(library(scriptName))

# Load inputs
args <- commandArgs(trailingOnly=T)
if (length(args) != 3) {
	stop("Usage:
	Rscript ", current_filename(), " TPM_matrix metadata_pca distinct1,distinct2"
	)
}

tpm_matrix = read_delim(args[1])

metadata = read_delim(args[2])

colDistinct = str_split_1(args[3], ",")

if (!all(colDistinct %in% names(metadata))) {
	stop("Choose columns from metadata to color the points in the plot")
}

theme_set(theme_classic())
# Transform tpm to log10(tpm+1)
Exp_table_long <- tpm_matrix %>% 
  pivot_longer(cols = !gene_ID, names_to = "replicateName", values_to = "tpm") %>% 
  mutate(logTPM = log10(tpm + 1)) 


Exp_table_log_wide <- Exp_table_long %>% 
  dplyr::select(gene_ID, replicateName, logTPM) %>% 
  pivot_wider(names_from = replicateName, 
              values_from = logTPM, 
              id_cols = gene_ID)

# Compute PCA
# remove all-zero rows
allZeros = rowSums(Exp_table_log_wide[-1]) == 0
data_PCA <- t(Exp_table_log_wide[!allZeros, -1])
plot_pca <- function(dataPCA, metadata, color_cols) {
	my_pca <- prcomp(data_PCA)
	pc_importance <- as.data.frame(t(summary(my_pca)$importance))


	PCA_coord <- my_pca$x[, 1:2] %>% 
	  as.data.frame() %>% 
	  mutate(replicateName = row.names(.)) %>% 
	  inner_join(metadata, 
	            by = "replicateName")
	xlab=paste("PC1 (", pc_importance[1, 2] %>% signif(3)*100, "% of Variance)", sep = "")
	ylab=paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", sep = "")
	for (col in color_cols) {
	 
		base_plot = PCA_coord %>%
		  ggplot(aes(x = PC1, y = PC2)) 
		if (length(PCA_coord[[col]] %>% unique) <= 6) {
			base_plot = 
				base_plot + 
				geom_point(aes(shape = .data[[col]]),
        	                	       size = 3)  
		} else {
			base_plot = 
				base_plot +
	 		 	geom_point(aes(color = .data[[col]]),
		               			size = 3) +
				  scale_color_manual(values = brewer.pal(length(PCA_coord[[col]] %>% unique), "Set1"))
		}

		base_plot + 
		  labs(x = xlab, y = ylab, color=NULL, shape=NULL) +
		  theme(
		#    aspect.ratio = .7,
		    text = element_text(size= 14),
		    axis.text = element_text(color = "black"),
		    legend.position = "bottom"
		  )
		ggsave(paste0("plots/PCA06_", col, ".svg"),
			height=6, width=6, dpi=1200)
	}
}
plot_pca(data_PCA, metadata, colDistinct)
