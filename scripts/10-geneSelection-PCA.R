#!/usr/bin/env Rscript
## run pca within script 10-geneSelection.R


my_pca <- prcomp(t(Exp_table_log_wide[, -1]))
pc_importance <- as.data.frame(t(summary(my_pca)$importance))
head(pc_importance, 20)

PCA_coord <- my_pca$x[, 1:10] %>% 
  as.data.frame() %>% 
  mutate(Replicate = row.names(.)) %>% 
  full_join(metadata, by = join_by("Replicate" == "replicateName"))

axis_titles = sapply(1:2, \(x) {
  paste("PC", x, " (", 
        pc_importance[x, 2] %>% signif(3)*100, 
        "% of Variance)", sep = "")
})

PCA_by_tissue <- PCA_coord %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = tissue,
                 shape = BioProject),
             size = 3, alpha = 0.8) +
  scale_fill_manual(values = brewer.pal(n = 8, "Accent")) +
  labs(x = axis_titles[1], 
       y = axis_titles[2],
       color = NULL, shape = NULL) +  
  theme_bw() +
  theme(text = element_text(size= 12),
        legend.text = element_text(size= 10),
        axis.text = element_text(color = "black"),
        #legend.position = "bottom"
  )

PCA_by_tissue

#ggsave("plots/PCA.svg", height = 5, width = 5.5, bg = "white")
ggsave("plots/MainAnalysis/PCA_filtered_data.svg",
       height = 5, width = 8, bg = "white")


PCA_by_part <- PCA_coord %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Part,
                 shape = BioProject),
             size = 3, alpha = 0.8) +
  scale_fill_manual(values = brewer.pal(n = 8, "Accent")) +
  labs(x = axis_titles[1], 
       y = axis_titles[2],
       color = NULL, shape = NULL) +  
  theme_bw() +
  theme(text = element_text(size= 12),
        legend.text = element_text(size= 10),
        axis.text = element_text(color = "black"),
        #legend.position = "bottom"
  )

PCA_by_part

#ggsave("plots/PCA.svg", height = 5, width = 5.5, bg = "white")
ggsave("plots/MainAnalysis/PCA_filtered_parts.svg",
       height = 5, width = 8, bg = "white")
