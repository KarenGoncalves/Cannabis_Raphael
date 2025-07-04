#!/usr/bin/env Rscript

pkgs=c("tidyverse", "igraph", "ggraph",
       "readxl", "patchwork", "RColorBrewer",
       "viridis", "parallel", "doParallel",
       "furrr", "future"
)

for (pkg in pkgs) suppressMessages(library(pkg, character.only = T))

### ggplot theme ###
theme_optimization = 
  theme_classic() +
  theme(text = element_text(size = 14),
        axis.text = element_text(color = "black"),
        axis.title.y.right = element_text(color = "violetred2"),
        axis.title.y.left = element_text(color = "dodgerblue2")
  )

# Z-score function #
zscore <- function(numeric_vector) {
  (numeric_vector - mean(numeric_vector)) / 
    sd(numeric_vector)
}

# relative variation
relative.variation <- 
  function(numeric_vector, na.rm=T, coef.var = T) {
    if (coef.var) {
      return(
        sd(numeric_vector, na.rm=na.rm) / 
          mean(numeric_vector, na.rm=na.rm)
      )
    } else {
      return(
        var(numeric_vector, na.rm=na.rm) / 
          mean(numeric_vector, na.rm=na.rm)
      )
    }
  }

# optimize_resolution function #
# Method to find best resolution for cluster detection
optimize_resolution <- function(network, resolution, minGenes){
  modules = network %>% 
    cluster_leiden(
      resolution = resolution,
      objective_function = "modularity")
  
  parsed_modules = data.frame(
    gene_ID = names(membership(modules)),
    module = as.vector(membership(modules)) 
  )
  
  # counts the number of modules found with more than 5 genes
  num_module_minGenes = parsed_modules %>% 
    group_by(module) %>% 
    summarize(n = n()) %>% 
    arrange(-n) %>% 
    filter(n >= minGenes) %>% 
    nrow() %>% 
    as.numeric()
  
  # Counts the number of genes contained in all modules with more than 5 genes
  num_genes_contained = parsed_modules %>% 
    group_by(module) %>% 
    summarize(n = n()) %>% 
    arrange(-n) %>% 
    filter(n >= minGenes) %>% 
    ungroup() %>% 
    summarise(sum = sum(n)) %>% 
    as.numeric()
  
  return(c(num_module_minGenes, num_genes_contained))
  
} 

get_optimization_table <- function(network, 
                                   min_res, 
                                   max_res, 
                                   step,
                                   minGenes) {
  future_map(.options = furrr_options(seed = TRUE),
             .x = seq(from = min_res, to = max_res, by = step),
             \(x) optimize_resolution(network = network, 
                                      resolution = x,
                                      minGenes)) %>% 
    lapply(\(x) data.frame(V1 = x[1], V2 = x[2])) %>%
    list_rbind %>% 
    mutate(resolution = seq(from = min_res, to = max_res, by = step)) %>% 
    rename(num_module = V1,
           num_contained_gene = V2)
}

select = dplyr::select


cluster_leiden_mod <- function(network_data, r) {
  return(cluster_leiden(network_data, 
                        resolution_parameter = r, 
                        objective_function = "modularity"
  ))
}

get_T_value <- function(r, nreps) {
  t = r*sqrt( (nreps-2) / (1-r^2) )
  return(t)
}
get_p_value <- function(t, nreps) {
  pvalue = case_when(
    t > 0 ~ pt(t, df = nreps-2, lower.tail = F),
    t <=0 ~ pt(t, df = nreps-2, lower.tail = T)
  )
  return(pvalue)
}

heatmap_theme <- theme_classic() +
  theme(text = element_text(size = 14),
        axis.text.y = element_text(color = "black"),
        #axis.text.x = element_blank(),
        legend.position = "top",
        strip.text = element_blank(),
        panel.spacing = unit(0, "lines")
  )
annotation_theme <- theme_void() +
  theme(axis.title.x = element_text(vjust = 8))

##### Heatmap of all genes in the module #####
heatmap_all_genes <- function(cor_matrix, 
                              Exp_table_long_averaged_z,
                              caption, title) {
  gene_clustering <- as.dist(1 - cor_matrix) %>%
    hclust(method = "complete")
  Exp_table_long_averaged_z$gene_clustering_order <- 
    factor(Exp_table_long_averaged_z$gene_ID, 
           levels = rownames(cor_matrix)[
             gene_clustering$order]
    )
  breaks = with(Exp_table_long_averaged_z,
                c(min(z.score), 0, max(z.score))) 
  heatmap_genes <-
    Exp_table_long_averaged_z %>%
    ggplot(aes(x = Sample_name, y = gene_clustering_order)) +
    geom_tile(aes(fill = z.score), color = NA) +
    scale_fill_gradient2(mid = "white",
                         high = "#67001F",
                         low = "#053061",
                         breaks = breaks,
                         labels = breaks %>% round(2)
    ) +
    labs(x = NULL, y = "Genes", fill = "z-score",
         title = title,
         caption = paste0(caption, "\nGenes in module: ",
                          length(gene_clustering$labels))
    )
  return(heatmap_genes)
  
}
