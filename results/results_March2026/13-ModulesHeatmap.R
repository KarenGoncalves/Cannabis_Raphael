#!/usr/bin/env Rscript

#### Heatmap of modules ####

##### Set environment ######
gc()
source("scripts/FUNCTIONS.R") # loads packages too
Baits=readxl::read_excel("BaitGenes_RNA_CoExprAnalysis.xlsx", 
                         sheet = "Bait_genes_Ensembl_ID") %>% 
  rename(#new_name=old_name
    Description = Name,
    gene_ID = `Ensembl Plant ID`)
args = commandArgs(trailingOnly=T)
cores = ifelse(length(args) != 1, detectCores(), as.numeric(args[1]))
basedir = "results/results_March2026/"
funct_anno = read_delim("results/Annotation_expressed_genes.txt")

# Order tissues the way they should appear in plot
tissue_order <- c("Trichome", "Flower", "Meristem", "Leaf", "Stem", "Root")

metadata <-
  read_delim("metadata/metadata_pca.txt", 
             col_names = T, delim = "\t") %>%
  mutate(TissueOrdered = factor(tissue, levels = tissue_order))

##### Load correlation results #####
###### Find r cutoff and resolution ####### 
# File names contain the correlation cutoff and the resolution
r_cutoffs <- grep("MainAnalysis_Cor",
                  list.dirs(basedir), value = T) %>%
  gsub(pattern = ".+MainAnalysis_Cor(0.\\d+)",
       replacement = "\\1")

resolutions <-  
  sapply(r_cutoffs, \(x) {
    list.files(path = paste0(basedir, "MainAnalysis_Cor", x),
               pattern = "ModuleData") %>%
      gsub(pattern = ".+_res([0-9\\.]+).RData", 
           replacement = "\\1") %>% unique
  })

cat("Resolution =", resolutions, "; Cut off = ", r_cutoffs)
 
###### Load files ######
for (cur_rCutoff in r_cutoffs) {
  for (curResolution in resolutions[1]) {
    rdata = paste0(basedir, "/MainAnalysis_Cor", cur_rCutoff, 
                   "/ModuleData_forPlot",
                   "_res", curResolution, ".RData")
    load(rdata)

    heatmap_data = 
     modules_mean_z %>%
      inner_join(metadata %>% select(!replicateName),
                 by = c("tissue", "Part")
      )
    
    # genes per module
    genes_per_module = Expr_averaged_z_high_var_modules %>% 
      select(module, gene_ID) %>% 
      unique() %>%
      group_by(module) %>% 
      summarize(n=n())
    
    heatmap_data = full_join(heatmap_data, genes_per_module,
                             by = "module") %>%
      arrange(desc(n)) %>%
      mutate(ordered_modules = factor(module,
                                      levels = unique(.$module))
      )
    
    ### Heatmap of module mean expression z-score
    module_heatmap <- heatmap_data %>% 
      ggplot(aes(x = TissueOrdered,
                 y = ordered_modules)) +
      geom_tile(aes(fill = mean.z), color = "grey80") +
      scale_fill_gradient2(mid = "white",
                           high = "#67001F",
                           low = "#053061",
                           breaks = c(-1.5, 0, 1.5), 
                           labels = c("< -1.5", "0", "> 1.5")) +
      labs(x = NULL, y = "Module", fill = "z-score",
           # caption = paste0("r threshold = ", cur_rCutoff, 
           #                  "\nresolution = ", curResolution)
      ) +
      heatmap_theme
    
    module_heatmap_nGenes = heatmap_data %>%
      select(ordered_modules, n) %>%
      unique %>%
      ggplot(aes(y = ordered_modules, x = "",
                 label = n)) +
      geom_tile(fill = "white", color = "white") +
      geom_text() + labs(x = "Genes\nin\nmodule") +
      annotation_theme
    
    wrap_plots(module_heatmap, 
               module_heatmap_nGenes, 
               ncol = 2, 
               widths = c(1, 0.08))
    plotName = paste0(basedir, "/Heatmap_Modules_",
                      "r", cur_rCutoff, "_res", curResolution, 
                      "highVarTop20pct.svg")
    ggsave(plotName, height = 7, width = 8)
  }
}

left_join(Baits, my_network_modules, by = "gene_ID") %>%
  arrange(module) %>% 
  write_delim(file = paste0(basedir, "/Baits_in_network.tsv"),
              quote = "none", append=F, delim="\t")

left_join(Baits, my_network_modules, by = "gene_ID") %>%
  arrange(module) %>% View()
