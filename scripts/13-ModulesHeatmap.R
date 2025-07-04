#!/usr/bin/env Rscript

#### Heatmap of modules ####

##### Set environment ######
gc()
source("scripts/FUNCTIONS.R") # loads packages too

Baits = read_delim("metadata/Baits_ensembl_ids.txt")
funct_anno = read_delim("results/Annotation_expressed_genes.txt")

# Order tissues the way they should appear in plot
tissue_order <- c("Trichome", "Flower", "Meristem", "Leaf", "Petiole", "Stem", "Root")

metadata <-
  read_delim("metadata/metadata_pca.txt", 
             col_names = T, delim = "\t") %>%
  mutate(TissueOrdered = factor(tissue, levels = tissue_order))

##### Load correlation results #####
###### Find r cutoff and resolution ####### 
# File names contain the correlation cutoff and the resolution
r_cutoffs <- grep("MainAnalysis_Cor",
                  list.dirs("RDATA"), value = T) %>%
  gsub(pattern = "RDATA/MainAnalysis_Cor(0.\\d+)",
       replacement = "\\1")

resolutions <-  
  sapply(r_cutoffs, \(x) {
    list.files(path = paste0("RDATA/MainAnalysis_Cor", x),
               pattern = "ModuleData") %>%
      gsub(pattern = ".+_res([0-9\\.]+).RData", 
           replacement = "\\1") %>% unique
  })

cat("Resolution =", resolutions, "; Cut off = ", r_cutoffs)
 
###### Load files ######
for (cur_rCutoff in r_cutoffs) {
  for (curResolution in resolutions[1]) {
    rdata = paste0("RDATA/MainAnalysis_Cor", cur_rCutoff, 
                   "/ModuleData_forPlot",
                   "_res", curResolution, ".RData")
    load(rdata)
    # annotationPath <- paste0("results/Main_analysis_annotations_Cor",
    #                          cur_rCutoff, "/Resolution",
    #                          curResolution)
    # dir.create(annotationPath, recursive = T)
    # sapply(Expr_averaged_z_high_var_modules$module %>% unique,
    #        \(curModule) {
    #          Expr_averaged_z_high_var_modules %>%
    #            dplyr::select(!all_of(c("tissue", "mean.logTPM", "z.score.TPM"))) %>%
    #            unique %>%
    #            filter(module == curModule) %>% 
    #            write_delim(
    #              file = paste0(annotationPath, "/",
    #                            "Module_", curModule, ".tsv"),
    #              delim = "\t", na = "", append = F,
    #              quote = "none",
    #              col_names = T
    #            )
    #        })

    
    heatmap_data = 
     modules_mean_z %>%
      inner_join(metadata %>% select(!replicateName),
                 by = c("tissue", "Part")
      )
    
    # genes per module
    genes_per_module = Expr_averaged_z_high_var_modules %>% 
      group_by(module) %>% 
      summarize(gene_ID = unique(gene_ID)) %>%
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
           caption = paste0("r threshold = ", cur_rCutoff, 
                            "\nresolution = ", curResolution)
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
    plotName = paste0("plots/MainAnalysis/Heatmap_Modules_",
                      "r", cur_rCutoff, "_res", curResolution, 
                      "highVarTop20pct.svg")
    ggsave(plotName, height = 7, width = 8)
  }
}

left_join(Baits, my_network_modules, by = "gene_ID") %>%
  arrange(module) %>% 
  write_delim(file = "results/Baits_in_network.tsv",
              quote = "none", append=F, delim="\t")
