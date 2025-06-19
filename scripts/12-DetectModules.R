#!/usr/bin/env Rscript

## Detect gene coexpression modules ##

#### Prepare environment ####
source("scripts/FUNCTIONS.R") # loads packages too
args = commandArgs(trailingOnly=T)
cores = ifelse(length(args) != 1, 1, as.numeric(args[1]))


load("RDATA/GeneSelection_objects.RDATA")
Baits = read_delim("metadata/Baits_ensembl_ids.txt")
metadata <- read_delim("metadata/metadata_pca.txt",
                       col_names = c("replicateName", "BioProject", "Layout", "Length_type", "tissue", "Part"),
                        skip = 1) %>%
  mutate(SampleName = tissue)

## Minimum number of genes to form a module ##
minGenes = 5
# future::plan() - sets parallel session (using multisession())
future::plan(multisession, workers = cores)

#### Find gene co-expression modules ####
# with Leiden algorithm (graph-based method)
edgeTableFiles = 
  paste0("results/Main_analysis_edge_table/",
         list.files(pattern = "edge_table_r",
                    path = "results/Main_analysis_edge_table/")
  )

edge_table_select = sapply(edgeTableFiles, simplify = F, \(x) {
  read_delim(x, col_names = T, delim = "\t")
}) %>% list_rbind %>%
  dplyr::select(!group)

dim(edge_table_select)

# Load gene annotation table
funct_anno = read_delim("results/Annotation_expressed_genes.txt")

# Merge annotation and edge tables
node_table <- data.frame(
  gene_ID = c(edge_table_select$from, 
              edge_table_select$to) %>% unique()
) %>% 
  left_join(funct_anno[, 1:2], # node table cannot have duplicated transcript ids
            # so we keep only the transcript id and description, not the GO annotation
            by = join_by("gene_ID" == "ensembl_gene_id")) %>%
  unique()

head(node_table)
dim(node_table)

# Create graph
my_network <- graph_from_data_frame(
  edge_table_select,
  vertices = node_table,
  directed = F
)

#' Too low resolution leads to forcing genes with 
#' different expression patterns into the same module.
#' Too high resolution leads to many genes not contained in any one module.

optimization_results <-
  future_map(.options = furrr_options(seed = TRUE),
             .x = seq(from = 0.25, to = 5, by = .25),
             \(x) optimize_resolution(network = my_network,
                                      resolution = x,
                                      minGenes = minGenes)) %>%
  lapply(\(x) data.frame(num_module = x[1], num_contained_gene = x[2])) %>%
  list_rbind %>%
  mutate(resolution = seq(from = 0.25, to = 5, by = 0.25)) 

head(optimization_results)
## Plot optimization results ##
transform_factor =
  with(optimization_results,
       max(num_contained_gene)/max(num_module)
  )

optimization_results %>% 
  ggplot(aes(x = resolution)) +
  geom_line(aes(y = num_module),
            linewidth = 1.1, alpha = 0.8, color = "dodgerblue2") +
  geom_point(aes(y = num_module),
             size = 3, alpha = 0.7) +
  geom_line(aes(y = num_contained_gene/transform_factor),
            linewidth = 1.1, alpha = 0.8, color = "violetred2") +
  geom_point(aes(y = num_contained_gene/transform_factor),
             size = 3, alpha = 0.7) +
  scale_y_continuous(
    sec.axis = sec_axis(
      ~. * transform_factor,
      name = paste("Genes in modules")
    )) +
  labs(x = "Resolution",
       y = paste("Modules")) +
  theme_optimization +
  theme(axis.title = element_text(face="bold"))

ggsave("plots/MainAnalysis/Optimize_resolution.svg", 
       height = 5, width = 8, bg = "white")

#### Select resolution and detect modules ####
# Here highest resolution that still maintains all the genes
# from the node_table
# lose at max 0 genes
ngenes_threshold = max(optimization_results$num_contained_gene) - 0
resolutions = 
  (optimization_results %>%
  filter(num_contained_gene >= ngenes_threshold))$resolution %>%
  max
  
r_cutoff = edgeTableFiles[1] %>%
  basename() %>%
  gsub(pattern = ".+r(0.\\d+)_.+$", 
       replacement = "\\1")

for (curResolution in resolutions) {
  modules_ <- cluster_leiden(
    my_network, 
    resolution_parameter = curResolution, 
    objective_function = "modularity"
  )
  
  # Merge edge and node tables
  my_network_modules <- data.frame(
    gene_ID = names(membership(modules_)),
    module = as.vector(membership(modules_))) %>% 
    inner_join(node_table, by = "gene_ID")
  
  # Continue only with modules with 5 or more genes
  module_minGenes <- my_network_modules %>% 
    group_by(module) %>% 
    summarize(n=n()) %>% 
    arrange(-n) %>% 
    filter(n >= minGenes)
  
  my_network_modules <- my_network_modules %>% 
    filter(module %in% module_minGenes$module)
  
  ## Peak expression of modules
  Expr_averaged_z_high_var_modules <- 
    Exp_table_long_averaged_z_high_var %>% 
    inner_join(metadata[, c("tissue", "Part")],
               by = "tissue") %>%
    inner_join(my_network_modules, by = "gene_ID")
  
  modules_mean_z <- Expr_averaged_z_high_var_modules %>% 
    group_by(module, tissue, Part) %>% 
    summarise(mean.z = mean(z.score.TPM)) %>% 
    ungroup()
  
  # This will create a table with the module number and the Sample in which the module is most expressed.
  module_peak_exp <- modules_mean_z %>% 
    group_by(module, tissue, Part) %>% 
    slice_max(order_by = mean.z, n = 1)
  
  paste0("RDATA/MainAnalysis_Cor", r_cutoff) %>% dir.create()
  save(Expr_averaged_z_high_var_modules,
       module_peak_exp, modules_mean_z, my_network_modules,
       file = paste0(
         "RDATA/MainAnalysis_Cor", r_cutoff, "/ModuleData_forPlot_res",
         curResolution, ".RDATA")
  )
}
