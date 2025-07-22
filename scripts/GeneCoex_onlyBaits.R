source("scripts/FUNCTIONS.R") # loads packages too

Exp_table_long_averaged_z = "results/Exp_table_long_averaged_z_0.tsv" %>% 
    read_delim
Baits = read_delim("metadata/Baits_ensembl_ids.txt") %>%
    filter(gene_ID %in% Exp_table_long_averaged_z$gene_ID)
metadata = read_delim("metadata/metadata_pca.txt")

# Order tissues the way they should appear in plot
tissue_order <- c("Trichome", "Flower", "Meristem", "Leaf", "Petiole", "Stem", "Root")

metadata <- metadata %>%
    mutate(TissueOrdered = factor(tissue, levels = tissue_order))

z_score_wide <- Exp_table_long_averaged_z %>%
    mutate(tissue = SampleName) %>%
    select(gene_ID, tissue, z.score.TPM) %>% 
    pivot_wider(names_from = tissue, values_from = z.score.TPM) %>% 
    as.data.frame()
nreps = ncol(z_score_wide) - 1

correlation_against_baits <-
    sapply(Baits$gene_ID, simplify = F, \(x) {
        dBait = (z_score_wide %>% filter(gene_ID == x))[-1] %>%
            as.matrix %>% t
        dOthers = (z_score_wide %>% filter(gene_ID != x))[-1] %>%
            as.matrix %>% t
        data.frame(correlation = cor(dBait, dOthers) %>% t) %>% 
            mutate(Gene1 = x,
                   Gene2 = (z_score_wide %>% filter(gene_ID != x))$gene_ID,
                   r2 = correlation^2) %>%     
            filter(Gene1 != Gene2) %>% 
            mutate(t = correlation*sqrt( (nreps-2) / (1-r2) ) ) %>% 
            mutate(p.value = case_when(
                t > 0 ~ pt(t, df = nreps-2, lower.tail = F),
                t <=0 ~ pt(t, df = nreps-2, lower.tail = T)
            )) %>% 
            mutate(FDR = p.adjust(p.value, method = "fdr"),
                   significant = ifelse(FDR < 0.01, T, F)) 
    }) %>% list_rbind



correlation_against_baits %>% 
    ggplot() +
    geom_histogram(aes(x = correlation),
                   bins = 100) +
    geom_vline(xintercept = .7, 
               color = "black",
               linewidth = 0.85) +
    # check at which r value, the number of correlations starts to drop fast
    labs(title = "Gene correlations - Bait genes") +
    scale_x_continuous(
        breaks = seq(-1, 1, .4)) +
    theme_classic() +
    theme(text = element_text(size = 14),
          axis.text = element_text(color = "black")
    )

r_cutoff = .7

edge_table = correlation_against_baits %>%
    filter(correlation > r_cutoff) %>%
    rename(from = Gene1,
           to = Gene2) %>% 
    relocate(from, to, .before = everything()) 

##### Detect modules #####
## Minimum number of genes to form a module ##
minGenes = 3
# future::plan() - sets parallel session (using multisession())
future::plan(multisession, workers = 1)

#### Find gene co-expression modules ####
# with Leiden algorithm (graph-based method)
# Load gene annotation table
funct_anno = read_delim("results/Annotation_expressed_genes.txt")

# Merge annotation and edge tables
node_table <- data.frame(
    gene_ID = c(edge_table$from, 
                edge_table$to) %>% unique()
) %>% 
    left_join(funct_anno[, 1:2], # node table cannot have duplicated transcript ids
              # so we keep only the transcript id and description, not the GO annotation
              by = join_by("gene_ID" == "ensembl_gene_id")) %>%
    unique()

head(node_table)
dim(node_table)

# Create graph
my_network <- graph_from_data_frame(
    edge_table,
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

# ggsave("plots/MainAnalysis/Optimize_resolution.svg", 
#        height = 5, width = 8, bg = "white")

#### Select resolution and detect modules ####
# Here highest resolution that still maintains all the genes
# from the node_table
# lose at max 0 genes
ngenes_threshold = max(optimization_results$num_contained_gene) - 0
resolution = 
    (optimization_results %>%
         filter(num_contained_gene >= ngenes_threshold))$resolution %>%
    max



modules_ <- cluster_leiden(
    my_network, 
    resolution_parameter = resolution, 
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
    Exp_table_long_averaged_z %>% 
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


#### Heatmap #####
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
         caption = paste0("r threshold = ", r_cutoff, 
                          "\nresolution = ", resolution)
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
ggsave("plots/GeneCoex_onlyBaits_r07_res175.svg")
ggsave("plots/GeneCoex_onlyBaits_r07_res175.png")

left_join(Baits, my_network_modules, by = "gene_ID") %>%
    arrange(module) %>%
    select(-Entrez, -Description, -description) %>%
    write_delim("results/GeneCoex_onlyBaits_moduleBaits.txt", 
                delim = "\t")

Expr_averaged_z_high_var_modules %>%
    select(gene_ID, tissue, module, description) %>%
    left_join(Baits[, c(1,4)], by = "gene_ID") %>%
    write_delim("results/GeneCoex_onlyBaits_geneAnnotations.txt", 
                delim = "\t")
