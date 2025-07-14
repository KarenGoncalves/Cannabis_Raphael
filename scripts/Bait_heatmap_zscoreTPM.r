library(tidyverse)
source("scripts/FUNCTIONS.R")
baits = read_delim("metadata/Baits_ensembl_ids.txt")
tpm_baits = read_delim("results/Exp_table_long_averaged_z_0.tsv")


baits_in_network = read_delim("results/Baits_in_network.tsv") %>% 
    left_join(baits, by = c("gene_ID", "Description")) %>%
    select(-starts_with("Entrez")) %>%
    arrange(module)
baits_in_network$Formatted.description =
    factor(baits_in_network$Formatted.description,
           levels = unique(baits_in_network$Formatted.description))

tpm_baits %>%
    mutate(SampleName = factor(SampleName,
                               levels = c("Meristem", "Trichome",
                                          "Flower", "Leaf", "Stem", "Root"))) %>% 
    inner_join(baits_in_network, by = "gene_ID") %>%
    ggplot(aes(SampleName, Formatted.description,
               fill = z.score.TPM)) +
    geom_tile() +
    heatmap_theme +
    scale_fill_gradient2(mid = "white",
                         high = "#67001F",
                         low = "#053061",
                         breaks = seq(-2,2,1),
                         labels = seq(-2,2,1)
    ) +
    labs(x = NULL, y = NULL, fill = "z-score TPM"
    )

ggsave("plots/Baits_heatmap.png",
       height = 5.86, width = 8,
       dpi=600)


ggsave("plots/Baits_heatmap.svg",
       height = 5.86, width = 8)


