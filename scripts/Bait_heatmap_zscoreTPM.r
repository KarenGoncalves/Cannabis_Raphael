library(tidyverse)
library(extrafont)
font_import()

setwd("D:/Work/Cannabis_Raphael/")
source("scripts/FUNCTIONS.R")
baits = read_delim("metadata/Baits_ensembl_ids2.txt") %>% 
    arrange(Order, Symbol) %>% 
    mutate(Symbol =
               factor(Symbol,
                      levels = unique(Symbol))
    )

tpm_baits = read_delim("results/Exp_table_long_averaged_z_0.tsv") %>%
    select(-tissue, -group, -mean.TPM,
           -mean.logTPM, -z.score.logTPM) %>%
    mutate(SampleName = 
               factor(SampleName,
                      levels = c("Meristem", "Trichome",
                                 "Flower", "Leaf", 
                                 "Stem", "Root")
                      )
           ) %>% 
    inner_join(baits, by = "gene_ID")


# baits_in_network = read_delim("results/Baits_in_network.tsv") %>% 
#     select(-Entrez, -description, -Description) %>%
#     left_join(baits, by = c("gene_ID")) %>%
#     select(-Entrez, -Formatted.description)
tpm_baits %>% pull(Symbol) %>% as.character %>% unique %>% length

tpm_baits %>% group_by(Symbol) %>% 
    slice_max(order_by = z.score.TPM, n = 1) %>% 
    select(Symbol, SampleName) %>% 
    ungroup %>% 
    unique %>% group_by(SampleName) %>% 
    count() 

tpm_baits  %>% 
    ggplot(aes(SampleName, Symbol %>% fct_rev(),
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
    ) + 
    theme(text = element_text(family = "Calibri"),
          legend.position = "bottom", 
          legend.key.height  = unit(.1, "in"),
          legend.ticks.length = unit(.02, "in"),
          axis.text.x = element_text(size=10),
          legend.text = element_text(size=10),
          legend.title.position = "top",
          legend.title = element_text(size=10))

ggsave("plots/Baits_heatmap_byName.png",
       height = 6.2, width = 4.4,
       dpi=600)

ggsave("plots/Baits_heatmap_byName.svg",
       height = 6.2, width = 4.4)

tpm_baits$Pathway %>% table
tpm_baits  %>% 
    ggplot(aes(SampleName, Symbol %>% fct_rev(),
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
    ) + 
    facet_grid(rows=vars(Pathway), scales="free_y", space="free_y" ) +
    theme(text = element_text(family = "Calibri"),
          legend.position = "bottom", 
          legend.key.height  = unit(.1, "in"),
          legend.ticks.length = unit(.02, "in"),
          axis.text.x = element_text(size=10),
          legend.text = element_text(size=10),
          legend.title.position = "top",
          legend.title = element_text(size=10))

ggsave("plots/Baits_heatmap_Pathway.svg",
       height = 6.2, width = 4.4)
ggsave("plots/Baits_heatmap_Pathway.png",
       height = 6.2, width = 4.4,
       dpi=1200)
