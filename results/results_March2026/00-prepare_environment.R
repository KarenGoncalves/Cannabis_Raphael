#!/usr/bin/env Rscript

# Prepare environment by creating right folders
options(repos=c("https://packagemanager.rstudio.com/all/__linux__/focal/latest",
                "https://cloud.r-project.org/"))

#folders = c("plots", "counts", "results", "RDATA")

#for (folder in folders) dir.create(folder)

# then download the necessary packages 

pkgs=c("devtools", "tidyverse", "igraph", "ggraph",
       "readxl", "patchwork", "RColorBrewer",
       "viridis", "parallel", "doParallel", "scriptName",
       "furrr", "future", "svglite", "BiocManager"
)
install.packages(pkgs[!pkgs %in% rownames(installed.packages())])
 
sapply(pkgs, \(x)  library(x, character.only=T))

bioconductor=c("DESeq2", "GenomicFeatures", "biomaRt", "GenomicRanges", "rtracklayer")

BiocManager::install(bioconductor[bioconductor %in% rownames(installed.packages())])
sapply(bioconductor, \(x)  library(x, character.only=T))

# CustomSelection has to be installed from github, so it cannot be grouped with the others yet
if (!require("CustomSelection") ) {
  devtools::install_github("KarenGoncalves/CustomSelection")
  library(x, character.only=T)
}

