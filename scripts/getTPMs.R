#!/usr/bin/env Rscript

library(GenomicRanges)
library(rtracklayer)

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript script.R input.gtf output.csv")
}

gtf_file <- args[1]
output_file <- args[2]

# Function to calculate gene length based on exon overlaps
calculate_gene_length <- function(gtf_file) {
  # Read GTF file
  gtf <- import(gtf_file)
  
  # Filter for exons
  exons <- gtf[gtf$type == "exon", ]
  
  # Split exons by gene
  exons_by_gene <- split(exons, exons$gene_id)
  
  # Calculate gene length by merging overlapping exons
  gene_lengths <- sapply(exons_by_gene, function(exons) {
    reduced_exons <- reduce(exons)
    sum(width(reduced_exons))
  })
  
  return(data.frame(Gene = names(gene_lengths), Length = gene_lengths))
}

# Run function and save results
gene_lengths_df <- calculate_gene_length(gtf_file)
write.csv(gene_lengths_df, output_file, row.names = FALSE)

cat("Gene lengths saved to:", output_file, "\n")
cat("Wrote lenths of", nrow(gene_lengths_df, "genes\n")
