#!/usr/bin/env Rscript

for (i in c("GenomicRanges", "rtracklayer", "tidyverse") ) {
 suppressMessages(library(i, character.only=T))
}

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if ((args[1] == "lengths" & length(args) != 3) ||
    (toupper(args[1]) == "TPM" & length(args) != 4)) {
  stop("Usage: 
		Rscript lengths_TPM.R lengths input.gtf output.csv
	OR
		Rscript lengths_TPM.R TPM lengths.csv readCounts.csv TPM.csv
  ")
}


# Function to calculate gene length based on exon overlaps
calculate_gene_length <- function(gtf_file) {
  # Read GTF file
  gtf <- import(gtf_file)
  
  # Filter for exons and split them by gene
  exons_by_gene <- gtf[gtf$type == "exon", ] %>% split(.$gene_id)
  
  # Split exons by gene
#  exons_by_gene <- split(exons, exons$gene_id)
  
  # Calculate gene length by merging overlapping exons
  gene_lengths <- sapply(exons_by_gene, function(exons) {
    reduce(exons) %>% width %>% sum
  })
  
  return(data.frame(Gene = names(gene_lengths), Length = gene_lengths))
}

Counts_to_tpm <- function(counts, featureLength) {
  # Find matching row names
  common_names <- intersect(rownames(counts), names(featureLength))
  
  # Raise a warning if some features are missing
  if (length(common_names) < nrow(counts)) {
    warning("Some features in 'counts' do not have matching lengths in 'featureLength'. ",
            paste0("Only processing matching features: ", length(common_names), ".")
	)
  }

  # Subset data to only include matching features
  counts_ <- counts[common_names, , drop = FALSE]
  effLen <- featureLength[common_names]

  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts_), function(i) {
    rate = log(counts_[,i]) - log(effLen)
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))

  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts_)
  return(
    tpm %>% as.data.frame %>%
      mutate(Gene_ID = rownames(counts_), .before = everything()
    )
  )
}

# Run function and save results
if ( args[1] == "length" ) {
	# Read inputs
	gtf_file <- args[2]
	output_file <- args[3]

	# Run function
	gene_lengths_df <- calculate_gene_length(gtf_file)

	# Prepare output
	write.csv(gene_lengths_df, output_file, row.names = FALSE)
	cat("Gene lengths saved to:", output_file, "\n")
	cat("Wrote lenths of", nrow(gene_lengths_df), "genes\n")
	# Done with lengths
} else if ( toupper(args[1]) == "TPM") {
	# Read inputs
	## Prepare featureLength
	lengths_table <- args[2] %>% read.csv(header=T)
	featureLength <- lengths_table[[2]]
	names(featureLength) <- lengths_table[[1]]
	## Prepare counts
	reads <- args[3] %>% read_delim(delim="\t")
	counts <- reads %>% dplyr::select(-Gene_ID) %>%
			as.data.frame
	rownames(counts) <- reads$Gene_ID

	# Run function and output result
	Counts_to_tpm(counts, featureLength) %>% 
	write_delim(append = F, quote = "none", file = args[4])
}
