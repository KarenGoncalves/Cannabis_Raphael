suppressMessages(library(tidyverse))

args <- commandArgs(trailing.only = T)
inputFolder = args[1]
folderfiles <- list.files(inputFolder, pattern="*ReadsPerGene.out.tab", full.names=T)

outputFile = args[2]

data_csv <- folderfiles %>%
    set_names() %>%
    map_dfr(.f = read_delim,
            delim = "\t",
                        col_names=c("Gene_ID", "Unstranded", "First_stranded", "Second_stranded"),
            .id = "file_name") %>%
        mutate(Run = gsub("_ReadsPerGene.out.tab", "", basename(file_name))
        )

metadata = read_delim("metadata/Read_lengths.txt", col_names = c("Run", "Layout_type", "Read_length"), skip = 1)
stranded_info = read_delim("metadata/Stranded_data.txt", col_names = c("Run", "Layout_type", "Strand"))
data_csv$Gene_ID = data_csv$Gene_ID %>%
        factor(levels = unique(data_csv$Gene_ID)
)

counts = sapply(unique(data_csv$Run), simplify = F, \(run) {
        count_col = stranded_info$Strand[stranded_info$Run == run]
        (data_csv %>% filter(Run == run) %>%
        arrange(Gene_ID))[[count_col]]
        }) %>% as.data.frame %>%
        mutate(Gene_ID = levels(data_csv$Gene_ID),
        .before = everything())

write_delim(counts, delim = "\t", file=outputFile)
