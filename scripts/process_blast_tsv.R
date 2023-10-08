suppressPackageStartupMessages({
  library(tidyverse)
})


args <- commandArgs(trailingOnly = TRUE) 
blast_tsv <- args[1]
align_minimum <- as.double(args[2])


column_names <- colnames(df)
df <- read_tsv(blast_tsv)
tryCatch({
  df <- df %>%
    filter(percent_aligned >= align_minimum) %>%
    group_by(query, accession) %>%
      mutate(max_read_pair_length = max(hit_to, hit_from) - min(hit_to, hit_from)) %>%
    ungroup() %>%
    arrange(hit_rank)
  df <- df[!duplicated(df[c("query", "blast_species")]),]
  write_tsv(df, blast_tsv)
}, error = function(err) {
  df <- data.frame(matrix(ncol = length(column_names), nrow = 0))
  colnames(df) <- column_names
  write_tsv(df, blast_tsv)
})
