#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]

library(dplyr)
library(data.table)

## Read file
gwas_data <- read.table(input_file, header = TRUE, stringsAsFactors = FALSE)

## format file
formatted_gwas <- gwas_data %>%
  dplyr::select(SNP, A1, A2, A1_FREQ, OR, SE, P, N) %>%
  dplyr::rename(
    freq = A1_FREQ,
    b = OR,
    se = SE,
    p = P,
    N = N
  )


formatted_gwas <- formatted_gwas %>%
  dplyr::mutate(b = log(b))

## save file
write.table(formatted_gwas, "formatted_gwas.ma", quote = FALSE, row.names = FALSE, sep = "\t")
