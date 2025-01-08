#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]

library(dplyr)
library(data.table)

# Read GWAS
gwas_data <- read.table(input_file, header = TRUE, stringsAsFactors = FALSE)

# Kolonları yeniden adlandırma ve seçme
formatted_gwas <- gwas_data %>%
  dplyr::select(SNP, A1, A2, A1_FREQ, OR, SE, P, N) %>%
  dplyr::rename(
    freq = A1_FREQ,
    b = OR,
    se = SE,
    p = P,
    N = N
  )

# Odds Ratio'yu log Odds Ratio'ya dönüştürme
formatted_gwas <- formatted_gwas %>%
  dplyr::mutate(b = log(b))

# Yeni formattaki dosyayı dışa aktarma
write.table(formatted_gwas, "formatted_gwas.ma", quote = FALSE, row.names = FALSE, sep = "\t")
