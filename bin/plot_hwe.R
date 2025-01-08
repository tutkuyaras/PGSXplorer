args <- commandArgs(trailingOnly = TRUE)
hwe_file <- args[1]

library(ggplot2)
library(data.table)
library(dplyr)

# Read data into R
hwe <- read.table(hwe_file, header = TRUE)

# Generate plot
pdf("hwe.pdf")
hist(hwe[, 9], main = "Histogram HWE", col = "#f4a582")
dev.off()
