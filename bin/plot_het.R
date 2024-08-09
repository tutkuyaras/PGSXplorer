args <- commandArgs(trailingOnly = TRUE)
het_file <- args[1]

library(ggplot2)
library(data.table)
library(dplyr)

# Read data into R
het <- read.table(het_file, header = TRUE)
het$HET_RATE = (het$"N.NM." - het$"O.HOM.") / het$"N.NM."

# Generate plot
pdf("het.pdf")
hist(het$HET_RATE, xlab = "Heterozygosity Rate", ylab = "Frequency", main = "Heterozygosity Rate", col = "#92c5de")
dev.off()
