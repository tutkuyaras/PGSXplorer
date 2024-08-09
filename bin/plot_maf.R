library(ggplot2)
library(data.table)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
maf_freq_file <- args[1]
# Read data into R
maf_freq <- read.table("MAF_check.frq", header =TRUE, as.is=T)

# Generate plot
pdf("MAF_dist.pdf")
MAF <- hist(maf_freq[,5], main = "MAF distribution", xlab = "MAF", col="#d6604d")
dev.off()
