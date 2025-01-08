args <- commandArgs(trailingOnly = TRUE)
ibd_file <- args[1]

library(ggplot2)
library(data.table)
library(dplyr)

# Read data into R
ibd_data <- read.table(ibd_file, header = TRUE)

# Generate plot
pdf("ibd_plot.pdf")
hist(ibd_data[, "PI_HAT"], main = "Histogram of IBD (PI_HAT)", col = "#66c2a4")
dev.off()
