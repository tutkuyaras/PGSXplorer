   #!/usr/bin/env Rscript
    library(ggplot2)
    library(data.table)
    args <- commandArgs(trailingOnly = TRUE)
    maf_freq <- read.table(args[1], header =TRUE, as.is=T)
    pdf("Target_MAF distribution.pdf")
    MAF <- hist(maf_freq[,5],main = "MAF distribution", xlab = "MAF",col="#d6604d")
    dev.off()
    