   #!/usr/bin/env Rscript
    library(data.table)
    args <- commandArgs(trailingOnly = TRUE)
    het <- read.table(args[1], header =TRUE, as.is=T)
    het_fail = subset(het, (het$HET_RATE < mean(het$HET_RATE)-3*sd(het$HET_RATE)) | (het$HET_RATE > mean(het$HET_RATE)+3*sd(het$HET_RATE)));
    het_fail$HET_DST = (het_fail$HET_RATE-mean(het$HET_RATE))/sd(het$HET_RATE);
    write.table(het_fail, "target_fail-het-qc.txt", row.names=FALSE)
    