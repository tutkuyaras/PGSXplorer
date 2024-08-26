# Set CRAN mirror
options(repos = c(CRAN = "https://cran.rstudio.com"))
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

# Install necessary packages if not already installed
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("bigsnpr", quietly = TRUE)) remotes::install_github("privefl/bigsnpr") 
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table") 
if (!requireNamespace("magrittr", quietly = TRUE)) install.packages("magrittr")

library(bigsnpr)
library(dplyr)
library(R.utils)
library(data.table)
library(magrittr)
library(fmsb)

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
phenotype_file <- args[1]
eigenvec_file <- args[2]
gwas_sumstat <- args[3]
bed_file <- args[4]

# 1. Read in the phenotype and covariate files
phenotype <- fread(phenotype_file)
pcs <- fread(eigenvec_file)
colnames(pcs) <- c('FID', 'IID', paste0('PC', 1:10))
pheno <- merge(phenotype, pcs)
pheno$Phenotype <- ifelse(pheno$Phenotype == 2, 1, 0)

# 2. Information about the HapMap3+ variants
info <- readRDS(runonce::download_file(
  'https://figshare.com/ndownloader/files/37802721',
  dir = 'tmp-data2', fname = 'map_hm3_plus.rds'))
str(info)

# 3. Read external summary statistics
sumstats <- bigreadr::fread2(gwas_sumstat)
sumstats <- sumstats %>% rename('chr' = 'CHR', 'pos' = 'BP', 'a1' = 'A1', 'a0' = 'A2', 'A1_FREQ' = 'A1_FREQ', 'n_eff' = 'N', 'OR' = 'OR', 'beta_se' = 'SE', 'L95' = 'L95', 'U95' = 'U95', 'Z_STAT' = 'Z_STAT', 'p' = 'P', 'rsid' = 'SNP')
sumstats <- sumstats %>% select(chr, pos, rsid, a1, a0, n_eff, beta_se, p, OR)
sumstats$beta <- log(sumstats$OR)
sumstats <- sumstats[sumstats$rsid %in% info$rsid,]
fwrite(sumstats, 'ldpred2auto_output.txt')

# Calculate LD Matrix
NCORES <- nb_cores()
tmp <- tempfile(tmpdir = 'tmp-data2')
on.exit(file.remove(paste0(tmp, '.sbk')), add = TRUE)
corr <- NULL
ld <- NULL
fam.order <- NULL

# Debug messages
cat("Reading bed file...\n")
# Read from bed/bim/fam, it generates .bk and .rds files.
snp_readBed(bed_file)
cat("Finished reading bed file.\n")
# Debug message
cat("Attaching bigSNP object...\n")
# Attach the "bigSNP" object in R session
rds_file <- sub(".bed$", ".rds", bed_file)
cat("RDS file path:", rds_file, "\n")
obj.bigSNP <- snp_attach(rds_file)
# extract the SNP information from the genotype
map <- setNames(obj.bigSNP$map[-3], c('chr', 'rsid', 'pos', 'a1', 'a0'))
# perform SNP matching
df_beta <- snp_match(sumstats, map)

# Get aliases for useful slots
G <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
# Change affection column values into 0 and 1 
obj.bigSNP$fam$affection[obj.bigSNP$fam$affection == 1] <- 0
obj.bigSNP$fam$affection[obj.bigSNP$fam$affection == 2] <- 1
y <- obj.bigSNP$fam$affection

#correlation (take 3cm)
# To convert physical positions (in bp) to genetic positions (in cM), use
POS2 <- snp_asGeneticPos(CHR, POS, dir = "tmp-data2", ncores = NCORES)
# Before computing the LD matrices, let us filter out variants with small minor allele frequencies (MAFs):
ind.row <- rows_along(G)
maf <- snp_MAF(G, ind.row = ind.row, ind.col = df_beta$`_NUM_ID_`, ncores = NCORES)
maf_thr <- 1 / sqrt(length(ind.row))  # threshold I like to use
df_beta <- df_beta[maf > maf_thr, ]
#df_beta <- snp_match(sumstats, map, join_by_pos = FALSE)  # use rsid instead of pos

## Prepare data as validation and test
set.seed(1)
# Calculate the number of samples for validation data
n_val <- floor(0.7 * nrow(G))
# Sample indices for validation data
ind.val <- sample(nrow(G), n_val)
# Remaining indices for test data
ind.test <- setdiff(rows_along(G), ind.val)
# Debug message
cat("Split data as validation (70%) and test (30%)\n")

# Debug message
cat("Calculating LD matrix...\n")
# calculate LD
for (chr in 1:22) {

  # print(chr)
  # indices in 'df_beta'
  ind.chr <- which(df_beta$chr == chr)

  # indices in 'G'
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  
  # here we compute LD matrices ourselves, BUT
  # we recall that we provide pre-computed LD matrices that are 
  # usually much better (larger N, LD blocks for robustness, etc)
  corr0 <- snp_cor(G, ind.col = ind.chr2, size = 3 / 1000,
                   infos.pos = POS2[ind.chr2], ncores = NCORES)
  
  if (chr == 1) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp, compact = TRUE)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}

# We assume the fam order is the same across different chromosomes
fam.order <- as.data.table(obj.bigSNP$fam)
# Rename fam order
setnames(fam.order,c("family.ID", "sample.ID"), c("FID", "IID"))

# Debug message
cat("Perform LD score regression\n")

# Perform LD score regression
ldsc <- snp_ldsc(   ld, 
                    length(ld), 
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff, 
                    blocks = NULL)
h2_est <- ldsc[["h2"]]

# Calculate the null R2 (binary trait
# Reformat the phenotype file such that y is of the same order as the 
# sample ordering in the genotype file
a <- pheno[fam.order, on = c("FID", "IID")]
# Calculate the null R2
# use glm for binary trait 
# (will also need the fmsb package to calculate the pseudo R2)
null.model <- paste("PC", 1:10, sep = "", collapse = "+") %>%
  paste0("Phenotype~", .) %>%
  as.formula %>%
  glm(data = a, family=binomial) 
null.r2 <- fmsb::NagelkerkeR2(null.model)

### LD-Pred:auto
coef_shrink <- 0.95  # reduce this up to 0.4 if you have some (large) mismatch with the LD ref
set.seed(1)  # to get the same result every time
                          
 # Get adjusted beta from the auto model
 multi_auto <- snp_ldpred2_auto(
 corr,
 df_beta,
 h2_init = h2_est,
 vec_p_init = seq_log(1e-4, 0.9, length.out = 30),
 ncores = NCORES,  allow_jump_sign = FALSE, shrink_corr = coef_shrink)

 str(multi_auto, max.level = 1)

 str(multi_auto[[1]], max.level = 1)
# Debug message
cat("Visualization of auto model \n")

# Visualization of auto
library(ggplot2)
auto <- multi_auto[[1]]  # first chain

pdf("LDpred2 auto graph.pdf")
plot_grid(
qplot(y = auto$path_p_est) + 
theme_bigstatsr() + 
geom_hline(yintercept = auto$p_est, col = "blue") +
scale_y_log10() +
labs(y = "p"),
qplot(y = auto$path_h2_est) + 
theme_bigstatsr() + 
geom_hline(yintercept = auto$h2_est, col = "blue") +
labs(y = "h2"),
ncol = 1, align = "hv")
dev.off()


(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
(keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE))))

# Debug message
cat("Obtain Model PRS \n")

# (Obtain Model PRS)
beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
pred_auto <- big_prodVec(G, beta_auto, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
pcor(pred_auto, y[ind.test], NULL)

# Get final effects(for all) 

# Get final effects (for all data)
beta_auto_all <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
pred_auto_all <- big_prodVec(G, beta_auto, ind.row = ind.row, ind.col = df_beta[["_NUM_ID_"]])

reg.formula <- paste("PC", 1:10, sep = "", collapse = "+") %>%
paste0("Phenotype~", .) %>%
as.formula
reg.dat <- a
reg.dat$PRS <- pred_auto_all
auto.model <- lm(reg.formula, dat=reg.dat) %>%
summary
(result <- data.table(
  auto = auto.model$r.squared - null.r2$R2,
  null = null.r2$R2
))

write.table(reg.dat, "LDpred2_auto.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

# Debug message
cat("Auto results was written as LDpred2_auto.txt file \n")

# Debug message
cat("Some cleaning \n")
# Some cleaning
rm(corr); gc(); file.remove(paste0(tmp, ".sbk"))
