# PGSXplorer
PGSXplorer is a bioinformatics workflow designed to calculate polygenic scores by processing genomic data through quality control steps, phasing and imputation. Optionally, it can utilize tools such as [PLINK](https://www.cog-genomics.org/plink/), [PRSice-2](https://choishingwan.github.io/PRSice/), [LD-Pred2 (grid)](https://privefl.github.io/bigsnpr/articles/LDpred2.html), [LD-Pred2 (auto)](https://privefl.github.io/bigsnpr/articles/LDpred2.html), [Lassosum2](https://privefl.github.io/bigsnpr/articles/LDpred2.html#lassosum2-grid-of-models), [MegaPRS](https://dougspeed.com/megaprs/), [PRS-CSx](https://github.com/getian107/PRScsx), and [MUSSEL](https://github.com/Jin93/MUSSEL). The workflow requires genomic files in PLINK format (.bed, .bim, .fam) and GWAS summary statistics for two different populations as input to complete the analysis.

You can try PGSXplorer with the sample data set in the [Target](https://drive.google.com/drive/folders/1u6iAEZaDpq9U-EfrRNv2fZLEW-WXBPM2?usp=drive_link) folder and compare the output with your own data. We also added examples of the necessary reference files for chromosomes 1 and 2 to this folder, which will allow users to experiment with sample data sets. This dataset is generated using [HAPNEST](https://github.com/intervene-EU-H2020/synthetic_data) with default parameters and is for a European population of 1000 people.

![PGSXplorer Diagram](https://github.com/tutkuyaras/PGSXplorer/blob/PGSXplorer/images/PGSExplorer%20Workflow.drawio.png)
## Workflow Overview

PGSXplorer includes a comprehensive pipeline that begins with rigorous quality control (QC) measures to ensure the integrity of genomic data. Following the completion of the QC module, users can optionally execute various polygenic score (PGS) calculation tools. The steps are as follows:

## Workflow Overview

PGSXplorer includes a comprehensive pipeline that begins with rigorous quality control (QC) measures to ensure the integrity of genomic data. Following the completion of the QC module, users can optionally execute various polygenic score (PGS) calculation tools. The steps are as follows:  
**1. GWAS QC:** Applying standard GWAS QC steps for GWAS summary statistic file.    
**2. Filtering Missing SNPs:** Identify and remove SNPs with missing genotype data.    
**3. Filtering Missing Individuals:** Exclude individuals with excessive missing genotype data.    
**4. Filtering by Minor Allele Frequency (MAF):** Retain SNPs above a certain MAF threshold.    
**5. Visualization of MAF Distributions:** Graphical representation of MAF across SNPs.  
**6. Filtering by Hardy-Weinberg Equilibrium (HWE):** Remove SNPs that deviate significantly from HWE.  
**7. Visualization of HWE Distributions:** Visual display of HWE p-values for SNPs.  
**8. Relatedness Checking:** Assess genetic relatedness between individuals.  
**9. Visualization of Identity by Descent (IBD):** Visualize IBD statistics to identify related individuals.  
**10. Heterozygosity Assessment:** Evaluate heterozygosity rates to identify potential outliers.  
**11. Visualization of Heterozygosity Distributions:** Plot distribution of heterozygosity rates.  
**11. Removal of Duplicate SNPs:** Eliminate duplicated SNPs to prevent redundancy.  
**12. Convert to VCF File:** Prepare genomic files for Phasing and Imputation  
**13. Phasing:** Genomic Data Phasing by using Eagle  
**14. Imputation:** Genomic Data Imputation by using Beagle  
**15. Post-Imputation QC:** Filtering VCF file by imputation score  
**16. Convert files into PLINK format:** Prepare files for rest of the pipeline  
**17. fastmixture:** Analyzing data for target ancestry inference   
**18. Calculating 10 Principal Components (PCA):** Perform PCA to capture population structure.  
**19. Pruning and Thresholding (P+T):** Execute P+T method using PLINK to calculate PGS.   
**20. PRSice-2:** Calculate and visualize PGS using PRSice-2.  
**21. LD-Pred2-grid:** Apply LD-Pred2 grid model for PGS estimation.  
**22. LD-Pred2-auto:** Apply LD-Pred2 auto model for PGS estimation.   
**23.Lassosum2:** Apply Lassosum2 for PGS estimation.     
**24.MegaPRS:** Apply MegaPRS for PGS estimation.   
**25. PRS-CSx:**  Implement PRS-CSx for multi-ancestry PGS.    
**26. MUSSEL:** Utilize MUSSEL for multi-ancestry PGS estimation.    

## Computational Requirements
Our tool is designed to work seamlessly on Linux-based local machines, as well as in HPC and cloud environments. To run the workflow, users need to have Nextflow and Docker installed on their systems.
For those who prefer not to use Docker, the required tools must be pre-installed on their system for the corresponding functions. Please refer to the table below for the list of tools needed for each step of the workflow. Please ensure that the versions of the tools installed are compatible with the workflow requirements. Detailed installation instructions for each tool can be found in their respective official documentation.  

```
|STEP         | FUNTIONALITY                       | TOOLS USED                          |  
|-------------|------------------------------------|-------------------------------------|  
|QC           | SNP and sample filtering           | PLINK, R, Shell                     |  
|-------------|------------------------------------|-------------------------------------|  
|C+T          | LD Pruning and Thresholding        | PLINK                               |  
|-------------|------------------------------------|-------------------------------------|  
|Phasing      | Phasing Genotype Data              | Eagle                               |  
|-------------|------------------------------------|-------------------------------------|  
|Imputation   | Genotype Imputation                | Beagle                              |
|-------------|------------------------------------|-------------------------------------|  
|PGS          | Single Ancestry                    | LDPred2, PRSice2,Lassosum 2, MegaPRS|
|-------------|------------------------------------|-------------------------------------| 
|PGS          | Multi Ancestry                     | MUSSEL, PRS-CSx                     |
|-------------|------------------------------------|-------------------------------------| 
|Visualization| Distribution and performance plots | R (custom scripts)                  |
|-------------|------------------------------------|-------------------------------------|
```

## Usage
This repository hosts a Dockerized version of the PGSXplorer pipeline, making it easy to run the entire analysis environment in a consistent and repeatable manner. With Docker, users can quickly get the pipeline up and running without worrying about software dependencies or compatibility issues. To get started, make sure Docker is installed on your system, pull the Docker image, and then run the Nextflow pipeline using the Docker profile. The Docker image includes all the necessary tools and configurations for smooth execution.  

First of all, clone the repository

```
git clone https://github.com/tutkuyaras/PGSXplorer.git
cd /PGSXplorer
```
Then, pull the docker image

```
docker pull tutkuyaras/pgsxplorer_image:v2
```
Then, you can run basically with following command:

```
nextflow run main.nf -profile docker
```

In order for PGSXplorer to work correctly, you must first examine the parameters of the quality control steps and the polygenic score calculation tools you want to run and start the analysis with the appropriate command. You can also access the parameters given below by running the command 

```
nextflow run main.nf --help
```

### Parameters
```
    // Help Flag
    --help = false

    Required Arguments:
    --vcf_to_plink = false
    --target = "$PWD/target/"
    --target_prefix = "$PWD/target/target"
    --target_qc_prefix = "$PWD/outputs/target_9"
    --pheno_file = "$PWD/target/phenotype_file_t2.txt"  // Default phenotype file
    --gwas_sumstat = "$PWD/target/GWAS_sumstat_t1.txt"  // Default GWAS summary statistics file
    --mega_summaries = "$PWD/outputs/quant.summaries"
    
    Optional Arguments: 
    // Quality Control Parameters
    --mind                  Threshold value to delete SNPs with missingness (Default: 0.02)
    --geno                  Threshold value to delete individuals with missingness (Default: 0.02)
    --maf                   Threshold value for minimum MAF frequencies of SNPs (Default: 0.05)
    --hwe_case              Threshold value for Hardy-Weinberg Equilibrium for controls (Default: 1e-6)
    --hwe_ctrl              Threshold value for Hardy-Weinberg Equilibrium for cases (Default: 1e-10)
    --indep_window_size     The window size (Default: 100)
    --indep_step_size       The number of SNPs to shift the window at each step (Default: 5)
    --indep_threshold       The multiple correlation coefficient for a SNP being regressed on all other SNPs simultaneously. (Default: 0.2)
    --pihat                 The default threshold 0.1875 represents the half-way point between 2nd and 3rd degree relatives (Default: 0.1875)
    --relatedness           The same threshold with pihat value (Default: 0.1875)
    --pca                   Number of principal components to compute (Default: 10 )
    --pheno_file            Name of the phenotype file located under the target folder

    // Phasing - Imputation
    ref_phase = "$PWD/1KG_hg38_BCF/"
    genetic_map = "$PWD/genetic_map_hg38_withX.txt.gz"
    ref_imp = "$PWD/1KG_hg38_bref3/" 
    g_map = "$PWD/plink.GRCh38.map/" // Genetic map directory

   // PGS Parameters
    --run_plink             Run the PLINK part of the workflow if set to true
    --run_prsice            Run the PRSice-2 part of the workflow if set to true
    --run_pca               Run the PCA part of  the workflow if set to true, cretaes ".eigenvec" file
    --run_LDpred2grid       Run the LDPred2 Grid Model of the workflow if set to true
    --run_LDpred2auto       Run the LDPred2 Auto Model of the workflow if set to true
    --run_Lassosum2         Run the Lassosum2 model of the workflow, default = true
    --run_megaprs           Run the MegaPRS model of the workflow, default = true
    --run_prscsx            Run the PRScsx of the workflow, default = true
    --run_mussel            Run the MUSSEL of the workflow if set to true, default = false

    // MegaPRS Parameters
    ldak_executable = "$PWD/bin/ldak6.mac" //Path to the LDAK executable (e.g., ldak6.mac or ldak6.linux)
    mega_model = "bayesr"  // Model type for MegaPRS. Options: lasso, ridge, bayesr, etc. Default = bayesr

    // PRS-Csx parameters
    --prsice_script         Path to the PRSice R script
    --prsice_executable     Path to the PRSice executable
    --prscsxref_dir         Path to PRScsx reference panel, could be UKBB or 1KG
    --prscsx_gwas1          1st GWAS sum stat for PRScsx analysis
    --prscsx_gwas2          2nd GWAS sum stat for PRScsx analysis
    --n_gwas1               Sample size for GWAS1
    --n_gwas2               Sample size for GWAS2
    --pop1                  Ancestry of 1st Gwas sum stat, could be AFR,AMR,EAS,EUR,SAS
    --pop2                  Ancestry of 2nd Gwas sum stat, could be AFR,AMR,EAS,EUR,SAS
    --phi                   Global shrinkage parameter phi, fixing phi to 1e-2(for highly polygenic traits) or 1e-4(for less polygenic tratits)
    --meta                  Return combined SNP effect sizes across populations using  inverse variance weighted meta-analysis of population-specific posterior effect size estimates. Default is True.
    --out_name              Specify output prefix, e.g 'sample_ukbb_meta'

    // MUSSEL parameters
    --pack                  Path to MUSSEL folder
    --data                  Path to MUSSEL data
    --LDref                 Path to LDref folder for MUSSEL module
    --sst                   Path to summary statistic data files for MUSSEL module
    --pop                   Used populations for MUSSEL module, could be EUR, AFR, AMR, EAS or SAS
    --mussel_chrom          Specify the chromosomes to be analyzed by MUSSEL module , by default all chromosomes are analyzed
    --bfile_tuning          Path to PLINK binary input file prefix for tuning of MUSSEL module
    --pheno_tuning          Path to phenotype file (PLINK format) for tuning of MUSSEL module
    --covar_tuning          Path to quantitative covariates (PLINK format) for tuning
    --bfile_testing         Path to PLINK binary input file prefix for testing of MUSSEL module
    --pheno_testing         Path to phenotype file (PLINK format) for testing of MUSSEL module
    --covar_testing         Path to quantitative covariates (PLINK format) for testing
    --trait_type            Type of phenotype, continuous or binary for MUSSEL module. Default: continuous
    --NCORES                How many cores to use for MUSSEL modules
    --plink                 path to plink2 for MUSSEL module

    // Output directories
    outdir = "$PWD/outputs"
    graphs = "$PWD/QC_graphs"
    
    
```
## Running Examples
To execute the workflow exclusively for Quality Control (QC), phasing, and imputation, use the following command when the target data is in PLINK format (e.g., .bed, .bim, .fam) with a prefix named target:  

```
nextflow run main.nf --run_plink false --run_prsice false --run_LDpred2grid false --run_LDpred2auto false --run_Lassosum2 false --run_prscsx false --run_megaprs false
```
If the input data is in VCF format, use the following command to convert the VCF file into PLINK format before proceeding with QC, phasing, and imputation:    
```
nextflow run main.nf --run_prsice false --run_LDpred2grid false --run_LDpred2auto false --run_Lassosum2 false --run_prscsx false --vcf_to_plink true --vcf ./target/target.vcf
```
For files with different names or to provide custom paths for reference files, use the --help flag to list all available parameters and their descriptions:  
```
nextflow run main.nf --help
```
This flexibility allows users to adjust parameters as needed, ensuring compatibility with different file names, formats, and reference paths.  

## Quality Control 

Quality control modules consist of seven basic steps which are not optional. It consists of SNP filtering, individual filtering, filtering by MAF, HWE, relatedness, heterozygosity and elimination of duplicate SNPs.   
The files obtained after each step are saved in the **outputs** folder. In addition, the distribution graphs of MAF, HWE, heterozygosity and kinship relationships are also saved in the **QC_graphs** folder.
  
Only for quality control module, you can run the pipeline using: 

```
nextflow run main.nf --run_prsice false --run_LDpred2grid false --run_LDpred2auto false --run_Lassosum2 false --run_prscsx false --run_mussel false
```
## Phasing 
Involves determining the arrangement of alleles on homologous chromosomes. This process is crucial for identifying haplotypes, which are essential for downstream analyses like imputation. Phasing reconstructs the sequence of genetic variants by aligning them to a reference panel, ensuring accurate representation of genetic variation.  
For the phasing step, the reference files were downloaded from the [1000 Genomes Project FTP site](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/). The reference data includes biallelic SNVs and INDELs, which are critical for accurate phasing.    
The downloaded .vcf.gz files were processed as follows:  
1. Conversion to BCF format:  
Using bcftools view:  
```
bcftools view -O b -o output_file.bcf input_file.vcf.gz
```
2. Indexing:
The converted BCF files were indexed using:

```
bcftools index output_file.bcf
```
The phased haplotypes were generated using Eagle, a state-of-the-art phasing algorithm optimized for speed and accuracy. Eagle enables efficient large-scale phasing and is particularly effective when paired with high-quality reference panels. For detailed information about Eagle v2.4.1, refer to the official documentation at [Broad Institute's Eagle page](https://alkesgroup.broadinstitute.org/Eagle/).  

## Imputation 
Predicts missing genotypes in genomic datasets by leveraging a reference panel of well-characterized genetic data. It enhances data completeness and resolution, allowing for more comprehensive genome-wide association studies (GWAS) and increasing the power to detect genetic associations. By filling in missing data points, imputation improves the accuracy of genetic analyses, especially in studies with sparse or incomplete datasets.  
For the imputation step, the reference files were downloaded from the [1000 Genomes Project FTP site](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/). These files are converted to bref3 files. For detailed information about Beagle v5.4, refer to the official documentation at [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html).  


## Single Ancestry PGS Tools
Default parameters were used for the integration of tools that calculate polygenic scores using summary statistics for a single population.
PGS calculations for PLINK, PRSice-2, LD-Pred2grid, LD-Pred2-auto and Lassosum2 were defined as default true.  
As input;  
- _Phenotype file_ (in binary format) should be given with the --pheno_file parameter, the file format should be as in the example


```
FID	IID	Phenotype
syn1	syn1	1
syn2	syn2	1
syn3	syn3	1
syn4	syn4	1
syn5	syn5	2
syn6	syn6	2
syn7	syn7	1
syn8	syn8	1
syn9	syn9	2
syn10	syn10	2
syn11	syn11	1
syn12	syn12	2
syn13	syn13	1
syn14	syn14	1
syn15	syn15	2
syn16	syn16	2
```

  
- _GWAS summary statistics_ should also be given with the --gwas_sumstat parameter, the file format should be as in the example

```
CHR BP 	   SNP     A2 A1 A1_FREQ N     OR       SE 	L95 	U95 	Z_STAT      P   INFO
1 793571 rs11240767 C T 0.19979 9978 1.08025 0.0356268 1.00739 1.15838 2.16671 0.0302568 1
1 817341 rs3131972 G A 0.468882 9978 0.973613 0.028382 0.920932 1.02931 -0.942194 0.346094 1
1 818802 rs3131969 G A 0.393065 9978 1.01116 0.0291602 0.954985 1.07063 0.380434 0.703623 1
1 818954 rs3131967 C T 0.131038 9978 1.0179 0.0418892 0.937665 1.10499 0.423476 0.671948 1
1 825532 rs1048488 T C 0.230908 9978 1.03046 0.0336042 0.964776 1.10061 0.892852 0.371936 1
1 833068 rs12562034 G A 0.0872419 9978 1.02131 0.0501883 0.925629 1.12688 0.420125 0.674394 1
1 841166 rs12124819 A G 0.467529 9978 1.00723 0.0283809 0.952734 1.06485 0.253914 0.799562 1
1 903175 rs4970383 C A 0.223642 9978 0.972941 0.0339555 0.910298 1.0399 -0.807874 0.419163 1
1 917392 rs28587382 G A 0.333834 9978 0.965077 0.0299153 0.910119 1.02335 -1.18827 0.234729 1
1 918574 rs1806509 A C 0.386901 9978 1.04957 0.0292797 0.991033 1.11156 1.65235 0.0984641 1
1 918870 rs7537756 A G 0.0792744 9978 1.08592 0.0521271 0.980451 1.20273 1.58122 0.113827 1
1 922448 rs6694982 A G 0.350972 9978 1.0003 0.0295217 0.944068 1.05989 0.0103257 0.991761 1
1 926428 rs13302982 G A 0.469683 9978 1.02356 0.0284034 0.968133 1.08215 0.819746 0.412361 1
1 927744 rs4040604 T G 0.337693 9978 0.994669 0.0300084 0.937854 1.05492 -0.17814 0.858613 1
```

> **NOTE:** For PGSXplorer to work correctly, file formats, column order and names must be the same.

### PCA Calculation
In order to use PRSice-2, LD-Pred2 grid, LD-Pred2 auto, LD-Pred2 auto, Lassosum2 PGS models, the .eigenvec file used as covariate is needed. This file is created with the pca module in PGSExplorer. The default pca value is 10. You can change this by using the **--pca** parameter.

### PLINK
PLINK tool was used to calculate PGS with Pruning and Thresholding method and [this](https://choishingwan.github.io/PRS-Tutorial/plink/) tutorial was followed. 

This section consists of the following modules.  
>- ClumpSNPs  
>- CreateValidSNPs  
>- CreateSNPpvalue  
>- CreateRangeList  
>- PruneSNPs  
>- CalculatePGS   

For pruning in PLINK module, indep_window_size = 100, indep_step_size = 5, indep_threshold = 0.2 are defined by default. For example, you can change one or all of them as follows.

```
nextflow run main.nf --indep_window_size 200  
```
 
### PRSice-2
PRSice-2 is a comprehensive software tool for calculating polygenic scores (PGS) by integrating genome-wide association study (GWAS) summary statistics with individual genotype data. It provides flexible options for clustering and thresholding to optimize PGS structure, enabling analysis of genetic susceptibility to various traits and diseases. You can check details from [here](https://github.com/choishingwan/PRSice)

### LD-Pred2(Grid) & LD-Pred2(Auto) and Lassosum2 
LDpred2 is a Bayesian polygenic score (PGS) tool that includes two primary models: LDpred2-grid and LDpred2-auto. LDpred2-grid explores a grid of hyperparameters to find the best-fit PGS model, while LDpred2-auto automatically adjusts parameters based on the data, eliminating the need for predefined hyperparameters. Lassosum2 is another PGS method that applies penalized regression techniques, effectively handling linkage disequilibrium (LD) patterns to improve the prediction accuracy of polygenic scores. You can check details from [here](https://privefl.github.io/bigsnpr/articles/LDpred2.html).  

### MegaPRS
MegaPRS is a powerful tool for calculating polygenic risk scores (PRS) by integrating GWAS summary statistics with individual-level genotype data. It uses advanced modeling approaches to account for linkage disequilibrium and population structure, improving PRS prediction accuracy across diverse ancestries. For more details, visit the [MegaPRS](https://dougspeed.com/megaprs/)

## Multi Ancestry PGS Tools  
In this section we integrate two tools that improve polygenic prediction power using GWAS data from multiple populations. 

### PRS-Csx
PRS-CSx is an advanced polygenic score (PGS) tool that extends [PRS-CS](https://github.com/getian107/PRScs) to accommodate trans-ancestry genetic data. By leveraging summary statistics from multiple populations, PRS-CSx improves the accuracy and generalizability of polygenic scores across diverse ancestries, making it particularly useful for studies involving multi-ethnic cohorts. You can check details from [here](https://github.com/getian107/PRScsx).

* For this tool, GWAS summary statistics for 2 different populations should be used. These populations could be AFR, AMR, EAS, EUR or SAS. These populations should be determined with **--prscsx_gwas1** and **--prscsx_gwas2** parameters. Populations should be determined with **pop1** and **pop2** parameters. Lastly populations' size should be determined with **n_gwas1** and **n_gwas2** parameters.
* A reference genome is needed to use this tool. You can access and download these reference genomes from PRS-Csx's [github page](https://github.com/getian107/PRScsx).

* Also you can choose reference panel for polygenic prediction as UKBB (The UK Biobank) or 1KG (1000 Genome Poject) by following parameter: 

```
nextflow run main.nf --ref_dir '$PWD/PRScsx/UKBB'
```

* The format of the GWAS summary statistics should be as follows:

```
SNP	        A1	A2	OR	         P
rs11240767	T	C	0.932267	0.0857218
rs3131972	A	G	0.995723	0.881993
rs1048488	C	T	0.935036	0.0291286
rs12124819	G	A	0.951011	0.113692
rs4040617	G	A	0.961266	0.438928
rs4970383	A	C	0.979289	0.608764
rs13302982	A	G	0.991736	0.775989
rs4040604	G	T	0.996435	0.926046
rs2340587	A	G	0.993043	0.836498
rs28576697	C	T	0.981964	0.701685
rs1110052	T	G	0.994328	0.898852
rs7523549	T	C	0.958352	0.139684
rs3748592	A	G	0.980214	0.667305
rs3748593	A	C	0.999168	0.985898 
```
### MUSSEL
PGSXplorer integrates the MUSSEL tool, a software designed for multi-population polygenic score (PGS) estimation. MUSSEL enables the calculation of polygenic risk scores by utilizing summary statistics from multiple populations and provides improved prediction accuracy across different genetic backgrounds. This pipeline including MUSSEL allows for robust PGS analysis, making it suitable for studies involving heterogeneous populations.
MUSSEL is defined as default:false in PGSXplorer due to the differences in the datasets it uses. The reference genome datasets you need to use for MUSSEL are also available on the MUSSEL [github page](https://github.com/Jin93/MUSSEL).

## Evaluation 
The evaluation of PGS performance is intentionally left flexible, allowing users to adapt this step to their specific research questions and analysis needs.
To help users evaluate PGS performance, we recommend the following commonly used metrics:
**AUC (Area Under the Curve):** A standard metric for case-control datasets, AUC evaluates the classification performance of PGS in distinguishing cases from controls.
**Pseudo R²:** This metric evaluates the agreement between phenotype and PGS in logistic regression models, which is particularly useful for binary outcomes.
**Odds Ratio (OR):** OR measures the increase in risk associated with specific PGS thresholds, making it useful for understanding relative risk.
**Percent Risk:** Identifies high-risk groups by categorizing individuals according to their PGS distributions and helps stratify populations for further analysis.
**Providing Phenotype Data**
PGSXplorer accepts phenotype files as an independent input via the --phenotype parameter. These files should contain binary (for example, 0/1) or continuous values, depending on the type of analysis currently being performed.
Below are sample scripts to help you calculate these metrics.

```
install.packages("pROC")
install.packages("tidyverse")
library(pROC)
library(tidyverse)
# Read files and merge
load_and_merge_data <- function(phenotype_file, pgs_file) {
  phenotype <- read.table(phenotype_file, header = TRUE)
  phenotype$Phenotype <- ifelse(phenotype$Phenotype == 1, 0, 1)
  pgs <- read.table(pgs_file, header = TRUE)
  merged_data <- inner_join(phenotype, pgs, by = "IID")
  return(merged_data)
}
# AUC ve ROC Curve calculation
calculate_auc <- function(merged_data) {
  roc_obj <- roc(merged_data$Phenotype, merged_data$PRS)
  auc_value <- auc(roc_obj)
  print(paste("AUC:", auc_value))
  plot(roc_obj, col = "blue", main = "ROC Curve")
  return(auc_value)
}
# Odds Ratio Calculation
calculate_odds_ratio <- function(merged_data) {
  merged_data$z_scores <- scale(merged_data$PRS)
  model <- glm(Phenotype ~ z_scores, data = merged_data, family = binomial)
  or <- exp(coef(model)["z_scores"])
  print(paste("Odds Ratio per standard deviation:", or))
  return(or)
}
# Determine High Risks
identify_high_risk_individuals <- function(merged_data) {
  risk_threshold <- 3 * median(merged_data$PRS)
  high_risk <- merged_data %>% filter(PRS > risk_threshold)
  print(paste("3 kat veya daha fazla risk taşıyan birey sayısı:", nrow(high_risk)))
  ggplot(merged_data, aes(x = PRS)) +
    geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
    geom_vline(xintercept = risk_threshold, col = "red", linetype = "dashed") +
    annotate("text", x = risk_threshold, y = 20, label = "3x Risk Threshold", color = "red", angle = 90, vjust = -0.5) +
    labs(title = "PGS Risk Distribution", x = "Polygenic Score (PRS)", y = "Frequency")
  return(high_risk)
}
phenotype_file <- "phenotype_file.txt"
pgs_file <- "target_PRSice2.best"
merged_data <- load_and_merge_data(phenotype_file, pgs_file)
auc_result <- calculate_auc(merged_data)
odds_ratio <- calculate_odds_ratio(merged_data)
increased_risk <- calculate_risk_percentiles(merged_data)
high_risk_individuals <- identify_high_risk_individuals(merged_data)
```



