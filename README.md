# PGSExplorer
PGSExplorer is a bioinformatics workflow designed to calculate polygenic scores by processing genomic data through quality control steps. Optionally, it can utilize tools such as [PLINK](https://www.cog-genomics.org/plink/), [PRSice-2](https://choishingwan.github.io/PRSice/), [LD-Pred2 (grid)](https://privefl.github.io/bigsnpr/articles/LDpred2.html), [LD-Pred2 (auto)](https://privefl.github.io/bigsnpr/articles/LDpred2.html), [Lassosum2](https://privefl.github.io/bigsnpr/articles/LDpred2.html#lassosum2-grid-of-models), [PRS-CSx](https://github.com/getian107/PRScsx), and [MUSSEL](https://github.com/Jin93/MUSSEL). The workflow requires genomic files in PLINK format (.bed, .bim, .fam) and GWAS summary statistics for two different populations as input to complete the analysis.


![PGSExplorer Diagram](https://github.com/tutkuyaras/PGSExplorer/blob/main/images/PGSExplorer%20Workflow.drawio.png)
## Workflow Overview

PGSExplorer includes a comprehensive pipeline that begins with rigorous quality control (QC) measures to ensure the integrity of genomic data. Following the completion of the QC module, users can optionally execute various polygenic score (PGS) calculation tools. The steps are as follows:

## Workflow Overview

PGSExplorer includes a comprehensive pipeline that begins with rigorous quality control (QC) measures to ensure the integrity of genomic data. Following the completion of the QC module, users can optionally execute various polygenic score (PGS) calculation tools. The steps are as follows:

**1. Filtering Missing SNPs:** Identify and remove SNPs with missing genotype data.  
**2. Filtering Missing Individuals:** Exclude individuals with excessive missing genotype data.  
**3. Filtering by Minor Allele Frequency (MAF):** Retain SNPs above a certain MAF threshold.  
**4. Visualization of MAF Distributions:** Graphical representation of MAF across SNPs.  
**5. Filtering by Hardy-Weinberg Equilibrium (HWE):** Remove SNPs that deviate significantly from HWE.  
**6. Visualization of HWE Distributions:** Visual display of HWE p-values for SNPs.  
**7. Relatedness Checking:** Assess genetic relatedness between individuals.  
**8. Visualization of Identity by Descent (IBD):** Visualize IBD statistics to identify related individuals.  
**9. Heterozygosity Assessment:** Evaluate heterozygosity rates to identify potential outliers.  
**10. Visualization of Heterozygosity Distributions:** Plot distribution of heterozygosity rates.  
**11. Removal of Duplicate SNPs:** Eliminate duplicated SNPs to prevent redundancy.  
**12. Calculating 10 Principal Components (PCA):** Perform PCA to capture population structure.  
**13. Pruning and Thresholding (P+T):** Execute P+T method using PLINK to calculate PGS.   
**14. PRSice-2:** Calculate and visualize PGS using PRSice-2.  
**15. LD-Pred2-grid :** Apply LD-Pred2 grid model for PGS estimation.  
**16. LD-Pred2-auto :** Apply LD-Pred2 auto model for PGS estimation.   
**17.Lassosum2 :** Apply Lassosum2 for PGS estimation.   
**18. PRS-CSx:**  Implement PRS-CSx for multi-ancestry PGS.    
**19. MUSSEL:** Utilize MUSSEL for multi-ancestry PGS estimation.    


## Usage
In order for PGSExplorer to work correctly, you must first examine the parameters of the quality control steps and the polygenic score calculation tools you want to run and start the analysis with the appropriate command. You can also access the parameters given below by running the command 
```
nextflow run main.nf --help
```
  
### Parameters
```
    Required Arguments:
    --target   prefix of plink files (.bed, .bim, .fam)

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
    
    // PGS Parameters
    --run_plink             Run the PLINK part of the workflow if set to true
    --run_prsice            Run the PRSice-2 part of the workflow if set to true
    --run_pca               Run the PCA part of  the workflow if set to true, cretaes ".eigenvec" file
    --run_LDpred2grid       Run the LDPred2 Grid Model of the workflow if set to true
    --run_LDpred2auto       Run the LDPred2 Auto Model of the workflow if set to true
    
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
    
    // MUSSEL parameters
    
    
```
## Quality Control Part

Quality control modules consist of seven basic steps which are not optional. It consists of SNP filtering, individual filtering, filtering by MAF, HWE, relatedness, heterozygosity and elimination of duplicate SNPs.   
The files obtained after each step are saved in the QC_outputs folder. In addition, the distribution graphs of MAF, HWE, heterozygosity and kinship relationships are also saved in the QC_graphs folder.
  
Only for quality control module, you can run the pipeline using: 

```
nextflow run main.nf --run_prsice false --run_LDpred2grid false --run_LDpred2auto false --run_Lassosum2 false --run_prscsx false --run_mussel false
```

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

> **NOTE:** For PGSExplorer to work correctly, file formats, column order and names must be the same.  


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
### PCA Calculation
In order to use PRSice-2, LD-Pred2 grid, LD-Pred2 auto, LD-Pred2 auto, Lassosum2 PGS models, the .eigenvec file used as covariate is needed. This file is created with the pca module in PGSExplorer. The default pca value is 10. You can change this by using the **--pca** parameter. 

### PRSice-2
PRSice-2 is a comprehensive software tool for calculating polygenic risk scores (PRS) by integrating genome-wide association study (GWAS) summary statistics with individual genotype data. It provides flexible options for clustering and thresholding to optimize PRS structure, enabling analysis of genetic susceptibility to various traits and diseases. You can check details from [here](https://github.com/choishingwan/PRSice)

### LD-Pred2 Grid & LD-Pred2 Auto and Lassosum2 
LDpred2 is a Bayesian polygenic risk score (PRS) tool that includes two primary models: LDpred2-grid and LDpred2-auto. LDpred2-grid explores a grid of hyperparameters to find the best-fit PRS model, while LDpred2-auto automatically adjusts parameters based on the data, eliminating the need for predefined hyperparameters. Lassosum2 is another PRS method that applies penalized regression techniques, effectively handling linkage disequilibrium (LD) patterns to improve the prediction accuracy of polygenic scores. Yo can check details from [here](https://privefl.github.io/bigsnpr/articles/LDpred2.html) 


## Multi Ancestry PGS Tools  

