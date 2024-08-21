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
- _Phenotype file_ (in binary format) should be given with the --pheno_file parameter, the file format should be as in the [example](https://github.com/tutkuyaras/PGSExplorer/blob/main/images/pheno_file.png)
  
- _GWAS summary statistics_ should also be given with the --gwas_sumstat parameter, the file format should be as in the example.
