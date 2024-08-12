# PGSExplorer
PGSExplorer is a bioinformatics workflow designed to calculate polygenic scores by processing genomic data through quality control steps. Optionally, it can utilize tools such as [PLINK](https://www.cog-genomics.org/plink/), [PRSice-2](https://choishingwan.github.io/PRSice/), [LD-Pred2 (grid)](https://privefl.github.io/bigsnpr/articles/LDpred2.html), [LD-Pred2 (auto)](https://privefl.github.io/bigsnpr/articles/LDpred2.html), [PRS-CSx](https://github.com/getian107/PRScsx), and [MUSSEL](https://github.com/Jin93/MUSSEL). The workflow requires genomic files in PLINK format (.bed, .bim, .fam) and GWAS summary statistics for two different populations as input to complete the analysis.


![PGSExplorer Diagram](https://github.com/tutkuyaras/PGSExplorer/blob/main/PGSExplorer%20Workflow.png)
## Workflow Overview

PGSExplorer includes a comprehensive pipeline that begins with rigorous quality control (QC) measures to ensure the integrity of genomic data. Following the completion of the QC module, users can optionally execute various polygenic score (PGS) calculation tools. The steps are as follows:

**1.  Filtering Missing SNPs:** Identify and remove SNPs with missing genotype data.  
**2.  Filtering Missing Individuals:** Exclude individuals with excessive missing genotype data.  
**3.  Filtering by Minor Allele Frequency (MAF):** Retain SNPs above a certain MAF threshold.  
**4.  Visualization of MAF Distributions:** Graphical representation of MAF across SNPs.  
**5.  Filtering by Hardy-Weinberg Equilibrium (HWE):** Remove SNPs that deviate significantly from HWE.  
**6.  Visualization of HWE Distributions:** Visual display of HWE p-values for SNPs.  
**7.  Relatedness Checking:** Assess genetic relatedness between individuals.  
**8.  Visualization of Identity by Descent (IBD):** Visualize IBD statistics to identify related individuals.    
**9.  Heterozygosity Assessment:** Evaluate heterozygosity rates to identify potential outliers.  
**10. Visualization of Heterozygosity Distributions:** Plot distribution of heterozygosity rates.  
**11. Removal of Duplicate SNPs:** Eliminate duplicated SNPs to prevent redundancy.  
**12. Calculating 10 Principal Components (PCA):** Perform PCA to capture population structure.  
**13. Pruning and Thresholding (P+T):** Execute P+T method using PLINK to calculate PGS.  
**14. PRSice-2:** Calculate and visualize PGS using PRSice-2.  
**15. LD-Pred2 :** Apply LD-Pred2 (grid and auto)  for PGS estimation.  
**16. PRS-CSx:**  Implement PRS-CSx for multi-ancestry PGS.  
**17. MUSSEL:** Utilize MUSSEL for multi-ancestry PGS estimation.  
