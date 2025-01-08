/*
 * HELP MESSAGE
 */

if (params.help) {
    log.info """
    How to use :

    Required Arguments:
    --target   prefix of plink files (.bed, .bim, .fam)

    Optional Arguments: 
    --mind                  Threshold value to delete SNPs with missingness
    --geno                  Threshold value to delete individuals with missingness 
    --maf                   Threshold value for minimum MAF frequencies of SNPs
    --hwe_case              Threshold value for Hardy-Weinberg Equilibrium for controls
    --hwe_ctrl              Threshold value for Hardy-Weinberg Equilibrium for cases
    --indep_window_size     The window size
    --indep_step_size       The number of SNPs to shift the window at each step
    --indep_threshold       The multiple correlation coefficient for a SNP being regressed on all other SNPs simultaneously.
    --pihat                 The default threshold 0.1875 represents the half-way point between 2nd and 3rd degree relatives
    --relatedness           The same threshold with pihat value
    --ref_phase             Path of Reference BCF folder for Phasing
    --genetic_map           Path of genetic_map file for Phasing
    --ref_imp               Path of Refenrence bref3 folder for Imputation
    --g_map                 Path of Reference genetic map folder for Imputation
    --pca                   Number of principal components to compute, default is 10 
    --pheno_file            Name of the phenotype file located under the target folder
    --run_plink             Run the PLINK part of the workflow, default = true
    --run_prsice            Run the PRSice-2 part of the workflow, default = true
    --run_pca               Run the PCA part of  the workflow if set to true, cretaes ".eigenvec" file
    --run_LDpred2grid       Run the LDPred2 Grid Model of the workflow , default = true
    --run_LDpred2auto       Run the LDPred2 Auto Model of the workflow, default = true
    --run_Lassosum2         Run the Lassosum2 model of the workflow, default = true
    --run_megaprs           Run the MegaPRS model of the workflow, default = true
    --ldak_executable       Path to the LDAK executable (For MegaPRS)
    --mega_model            Model type for MegaPRS. Options: lasso, ridge, bayesr, etc. Default = bayesr
    --run_prscsx            Run the PRScsx of the workflow, default = true
    --run_mussel            Run the MUSSEL of the workflow if set to true, default = false
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
    --meta                  Return combined SNP effect sizes across populations using inverse variance weighted meta-analysis of population-specific posterior effect size estimates. Default is True.             
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
   """
    exit 1
}

println """

    Quality Control Steps of Target Data

    =========================================
    target:           ${params.target}
    outdir:           ${params.outdir}
    visualization :   ${params.graphs}

    --help parameter could be used to see arguments

""".stripIndent()
include {ConvertVCFtoPLINK } from './modules/ConvertVCFtoPLINK'
include { GWAS_QC } from './modules/GWAS_QC'
include { QCSNPmissingness } from './modules/QCSNPmissingness'
include { QCindmissingness } from './modules/QCindmissingness'
include { FilterMAF } from './modules/FilterMAF'
include { PlotMAF } from './modules/PlotMAF'
include { RemoveLowMAF } from './modules/RemoveLowMAF'
include { CalculateHWE } from './modules/CalculateHWE'
include { PlotHWE } from './modules/PlotHWE'
include { InitialHWEFilter } from './modules/InitialHWEFilter'
include { SecondHWEFilter } from './modules/SecondHWEFilter'
include { IBDLDPruning } from './modules/IBDLDPruning'
include { GenerateIBD } from './modules/GenerateIBD'
include { PlotIBD } from './modules/PlotIBD'
include { RelCutoff } from './modules/RelCutoff'
include { FilterIndividuals } from './modules/FilterIndividuals'
include { CalculateHet } from './modules/CalculateHet'
include { PlotHet } from './modules/PlotHet'
include { ProcessHet } from './modules/ProcessHet'
include { CleanHetFailInd } from './modules/CleanHetFailInd'
include { RemoveHetOutliers } from './modules/RemoveHetOutliers'
include { FindDuplicateSNPs } from './modules/FindDuplicateSNPs'
include { RemoveDuplicateSNPs } from './modules/RemoveDuplicateSNPs'
include { ConvertPLINKtoVCF } from './modules/ConvertPLINKtoVCF'
include { SortVCF } from './modules/SortVCF.nf'
include { SplitVCF } from './modules/SplitVCF.nf'
include { PhaseVCF } from './modules/PhaseVCF'
include { Takersid } from './modules/Takersid'
include { ImputeVCF } from './modules/ImputeVCF'
include { MergedVCF } from './modules/MergedVCF'
include { fastmixture } from './modules/fastmixture'
include { VCFtoPLINK } from './modules/VCFtoPLINK'
include { FindDuplicates } from './modules/FindDuplicates'
include { RemoveDuplicates } from './modules/RemoveDuplicates'
include { ClumpSNPs } from './modules/ClumpSNPs'
include { CreateValidSNPs } from './modules/CreateValidSNPs'
include { CreateSNPpvalue } from './modules/CreateSNPpvalue'
include { CreateRangeList } from './modules/CreateRangeList'
include { PruneSNPs } from './modules/PruneSNPs'
include { CalculatePCA } from './modules/CalculatePCA'
include { CalculatePGS } from './modules/CalculatePGS'
include { PRSice2 } from './modules/PRSice2.nf'
include { LDpred2grid } from './modules/LDpred2grid.nf'
include { LDpred2auto } from './modules/LDpred2auto.nf'
include { Lassosum2 } from './modules/Lassosum2.nf'
include { PRSCSx } from './modules/PRSCSx.nf'
include { Step1_Run_LDpred2 } from './modules/Step1_Run_LDpred2.nf'
include { Step2_LDpred2_Tuning } from './modules/Step2_LDpred2_Tuning.nf'
include { Step3_Run_MUSS } from './modules/Step3_Run_MUSS.nf'
include { Step4_Combine_PRS_Models } from './modules/Step4_Combine_PRS_Models.nf'
include { MegaPRS } from './modules/MegaPRS.nf'
include { FormatGWAS } from './modules/FormatGWAS.nf'
include { Impute } from './modules/Impute.nf'
include { Sbayesrc } from './modules/Sbayesrc.nf'

workflow {
    // Choose File Type 

    if (params.vcf_to_plink) {
        println "Files are VCF, Lets convert to PLINK format.."
        // Take VCF file
        vcf_ch = Channel.fromPath("${params.vcf}").ifEmpty { error "VCF file not found: ${params.vcf}" }
        converted_plink_ch = ConvertVCFtoPLINK(vcf_ch)
        target_ch = converted_plink_ch
    } else {
        println "Files are PLINK format"
        // Kullanıcı PLINK formatında dosya belirttiğinde bu süreç çalışır
        target_ch = Channel.fromPath("${params.target_prefix}").ifEmpty { error "PLINK files not found: ${params.target_prefix}" }
    }

    // Define the input channel
    fail_het_qc_ch = Channel.fromPath("fail-het-qc.txt").ifEmpty { error "Input file not found: fail-het-qc.txt" }
    gwas_sumstat_ch = Channel.fromPath("${params.target}/GWAS_sumstat_t1.txt").ifEmpty { error "Input file not found: ${params.target}/GWAS_sumstat_t1.txt" }
    pheno_file_ch = Channel.fromPath("${params.pheno_file}").ifEmpty { error "Input file not found: ${params.pheno_file}"}
    mega_summaries_ch = Channel.fromPath("${params.outdir}/quant.summaries").ifEmpty { error "Input file not found: ${params.outdir}/quant.summaries" }
    ref_phase_ch = Channel.fromPath("${params.ref_phase}/Ref_chr${params.chr}_hg38.bcf.gz").ifEmpty { error "Input file not found: ${params.ref_phase}/Ref_chr${params.chr}_hg38.bcf.gz" }
    genetic_map_ch = Channel.fromPath("$PWD/genetic_map_hg38_withX.txt.gz").ifEmpty { error "Input file not found: $PWD/genetic_map_hg38_withX.txt.gz"}
    formatted_gwas_ch =  Channel.fromPath("${params.outdir}/formatted_gwas.ma").ifEmpty { error "Input file not found: ${params.outdir}/formatted_gwas.ma" }
    impute_out_ch =  Channel.fromPath("$PWD/sbayesrc.imputed.ma").ifEmpty { error "Input file not found: $PWD/sbayesrc.imputed.ma"}

    // Execute QC Part

    // Execute GWAS QC
    gwas_qc = GWAS_QC(gwas_sumstat_ch)

    // Execute QCSNPmissingness
    qcsnp_out = QCSNPmissingness(target_ch)
    
    // Execute QCindmissingness
    qcind_out = QCindmissingness(qcsnp_out) 

    // Execute FilterMAF
    maf_out = FilterMAF(qcind_out)
    
    // Separate the outputs
    merged_bfile_ch = maf_out.map { it[0] }
    maf_frq_ch = maf_out.map { it[1] }
    
    // Execute RemoveLowMAF
    remove_low_maf_out = RemoveLowMAF(merged_bfile_ch.first())
    
    // Execute PlotMAF
    plot_out = PlotMAF(maf_frq_ch)

    // Execute CalculateHWE
    hwe_out = CalculateHWE(remove_low_maf_out)

    // Execute PlotHWE
    plot_hwe_out = PlotHWE(hwe_out)

    // Execute InitialHWEFilter 
    initial_hwe_filter_out = InitialHWEFilter(remove_low_maf_out)

    // Execute SecondHWEFilter
    second_hwe_filter_out = SecondHWEFilter(initial_hwe_filter_out) 

    // Execute IBD-LD-pruning
    ibd_ld_pruning_out = IBDLDPruning(second_hwe_filter_out)

    // Extract the prune.in file
    prune_in_file_ch = ibd_ld_pruning_out.prune_in

    // Execute GenerateIBD
    generate_ibd_out = GenerateIBD(second_hwe_filter_out, prune_in_file_ch)

    // Execute PlotIBD
    plot_ibd_out = PlotIBD(generate_ibd_out)

    // Execute RelCutoff
    rel_cutoff_out = RelCutoff(second_hwe_filter_out, prune_in_file_ch)

    // Execute FilterIndividuals
    filter_individuals_out = FilterIndividuals(second_hwe_filter_out, rel_cutoff_out)

    // Execute CalculateHet
    calculate_het_out = CalculateHet(filter_individuals_out, prune_in_file_ch)

    // Execute PlotHet
    plot_het_out = PlotHet(calculate_het_out)

    // Execute ProcessHet
    process_het_out = ProcessHet(calculate_het_out)

    // Execute CleanHetFailInd
    clean_het_fail_ind_out = CleanHetFailInd(fail_het_qc_ch)

    // Execute RemoveHetOutliers
    remove_het_outliers_out = RemoveHetOutliers(filter_individuals_out, clean_het_fail_ind_out )

    // Execute FindDuplicatesSNPs
    find_duplicate_snps_out = FindDuplicateSNPs(remove_het_outliers_out)

    // Execute RemoveDuplicateSNPs
    remove_duplicate_snps_out = RemoveDuplicateSNPs(remove_het_outliers_out, find_duplicate_snps_out )
    
    // LD Pruning
    prune_snps_out = PruneSNPs(remove_duplicate_snps_out)

    // Convert PLINK to VCF
    vcf_out = ConvertPLINKtoVCF(remove_duplicate_snps_out)
    
    // Sort VCF Files
    sorted_vcf_ch = SortVCF(input_vcf = vcf_out)

    // Split VCF into chromosomes
    chr_ch = Channel.of(1..22)
    
    split_ch = SplitVCF(vcf_file = sorted_vcf_ch.first(),
                        chr = chr_ch.flatten())

    // Phase VCF file
    // Map split VCF files to corresponding reference files
    split_ref_ch = split_ch
        .map { splitted_vcf ->
        // Extract chromosome number from the split VCF filename
        def chr = splitted_vcf.name.replaceAll(/.*vcfbychrom_(\d+)\.vcf\.gz/, '$1') 
        def ref_bcf = file("${params.ref_phase}/Ref_chr${chr}_hg38.bcf.gz")
        [splitted_vcf, ref_bcf, chr]
    }
    .filter { tuple -> tuple[0] && tuple[1].exists() } // Ensure both files exist

    phased_vcf_ch = PhaseVCF(
        splitted_vcf = split_ref_ch.map { it[0] },
        ref_bcf = split_ref_ch.map { it[1] },
        genetic_map = file("${params.genetic_map}"),
        chr = split_ref_ch.map { it[2] }
    )
    phased_vcf_ch.view { "Processing (Phasing): Phased VCF=${it[0]}, Ref=${it[1]}, Chr=${it[2]}" }

    // Takersid
    takersid_input_ch = phased_vcf_ch
        .map { phased_vcf ->
            // Extract chromosome number from phased VCF filename
            def chr = phased_vcf.name.replaceAll(/.*phasedvcf_chr(\d+)\.vcf\.gz/, '$1')
            [phased_vcf, chr]
    }
    .filter { tuple -> tuple[0].exists() } 

    rsid_map_ch = Takersid (
        phased_vcf = takersid_input_ch.map { it[0] },
        chr = takersid_input_ch.map { it[1] }
    )
  
    // Imputation 
    // Map phased VCF files and reference files for imputations and combine rsid
    imputation_input_ch = phased_vcf_ch
        .map { phased_vcf ->
        def chr = phased_vcf.name.replaceAll(/.*phasedvcf_chr(\d+)\.vcf\.gz/, '$1')
        def ref_bref = file("${params.ref_imp}/Ref_chr${chr}_hg38.bref3")
        def g_map = file("${params.g_map}/plink.chr${chr}.GRCh38.map")
        [phased_vcf, ref_bref, g_map, chr]
    }
    .filter { tuple -> tuple[0] && tuple[1].exists() && tuple[2].exists() } // Ensure all files exist

    // Debug to view matched files
    imputation_input_ch.view { 
        "Processing (Imputation): Phased VCF=${it[0]}, Ref=${it[1]}, Map=${it[2]}, Chr=${it[3]}" 
    }

    // Run imputation process
    imputed_vcf_ch = ImputeVCF(
        phased_vcf = imputation_input_ch.map { it[0] },
        ref_bref = imputation_input_ch.map { it[1] },
        g_map = imputation_input_ch.map { it[2] },
        chr = imputation_input_ch.map { it[3] },
        rsid_map = rsid_map_ch
    )

    
    // Merge VCF Files 
    merged_vcf_ch = MergedVCF(imputed_vcf_ch.collect())

    // Convert Imputed VCF to PLINK format
    plink_out = VCFtoPLINK(merged_vcf_ch)

    // Remove Duplicates
    find_duplicates = FindDuplicates(plink_out)

    last_plink_out = RemoveDuplicates(plink_out, find_duplicates)

    // Target Ancestry Inference
    ancestry_inf = fastmixture(last_plink_out)
    
    // Conditionally execute PLINK processes
    if (params.run_plink) {
        clump_snps_out = ClumpSNPs(last_plink_out, gwas_qc)
        create_valid_snps_out = CreateValidSNPs(clump_snps_out)
        create_snp_pvalue_out = CreateSNPpvalue(gwas_qc)
        create_range_list_out = CreateRangeList()
        calculate_pgs_out = CalculatePGS(
            last = last_plink_out,
            gwas = gwas_qc,
            range_list = create_range_list_out,
            valid_snp = create_valid_snps_out,
            snp_pvalue = create_snp_pvalue_out
        )
    }

    if (params.run_pca) { 
        // Execute CalculatePCA and extract .eigenvec file
        calculate_pca_out = CalculatePCA(last_plink_out, prune_snps_out)

        // Separate the outputs
        eigenvec_ch = calculate_pca_out.map { it[0] }
        eigenval_ch = calculate_pca_out.map { it[1] }
    }

     // Conditionally execute PRSice-2 processes
    if (params.run_prsice) {
        prsice_out = PRSice2(
            last = last_plink_out,
            gwas = gwas_qc,
            pheno_file = pheno_file_ch,
            eigenvec = eigenvec_ch,
        )
    }

    // Conditionally execute LDpred grid process
    if (params.run_LDpred2grid) {
        bed_ch = last_plink_out.map { it[0] }

        LDpred2grid(
            pheno_file = pheno_file_ch,
            eigenvec = eigenvec_ch,
            gwas = gwas_qc,
            last = last_plink_out
             
        )
    }

     // Conditionally execute LDpred auto process
    if (params.run_LDpred2auto) {
        bed_ch = last_plink_out.map { it[0] }

        LDpred2auto(
            pheno_file = pheno_file_ch,
            eigenvec = eigenvec_ch,
            gwas = gwas_qc,
            last = last_plink_out,
             
        )
    }

    // Conditionally execute Lassosum2 process
    if (params.run_Lassosum2) {
        bed_ch = last_plink_out.map { it[0] }

        Lassosum2(
            pheno_file = pheno_file_ch,
            eigenvec = eigenvec_ch,
            gwas = gwas_qc,
            last = last_plink_out,
             
        )
    }
   
    // Conditionally execute PRS-CSx process
    if (params.run_prscsx) {
        PRSCSx(
            last = last_plink_out
        )
    }

    // Conditionally execute MUSSEL process 
    if (params.run_mussel) {
        ldpred2_outputs = Step1_Run_LDpred2()
        ldpred2_tuned_outputs = Step2_LDpred2_Tuning(ldpred2_outputs)
        muss_outputs = Step3_Run_MUSS()
        Step4_Combine_PRS_Models()
    }

    // Conditionally execute MegaPRS process

    if (params.run_megaprs) {
        mega_prs_out = MegaPRS(
            target_qc_prefix = last_plink_out,
            pheno_file = pheno_file_ch,
            mega_summaries = mega_summaries_ch)
    }

    if (params.run_sbayesrc) {
        formatted_gwas = FormatGWAS(gwas_qc)
        impute_out = Impute(
                ld_folder = file(params.ld_folder),
                ma_file = formatted_gwas_ch,
                out_prefix = params.out_prefix,
                threads = params.threads
        )
       sbayesr_out = Sbayesrc(
                ld_folder = file(params.ld_folder),
                imp_file = impute_out_ch,
                annot = file(params.annot_file),
                out_prefix = params.out_prefix,
                threads = params.threads
       )
    }
}




