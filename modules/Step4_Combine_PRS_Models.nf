process Step4_Combine_PRS_Models {
    tag "Step 4: Combining PRS models with SuperLearner"
    publishDir "${params.outdir}", mode: 'copy'

    script:
    """
    Rscript ${params.pack}/R/MUSSEL.R \\
        --PATH_package ${params.pack} \\
        --PATH_out ${params.outdir} \\
        --PATH_plink ${params.plink} \\
        --pop ${params.pop} \\
        --target_pop ${params.target_pop} \\
        --chrom ${params.mussel_chrom} \\
        --bfile_tuning ${params.bfile_tuning} \\
        --pheno_tuning ${params.pheno_tuning} \\
        --covar_tuning ${params.covar_tuning} \\
        --bfile_testing ${params.bfile_testing} \\
        --pheno_testing ${params.pheno_testing} \\
        --covar_testing ${params.covar_testing} \\
        --trait_type ${params.trait_type} \\
        --testing TRUE \\
        --NCORES 4
    """
}
