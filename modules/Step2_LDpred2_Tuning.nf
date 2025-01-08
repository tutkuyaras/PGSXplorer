process Step2_LDpred2_Tuning {
    tag "Step 2: LDpred2 tuning is running"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path ldpred2_outputs

    output:
    tuple path ("EUR_optim_params.txt"), 
          path ("AFR_optim_params.txt"), emit: ldpred2_tuned_outputs

    script:
    """
    Rscript ${params.pack}/R/LDpred2_tuning.R \\
        --PATH_package ${params.pack} \\
        --PATH_out ${params.outdir} \\
        --PATH_plink ${params.plink} \\
        --FILE_sst ${params.sst} \\
        --pop ${params.pop} \\
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

    # Dosyaları işlem dizinine kopyalayın
    cp ${params.outdir}/LDpred2/EUR_optim_params.txt EUR_optim_params.txt
    cp ${params.outdir}/LDpred2/AFR_optim_params.txt AFR_optim_params.txt
    """
}
