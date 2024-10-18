process Step1_Run_LDpred2 {
    tag "Step 1: LDpred2 is running"
    publishDir "${params.out}", mode: 'copy'

    output:
    path "*", emit: ldpred2_outputs

    script:
    """
    Rscript ${params.package}/R/LDpred2_jobs.R \\
        --PATH_package ${params.package} \\
        --PATH_data ${params.data} \\
        --PATH_LDref ${params.LDref} \\
        --PATH_out ${params.out} \\
        --FILE_sst ${params.sst} \\
        --pop ${params.pop} \\
        --bfile_tuning ${params.bfile_tuning} \\
        --NCORES 20

    while [ ! -f "${params.out}/EUR/tmp/beta_files/beta_in_all_settings/params.txt" ] || [ ! -f "${params.out}/AFR/tmp/beta_files/beta_in_all_settings/params.txt" ]; do
        echo "Files are writing..."
        sleep 60
    done

    echo "Files are ready, proceeding."
    
    """
}

process Step2_LDpred2_Tuning {
    tag "Step 2: LDpred2 tuning is running"
    publishDir "${params.out}", mode: 'copy'

    input:
    path ldpred2_outputs

    output:
    tuple path ("EUR_optim_params.txt"), 
          path ("AFR_optim_params.txt"), emit: ldpred2_tuned_outputs

    script:
    """
    Rscript ${params.package}/R/LDpred2_tuning.R \\
        --PATH_package ${params.package} \\
        --PATH_out ${params.out} \\
        --PATH_plink ${params.plink} \\
        --FILE_sst ${params.sst} \\
        --pop ${params.pop} \\
        --chrom ${params.chrom} \\
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
    cp ${params.out}/LDpred2/EUR_optim_params.txt EUR_optim_params.txt
    cp ${params.out}/LDpred2/AFR_optim_params.txt AFR_optim_params.txt
    """
}

process Step3_Run_MUSS {
    tag "Step 3: Running MUSS"
    publishDir "${params.out}", mode: 'copy'

    input:
    tuple path(eur_optim_params), path(afr_optim_params)


    output:
    tuple path ("beta_file_AFR.txt"), 
          path ("beta_file_EUR.txt"),
          path ("beta_file_all.txt"),
          path ("beta_settings.txt"), emit: muss_outputs

    script:
    """
    Rscript ${params.package}/R/MUSS_jobs.R \\
        --PATH_package ${params.package} \\
        --PATH_data ${params.data} \\
        --PATH_LDref ${params.LDref} \\
        --PATH_out ${params.out} \\
        --FILE_sst ${params.sst} \\
        --pop ${params.pop} \\
        --LDpred2_params ${eur_optim_params},${afr_optim_params} \\
        --chrom ${params.chrom} \\
        --bfile_tuning ${params.bfile_tuning} \\
        --NCORES 5
    
    while [ ! -f "${params.out}/MUSS/beta_file_AFR.txt" ] || [ ! -f "${params.out}/MUSS/beta_file_AFR.txt" ] || [ ! -f "${params.out}/MUSS/beta_file_all.txt" ]  ; do
        echo "Files are writing..."
        sleep 60
    done

    echo "Files are ready, proceeding."
    """
}


process Step4_Combine_PRS_Models {
    tag "Step 4: Combining PRS models with SuperLearner"
    publishDir "${params.out}", mode: 'copy'

    input:
    path muss_outputs

    script:
    """
    Rscript ${params.package}/R/MUSSEL.R \\
        --PATH_package ${params.package} \\
        --PATH_out ${params.out} \\
        --PATH_plink ${params.plink} \\
        --pop ${params.pop} \\
        --target_pop ${params.target_pop} \\
        --chrom ${params.chrom} \\
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

workflow {
    // Step 1: LDpred2 Çalıştırma
    ldpred2_outputs = Step1_Run_LDpred2()

    // Step 2: LDpred2 Tuning işlemi
    ldpred2_tuned_outputs = Step2_LDpred2_Tuning(ldpred2_outputs)

    // Step 3: MUSS işlemi
    muss_outputs = Step3_Run_MUSS(ldpred2_tuned_outputs)

    // Step 4: SuperLearner PRS Modellerini Birleştirme
    Step4_Combine_PRS_Models(muss_outputs)
}
