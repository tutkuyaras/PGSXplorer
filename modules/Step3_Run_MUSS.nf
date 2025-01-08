process Step3_Run_MUSS {
    tag "Step 3: Running MUSS"
    publishDir "${params.outdir}", mode: 'copy'

    output:
    path "MUSS_beta_in_all_settings_bychrom/*", emit: muss_outputs

    script:
    """
    Rscript ${params.pack}/R/MUSS_jobs.R \\
        --PATH_package ${params.pack} \\
        --PATH_data ${params.data} \\
        --PATH_LDref ${params.LDref} \\
        --PATH_out ${params.outdir} \\
        --FILE_sst ${params.sst} \\
        --pop ${params.pop} \\
        --LDpred2_params ${params.outdir}/LDpred2/EUR_optim_params.txt,${params.outdir}/LDpred2/AFR_optim_params.txt \\
        --chrom ${params.mussel_chrom} \\
        --bfile_tuning ${params.bfile_tuning} \\
        --NCORES 5

    while [ ! -f "${params.outdir}/tmp/MUSS_beta_in_all_settings_bychrom/settings_22.txt" ] ; do
        echo "Files are writing..."
        sleep 60
    done

    echo "Files are ready, proceeding."

    # Copy the output files to the current working directory
    mkdir -p MUSS_beta_in_all_settings_bychrom
    cp ${params.outdir}/tmp/MUSS_beta_in_all_settings_bychrom/*.txt MUSS_beta_in_all_settings_bychrom/
    """
}
