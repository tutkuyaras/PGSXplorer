process Step1_Run_LDpred2 {
    tag "Step 1: LDpred2 is running"
    publishDir "${params.out}", mode: 'copy'

    output:
    path "*", emit: ldpred2_outputs

    script:
    """
    Rscript ${params.pack}/R/LDpred2_jobs.R \\
        --PATH_package ${params.pack} \\
        --PATH_data ${params.data} \\
        --PATH_LDref ${params.LDref} \\
        --PATH_out ${params.outdir} \\
        --FILE_sst ${params.sst} \\
        --pop ${params.pop} \\
        --bfile_tuning ${params.bfile_tuning} \\
        --NCORES 20

    while [ ! -f "${params.outdir}/EUR/tmp/beta_files/beta_in_all_settings/params.txt" ] || [ ! -f "${params.outdir}/AFR/tmp/beta_files/beta_in_all_settings/params.txt" ]; do
        echo "Files are writing..."
        sleep 60
    done

    echo "Files are ready, proceeding."
    
    """
}