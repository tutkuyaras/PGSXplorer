process MegaPRS {
    tag "Run MegaPRS Workflow"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path target_qc_prefix
    path pheno_file
    path mega_summaries

    output:
    path "*"

    script:
    """
    # Run MegaPRS logistic regression
    ${params.ldak_executable} --logistic quant --bfile ${params.target_qc_prefix} --pheno ${params.pheno_file} 
    sleep 60
    ${params.ldak_executable} --calc-cors cors --bfile ${params.target_qc_prefix}
    ${params.ldak_executable} --mega-prs mega --model ${params.mega_model} --summary ${params.mega_summaries} --power -0.25 --cors cors --allow-ambiguous YES
 
    """
}









