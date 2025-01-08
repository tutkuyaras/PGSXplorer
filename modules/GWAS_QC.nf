process GWAS_QC {
    tag "GWAS QC Steps"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path gwas_sumstat

    output:
    path 'GWAS_QC.*'

    script:
    """
    echo "Input file: ${gwas_sumstat}"
    awk 'NR==1 || (\$6 > 0.01) && (\$14 > 0.8) {print}' < ${gwas_sumstat} |\
    awk '{seen[\$3]++; if(seen[\$3]==1){ print}}' > GWAS_QC.txt
    """

}





