process ConvertPLINKtoVCF {
    tag "Convert PLINK to VCF and index files"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path target_9

    output:
    path "*.vcf"

    script:
    """
    plink --bfile target_9 --recode vcf --out target_9
    """
}
