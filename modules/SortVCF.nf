process SortVCF {
    tag "Sort VCF File"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path input_vcf

    output:
    path "*"

    script:
    """
    bcftools sort ${input_vcf} -Oz -o sorted.vcf.gz
    tabix -p vcf sorted.vcf.gz
    """
}
