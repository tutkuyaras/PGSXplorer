process Takersid {
    tag "Takingrsid"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path phased_vcf
    val chr

    output:
    path "rsid_map_chr${chr}.txt"

    script:
    """
    zcat ${phased_vcf} | grep -v '^#' | awk '{print \$1, \$2, \$3}' > rsid_map_chr${chr}.txt
    """
}









