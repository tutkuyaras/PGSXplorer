process MergedVCF {
    tag "Convert Last QC ed VCF files into PLINK format"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path imputed_vcfs

    output:
    path "merged_output.vcf.gz"

    script:
    """
    for vcf in ${imputed_vcfs}; do
        echo "Indexing \$vcf"
        bcftools index \$vcf
    done
    
    bcftools concat ${imputed_vcfs.join(' ')} -Oz -o merged_output.vcf.gz
    """
}



