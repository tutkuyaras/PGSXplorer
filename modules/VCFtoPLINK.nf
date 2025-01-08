process VCFtoPLINK {
    tag "Convert Imputed VCF to PLINK"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path merged

    output:
    path "merged.*"

    script:
    """
    # Step 1: Index the VCF file
    bcftools sort ${merged} -Oz -o merged_sorted.vcf.gz
    tabix -p vcf merged_sorted.vcf.gz || { echo "Tabix indexing failed"; exit 1; }

    # Step 2: Convert VCF to PLINK format
    plink --vcf merged_sorted.vcf.gz --allow-extra-chr --make-bed --out merged
    """
}
