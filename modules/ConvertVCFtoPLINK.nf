process ConvertVCFtoPLINK {
    tag "Convert VCF to PLINK"
    publishDir "${params.target}", mode: 'copy'

    input:
    path vcf_file

    output:
    path "target.*"

    script:
    """
    # Step 1: Index the VCF file
    bcftools sort ${vcf_file} -Oz -o target.vcf.gz
    tabix -p vcf target.vcf.gz || { echo "Tabix indexing failed"; exit 1; }

    # Step 2: Convert VCF to PLINK format
    plink --vcf target.vcf.gz --allow-extra-chr --make-bed --out target || { echo "PLINK conversion failed"; exit 1; }
    """
}
