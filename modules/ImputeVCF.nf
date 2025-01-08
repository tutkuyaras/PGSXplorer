process ImputeVCF {
    tag "Imputation for chromosome ${chr} and QC"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path phased_vcf
    path ref_bref
    path g_map
    val chr
    path rsid_map

    output:
    path "imputed_chr${chr}_with_rsid.vcf.gz"
   
    script:
    """
    tabix -p vcf ${phased_vcf} 
    beagle ref=${ref_bref} gt=${phased_vcf} map=${g_map} out=imputed_chr${chr}
    bcftools index imputed_chr${chr}.vcf.gz
    bcftools view -O z -o imputed_chr${chr}_QC.vcf.gz -e 'INFO/DR2<0.8' imputed_chr${chr}.vcf.gz
    bcftools index imputed_chr${chr}_QC.vcf.gz

    awk 'BEGIN {OFS="\t"} {if (NR > 1) print \$1, \$2, \$3, ".", ".", ".", ".", "."}' ${rsid_map} > rsid${chr}.vcf
    bgzip rsid${chr}.vcf
    tabix -p vcf rsid${chr}.vcf.gz
    bcftools annotate -a rsid${chr}.vcf.gz -c CHROM,POS,ID -o imputed_chr${chr}_with_rsid.vcf -O v imputed_chr${chr}_QC.vcf.gz
    bgzip imputed_chr${chr}_with_rsid.vcf
    bcftools index imputed_chr${chr}_with_rsid.vcf.gz


    sleep 200
    echo "Files are ready, proceeding."
    """
}



