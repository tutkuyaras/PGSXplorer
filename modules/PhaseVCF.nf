process PhaseVCF {
    tag "Phasing of VCF file"
    tag { "chr${chr}" }
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path splitted_vcf
    path ref_bcf
    path genetic_map
    val chr

    output:
    path "phasedvcf_chr${chr}.vcf.gz"
   

    script:
    """
    bcftools index ${ref_bcf}
    tabix -p vcf ${splitted_vcf} 
    $PWD/eagle \
        --vcfTarget ${splitted_vcf} \
        --vcfRef ${ref_bcf} \
        --geneticMapFile ${genetic_map} \
        --vcfOutFormat z \
        --numThreads 24 \
        --outPrefix phasedvcf_chr${chr}


    
    """
}
