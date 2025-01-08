process SplitVCF {
    tag "Split VCF files into Chromosomes"
    tag { "chr${chr}" }
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path vcf_file 
    val chr
    
    output:
    path "vcfbychrom_${chr}.vcf.gz"
 

    script:
    """
    tabix -p vcf ${vcf_file}
    bcftools view ${vcf_file} --regions ${chr} -o vcfbychrom_${chr}.vcf.gz -Oz
    """
}




