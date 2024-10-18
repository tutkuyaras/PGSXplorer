process ClumpSNPs {
    tag "ClumpSNPs"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path target_9
    path gwas_sumstat

    output:
    path "target_9.clumped"
    
    script:
    """
    plink --bfile target_9 \
          --clump-p1 1 \
          --clump-r2 0.1 \
          --clump-kb 250 \
          --clump ${gwas_sumstat} \
          --clump-snp-field SNP \
          --clump-field P \
          --out target_9
    """
}
