process ClumpSNPs {
    tag "ClumpSNPs"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path last
    path gwas

    output:
    path "last.clumped"
    
    script:
    """
    plink --bfile last \
          --clump-p1 1 \
          --clump-r2 0.1 \
          --clump-kb 250 \
          --clump ${gwas} \
          --clump-snp-field SNP \
          --clump-field P \
          --out last
    """
}
