process PRSCSx {
    tag "PRS-CSx"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path target_9
    

    output:
    path "*.txt"

    script:
    """
    echo "Output directory: ${params.outdir}"
    python /sw/bip_apps/PRScsx/PRScsx.py \
        --ref_dir=${params.prscsxref_dir} \
        --bim_prefix=target_9 \
        --sst_file=${params.prscsx_gwas1},${params.prscsx_gwas2} \
        --n_gwas=${params.n_gwas1},${params.n_gwas2}  \
        --pop=${params.pop1},${params.pop2} \
        --chrom=${params.chrom} \
        --phi=${params.phi} \
        --meta=${params.meta} \
        --out_dir=./ \
        --out_name=${params.out_name}
    """
}
