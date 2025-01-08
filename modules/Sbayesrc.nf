process Sbayesrc {
    tag "Run SBayesR-C"
    debug true
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path ld_folder
    path imp_file
    path annot
    val out_prefix
    val threads

    script:
    """
    # Main: SBayesRC
    docker pull zhiliz/sbayesrc
    docker run -v ${PWD}:/data zhiliz/sbayesrc --ldm-eigen /data/${params.ld_folder} --gwas-summary /data/${imp_file}  --sbayes RC --annot /data/${params.annot_file}  --out /data/${params.out_prefix} --threads ${params.threads}
    """
}
