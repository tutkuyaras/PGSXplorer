process Impute {
    tag "Run SBayesR-C"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path ld_folder
    path ma_file
    val out_prefix
    val threads

    script:
    """
    # Main: SBayesRC
    docker pull zhiliz/sbayesrc
    docker run -v ${PWD}:/data zhiliz/sbayesrc --ldm-eigen /data/${ld_folder} --gwas-summary /data/${ma_file} --impute-summary --out /data/${out_prefix} --threads ${threads}
    """

} 