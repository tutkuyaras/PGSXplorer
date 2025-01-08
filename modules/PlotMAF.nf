process PlotMAF {
    tag "Generate MAF distribution plot"
    debug true
    publishDir "${params.graphs}", mode: 'copy'

    input:
    path maf_freq_file

    output:
    path 'MAF_dist.pdf'

    script:
    """
    Rscript $PWD/bin/plot_maf.R ${maf_freq_file}
    """
}
