process PlotIBD {
    tag "Plot IBD distribution"
    publishDir "${params.graphs}", mode: 'copy'

    input:
    path ibd_file

    output:
    path "ibd_plot.pdf"

    script:
    """
    Rscript $PWD/bin/plot_ibd.R ${ibd_file}
    """
}
