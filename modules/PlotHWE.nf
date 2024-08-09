process PlotHWE {
    tag "Plotting Hardy-Weinberg Equilibrium"
    debug true
    publishDir "${params.graphs}", mode: 'copy'

    input:
    path hwe_file

    output:
    path "hwe.pdf"

    script:
    """
    Rscript $PWD/bin/plot_hwe.R ${hwe_file}
    """
}
