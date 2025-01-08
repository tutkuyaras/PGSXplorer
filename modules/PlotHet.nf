process PlotHet {
    tag "Plot heterozygosity rate"
    publishDir "${params.graphs}", mode: 'copy'

    input:
    path het_file

    output:
    path "het.pdf"

    script:
    """
    Rscript $PWD/bin/plot_het.R ${het_file}
    """
}
