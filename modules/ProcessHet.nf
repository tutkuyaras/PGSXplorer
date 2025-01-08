process ProcessHet {
    tag "Process heterozygosity data"
    debug true
    publishDir "${params.graphs}", mode: 'copy'

    input:
    path het_file

    output:
    path "het.pdf"

    script:
    """
    Rscript $PWD/bin/process_het.R ${het_file}
    """
}
