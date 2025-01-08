process CreateRangeList {
    tag "CreateRangeList for P+T"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    output:
    path "range_list"

    script:
    """
    echo "0.001 0 0.001" > range_list 
    echo "0.05 0 0.05" >> range_list
    echo "0.1 0 0.1" >> range_list
    echo "0.2 0 0.2" >> range_list
    echo "0.3 0 0.3" >> range_list
    echo "0.4 0 0.4" >> range_list
    echo "0.5 0 0.5" >> range_list
    """
}
