process GENERATE_GENOME {
    tag "${meta} -> my.genome"
    label "process_small"

    input:
    tuple val( meta ), path( fai )

    output:
    tuple val( meta ), file( "my.genome" ),     emit: dotgenome

    script:
    """
    cut -f1,2 $fai | sort -k2,2 -nr > my.genome
    """
}
