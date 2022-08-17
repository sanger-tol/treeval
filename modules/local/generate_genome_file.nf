process GENERATE_GENOME_FILE {
<<<<<<< HEAD
    tag "${meta.id}"
=======
    tag "${meta}"
>>>>>>> 8740473 (Adding GENERATE_GENOME subworkflow to main)
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