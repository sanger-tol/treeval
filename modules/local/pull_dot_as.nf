process PULL_DOT_AS {
    tag "${meta.id}-${meta.type}"
    label 'process_small'

    input:
    tuple val( meta ), file( tsv_file )

    output:
    tuple val( meta ), path( '*as' ),       emit: dotas

    script:
    """
    cp $projectDir/assets/gene_alignment/assm_${meta.type}.as .
    """
}
