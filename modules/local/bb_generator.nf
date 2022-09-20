// Software does not generate version data

process BB_GENERATOR {
    tag "${meta.id}-${meta.type}"
    label "process_medium"

    input:
    tuple val( meta ), path( blast_out ), path( dotas ), path( genome )

    output:
    tuple val( meta ), path( "$meta.id-$meta.type-BB.bb" ),        emit: bb_out

    script:
    """
    $projectDir/bin/bedToBigBed -as=$dotas -type=bed6+2 -extraIndex=name,geneSymbol $blast_out $genome $meta.id-$meta.type-BB.bb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedToBigBed: SOMETHING - UPDATE
    END_VERSIONS
    """
}
