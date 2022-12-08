// Software does not generate version data

process BB_GENERATOR {
    tag "${meta.id}-${meta.type}"
    label "process_medium"

    conda (params.enable_conda ? "conda-forge::coreutils=9.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'ubuntu:20.04' }"

    input:
    tuple val( meta ), path( blast_out ), path( dotas ), path( genome )

    output:
    tuple val( meta ), path( "$meta.id-$meta.type-BB.bb" ),        emit: bb_out

    script:
    """
    $projectDir/bin/bedToBigBed -as=$dotas -type=bed6+2 -extraIndex=name,geneSymbol $blast_out $genome $meta.id-$meta.type-BB.bb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedToBigBed: BASH
    END_VERSIONS
    """
}
