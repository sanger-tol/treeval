process CONCATMUMMER {
    tag "${meta.id}.mummer"
    label "process_medium"

    conda (params.enable_conda ? "conda-forge::coreutils=9.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
            'ubuntu:20.04' }"

    input:
    tuple val(meta), path(coords)

    output:
    tuple val(meta), path("*.mummer"), emit: mummer
    path "versions.yml", emit: versions

    script:
    """
    cat $coords > ${meta.id}.mummer

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ubuntu: \$(ubuntu --version | sed 's/Ubuntu //g')
    END_VERSIONS
    """
}
