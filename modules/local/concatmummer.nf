process CONCATMUMMER {
    tag "${meta.id}.mummer"
    label "process_medium"

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the mummer2bed process. Please use docker or singularity containers."
    }
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
