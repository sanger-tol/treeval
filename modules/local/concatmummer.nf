process CONCATMUMMER {
    tag "${meta.id}.mummer"
    label "process_medium"

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
            'ubuntu:20.04' }"

    input:
    tuple val(meta), path(coords)

    output:
    tuple val(meta), path("*.mummer")   , emit: mummer
    path "versions.yml"                 , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    cat $coords > ${prefix}.mummer

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ubuntu: \$(ubuntu --version | sed 's/Ubuntu //g')
        coreutils: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "9.1"
    """
    touch ${prefix}.mummer

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ubuntu: \$(ubuntu --version | sed 's/Ubuntu //g')
        coreutils: $VERSION
    END_VERSIONS
    """
}
