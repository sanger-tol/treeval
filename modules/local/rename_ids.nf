process RENAME_IDS {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'docker.io/ubuntu:20.04' }"

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), file("*bed")   , emit: bed
    path "versions.yml"             , emit: versions

    shell:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    $/
    cat "${file}" \
    | sed 's/\./0/g' > ${prefix}_renamed.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils:     $VERSION
    END_VERSIONS
    /$

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "9.1"  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}_renamed.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rename_ids:   \$(rename_ids.sh -v)
        coreutils:      $VERSION
    END_VERSIONS
    """

}
