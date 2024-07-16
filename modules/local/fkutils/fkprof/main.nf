process FKUTILS_FKPROF {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'docker.io/ubuntu:20.04' }"

    input:
    tuple val(meta), path(reference)
    tuple val(meta2), path(ktab)

    output:
    tuple val(meta), file("*bed"),  emit: bed
    path "versions.yml",            emit: versions

    script:
    def args    = task.ext.args     ?: ""
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    def VERSION = "2f7b2583092fafb0a5abc654b5857642" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    ${projectDir}/bin/FKprof \\
        $args \\
        -t$task.cpus \\
        ${reference} \\
        ${ktab}/*ktab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fkprof:    $VERSION
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    def kmer    = 31
    def VERSION = "2f7b2583092fafb0a5abc654b5857642"  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}_${kmer}_read.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fkprof:    $VERSION
    END_VERSIONS
    """
}
