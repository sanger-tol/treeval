process GAP_LENGTH {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'ubuntu:20.04' }"

    input:
    tuple val( meta ), path( file )

    output:
    tuple val( meta ), file( "*bed" ),      emit: bed
    path "versions.yml"              ,      emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    add_len2gap.sh $file > pretext_${prefix}_gap.bedgraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        add_len2gap: \$(add_len2gap.sh -v)
        coreutils: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "9.1"
    """
    touch ${prefix}_gap.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        add_len2gap: \$(add_len2gap.sh -v)
        coreutils: $VERSION
    END_VERSIONS
    """
}
