process REFORMAT_INTERSECT {
    tag "${meta.id}"
    label "process_low"

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'ubuntu:20.04' }"

    input:
    tuple val( meta ), path( file )

    output:
    tuple val( meta ), file( "*.bed" ), emit: bed

    shell:
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    cat $file | reformat.sh - > ${meta.id}_fmt_INTERSECT.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        reformat:   \$(reformat.sh -v)
        coreutils:  $VERSION
    END_VERSIONS
    """

    stub:
    def VERSION = "9.1"
    """
    touch ${meta.id}_fmt_INTERSECT.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        reformat:   \$(reformat.sh -v)
        coreutils:  $VERSION
    END_VERSIONS
    """
}
