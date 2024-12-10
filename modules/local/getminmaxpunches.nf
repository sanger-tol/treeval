process GETMINMAXPUNCHES{
    tag "${meta.id}"
    label "process_single"

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'docker.io/ubuntu:20.04' }"

    input:
    tuple val(meta), path(bedfile)

    output:
    tuple val(meta), path ( '*zero.bed' )   , optional: true    , emit: min
    tuple val(meta), path ( '*max.bed' )    , optional: true    , emit: max
    path "versions.yml"                     , emit: versions

    shell:
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    $/
    cat "${bedfile}" \
    | awk '{ if ($4 == 0) {print $0 >> "zero.bed" } else if ($4 > 1000) {print $0 >> "max.bed"}}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    /$

    stub:
    def VERSION = "9.1"  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch max.bed
    touch min.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    """
}
