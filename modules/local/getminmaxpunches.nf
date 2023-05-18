process GETMINMAXPUNCHES{
    tag "${meta.id}"
    label "process_single"

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'ubuntu:20.04' }"

    input:
    tuple val(meta), path(bedfile)

    output:
    tuple val(meta), path ( '*zero.bed' )   , optional: true    , emit: min
    tuple val(meta), path ( '*max.bed' )    , optional: true    , emit: max
    path "versions.yml"                     , emit: versions

    script:
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    coverage_punchlist.sh $bedfile

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(echo \$(bash --version | grep -Eo 'version [[:alnum:].]+' | sed 's/version //'))
        coreutils: $VERSION
    END_VERSIONS
    """

    stub:
    def VERSION = "9.1"
    """
    touch max.bed
    touch min.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(echo \$(bash --version | grep -Eo 'version [[:alnum:].]+' | sed 's/version //'))
        coreutils: $VERSION
    END_VERSIONS
    """
}
