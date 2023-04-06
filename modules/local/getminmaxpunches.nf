process GETMINMAXPUNCHES{
    tag "${assembly_classT}"
    label "process_single"

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'ubuntu:20.04' }"

    input:
    tuple val(meta), path(bedfile)

    output:
    tuple val(meta), path ( '*zero.bed' ), emit: min
    tuple val(meta), path ( '*max.bed' ), emit: max
    path "versions.yml"          , emit: versions

    script:
    """
    coverage_punchlist.sh $bedfile

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(echo \$(bash --version | grep -Eo 'version [[:alnum:].]+' | sed 's/version //'))
    END_VERSIONS
    """
}
