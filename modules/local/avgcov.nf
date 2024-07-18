process AVGCOV {
    tag "${meta.id}"
    label "process_single"

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'docker.io/ubuntu:20.04' }"

    input:
    tuple val(meta), path(bedfile)
    path genomefile

    output:
    tuple val(meta), path("*.bed")  , emit: avgbed
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args     ?: ''
    def prefix  = task.ext.prefix   ?: "avgcov"
    """
    get_avgcov.sh $bedfile $genomefile ${prefix}.bed $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        get_avgcov: \$(get_avgcov.sh -v)
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix   ?: "avgcov"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        get_avgcov: \$(get_avgcov.sh -v)
    END_VERSIONS
    """
}
