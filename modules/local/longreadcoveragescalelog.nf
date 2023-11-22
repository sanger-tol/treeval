process LONGREADCOVERAGESCALELOG {
    tag "${meta.id}"
    label "process_single"

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(bedfile)

    output:
    tuple val(meta), path("*.bed")  , emit: bed
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args     ?: ''
    def prefix  = task.ext.prefix   ?: "log"
    """
    longread_cov_log.py -i $bedfile > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/^.*python //; s/Using.*\$//')
        longread_cov_log.py: \$(longread_cov_log.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix   ?: "log2"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/^.*python //; s/Using.*\$//')
        longread_cov_log2.py: \$(longread_cov_log2.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
