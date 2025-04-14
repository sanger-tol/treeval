process EXTRACT_ANCESTRAL {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(fulltable)
    path(ancestraltable)

    output:
    tuple val(meta), path("*buscopainter_complete_location.tsv")  , emit: comp_location
    path("*buscopainter_duplicated_location.tsv")                 , emit: dup_location
    path("*summary.tsv")                                          , emit: summary
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args     ?: ''
    def prefix  = task.ext.prefix   ?: "${meta.id}"

    """
    buscopainter.py -r $ancestraltable -q $fulltable

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/^.*python //; s/Using.*\$//')
        buscopainter.py: \$(buscopainter.py -v)
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    """
    touch ${prefix}_buscopainter_complete_location.tsv
    touch ${prefox}_buscopainter_duplicated_location.tsv
    touch ${prefix}_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/^.*python //; s/Using.*\$//')
        buscopainter.py: \$(buscopainter.py -v)
    END_VERSIONS
    """
}
