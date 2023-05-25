process EXTRACT_ANCESTRAL {
    tag "$meta.id"

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'quay.io/biocontainers/python:3.9' }"

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
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    sed -e '1,2d' $fulltable | sed 's/# //g' > edited_fulltable.tsv
    buscopainter.py -r $ancestraltable -q edited_fulltable.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/^.*python //; s/Using.*\$//')
        buscopainter.py: \$(buscopainter.py -v)
    END_VERSIONS
    """
}
