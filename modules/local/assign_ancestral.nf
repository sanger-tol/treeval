process ASSIGN_ANCESTRAL {
    tag "$meta.id"
    label "process_low"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(comp_location)
    path(fulltable)

    output:
    tuple val(meta), path("*bed")  , emit: assigned_bed
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args     ?: ''
    def prefix  = task.ext.prefix   ?: "${meta.id}"

    """
    assign_anc.py -l $comp_location -f $fulltable -c ${prefix}_assigned.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/^.*python //; s/Using.*\$//')
        pandas: \$(echo \$(pandas: python -c "import pandas as pd; print(pd.__version__)"))
        assign_ancestral.py: \$(assign_ancestral.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}