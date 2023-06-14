process EXTRACT_BUSCOGENE {
    tag "$meta.id"
    label "process_low"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(fulltable)

    output:
    tuple val(meta), path("*.csv")  , emit: genefile
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args      ?: ''
    def prefix  = task.ext.prefix    ?: "${meta.id}"
    """
    get_busco_gene.py -o ${prefix}_buscogene.csv -i $fulltable

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/^.*python //; s/Using.*\$//')
        pandas: \$(echo \$(python -c "import pandas as pd; print(pd.__version__)"))
    END_VERSIONS
    """
}