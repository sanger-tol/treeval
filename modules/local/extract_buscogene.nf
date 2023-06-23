process EXTRACT_BUSCOGENE {
    tag "$meta.id"
    label "process_low"

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'ubuntu:20.04' }"


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
    get_busco_gene.sh $fulltable > ${prefix}_buscogene.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        get_busco_gene: \$(get_busco_gene.sh -v)
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix    ?: "${meta.id}"
    """
    touch ${prefix}_buscogene.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        get_busco_gene: \$(get_busco_gene.sh -v)
    END_VERSIONS
    """
}
