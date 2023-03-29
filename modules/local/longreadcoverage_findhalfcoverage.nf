process LONGREADCOVERAGE_FINDHALFCOVERAGE {
    tag "$meta.id"

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'quay.io/biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(bedfile)
    path(my_genome)
    path(depthgraph)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    findHalfcoverage.py -c $bedfile -m $my_genome -d $depthgraph > ${prefix}_halfdepth.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/^.*python //; s/Using.*\$//')
        findHalfcoverage.py: \$(findHalfcoverage.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}