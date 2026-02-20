process FIND_HALF_COVERAGE {
    tag "${meta.id}"
    label "process_single"

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(bedfile)
    path(my_genome)
    path(depthgraph)

    output:
    tuple val(meta), path("*.bed")  , emit: bed
    tuple val("${task.process}"), val('findHalfcoverage.py'), eval("findHalfcoverage.py --version | sed 's/^.*.py //'"), topic: versions, emit: versions_cmap2bed
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/^Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix   ?: "halfcoverage"
    """
    findHalfcoverage.py -c $bedfile -m $my_genome -d $depthgraph > ${prefix}.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "halfcoverage"
    """
    touch ${prefix}.bed
    """
}
