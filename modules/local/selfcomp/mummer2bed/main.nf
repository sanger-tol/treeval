process SELFCOMP_MUMMER2BED {
    tag "${meta.id}"
    label "process_medium"

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(mummerfile)

    output:
    tuple val(meta), path("*.bed"), emit: bedfile
    tuple val("${task.process}"), val('mummer2bed.py'), eval("mummer2bed.py --version | cut -d' ' -f2"), topic: versions, emit: versions_mummer2bed
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/^Python //'"), topic: versions, emit: versions_python


    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"

    """
    mummer2bed.py $args -i $mummerfile -l 0 > ${prefix}.bed
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    """
}
