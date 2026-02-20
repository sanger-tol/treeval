process MAKECMAP_CMAP2BED {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(cmap)
    val enzyme

    output:
    tuple val(meta), path("*.bed"), emit: bedfile
    tuple val("${task.process}"), val('cmap2bed.py'), eval("python --version | sed 's/^.*python //; s/Using.*\$//'"), topic: versions, emit: versions_extractrepeat
    tuple val("${task.process}"), val('python'), eval("perl --version | sed -n 's/.*(v\\([0-9.]\\+\\)).*/\\1/p'"), topic: versions, emit: versions_perl

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    grep -v '#' $cmap > ${prefix}_${enzyme}_edited.cmap
    cmap2bed.py -t ${prefix}_${enzyme}_edited.cmap -z $enzyme | sort -k1,1 -k2,2n > ${enzyme}.bed
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${enzyme}.bed
    """
}
