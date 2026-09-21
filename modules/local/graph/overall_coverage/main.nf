process GRAPH_OVERALL_COVERAGE {
    tag "${meta.id}"
    label "process_single"

    conda "conda-forge::perl=5.26.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2' :
        'biocontainers/perl:5.26.2' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.part") , emit: part
    tuple val("${task.process}"), val('graph_overall_coverage.pl'), eval("graph_overall_coverage.pl --version"), topic: versions, emit: versions_graph_overall_coverage
    tuple val("${task.process}"), val('perl'), eval("perl --version | sed -n 's/.*(v\\([0-9.]\\+\\)).*/\\1/p'"), topic: versions, emit: versions_perl

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    graph_overall_coverage.pl $bed > ${prefix}.part
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.part
    """
}
