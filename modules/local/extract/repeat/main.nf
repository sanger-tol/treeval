process EXTRACT_REPEAT {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::perl=5.26.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2' :
        'biocontainers/perl:5.26.2' }"

    input:
    tuple val(meta), path(file)

    output:
    tuple val( meta ), path( "*.bed" )  , emit: bed
    tuple val("${task.process}"), val('extract_repeat.pl'), eval("echo '1.0.0'"), topic: versions, emit: versions_extractrepeat
    tuple val("${task.process}"), val('perl'), eval("perl --version | sed -n 's/.*(v\\([0-9.]\\+\\)).*/\\1/p'"), topic: versions, emit: versions_perl

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    extract_repeat.pl $file > ${prefix}_repeats.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_repeats.bed
    """
}
