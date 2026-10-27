process MAKECMAP_FA2CMAPMULTICOLOR {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::perl=5.26.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2' :
        'biocontainers/perl:5.26.2' }"

    input:
    tuple val(meta), path(fasta)
    val enzyme

    output:
    tuple val(meta), path("*.cmap"), emit: cmap
    path("*key.txt")               , emit: cmapkey
    tuple val("${task.process}"), val('fa2cmap_multi_color.pl'), eval("fa2cmap_multi_color.pl -v | sed 's/fa2cmap_multi_color.pl //g'"), topic: versions, emit: versions_fa2cmap_multi_color
    tuple val("${task.process}"), val('perl'), eval("perl --version | sed -n 's/.*(v\\([0-9.]\\+\\)).*/\\1/p'"), topic: versions, emit: versions_perl

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    """
    fa2cmap_multi_color.pl -i $fasta -e $enzyme 1 $args
    """

    stub:
    """
    touch ${enzyme}_key.cmap
    """
}
