process MAKECMAP_RENAMECMAPIDS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::perl=5.26.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2' :
        'biocontainers/perl:5.26.2' }"

    input:
    tuple val(meta), path(cmap)
    path keys

    output:
    tuple val(meta), path("*.cmap"), emit: renamedcmap
    tuple val("${task.process}"), val('rename_cmapids.pl'), eval("rename_cmapids.pl -version"), topic: versions, emit: versions_renamedcmap
    tuple val("${task.process}"), val('perl'), eval("perl --version | sed -n 's/.*(v\\([0-9.]\\+\\)).*/\\1/p'"), topic: versions, emit: versions_perl

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"

    """
    rename_cmapids.pl -cmapfile $cmap -idx_key $keys $args > ${prefix}_EDITED.cmap
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_EDITED.cmap
    """
}
