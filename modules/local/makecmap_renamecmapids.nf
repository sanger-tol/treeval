process MAKECMAP_RENAMECMAPIDS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::perl=5.26.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2' :
        'perl:5.26.2' }"

    input:
    tuple val(meta), path(cmap)
    path keys

    output:
    tuple val(meta), path("*.cmap"), emit: renamedcmap
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    rename_cmapids.pl -cmapfile $cmap -idx_key $keys $args > ${prefix}_EDITED.cmap

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/^.*perl //; s/Using.*\$//')
    END_VERSIONS
    """
}
