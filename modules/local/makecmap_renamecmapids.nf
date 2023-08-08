process MAKECMAP_RENAMECMAPIDS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::perl=5.26.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2' :
        'quay.io/biocontainers/perl:5.26.2' }"

    input:
    tuple val(meta), path(cmap)
    path keys

    output:
    tuple val(meta), path("*.cmap"), emit: renamedcmap
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"

    """
    rename_cmapids.pl -cmapfile $cmap -idx_key $keys $args > ${prefix}_EDITED.cmap

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/^.*perl //; s/Using.*\$//')
        rename_cmapids.pl: \$(rename_cmapids.pl -version)
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_EDITED.cmap

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/^.*perl //; s/Using.*\$//')
        rename_cmapids.pl: \$(rename_cmapids.pl -version)
    END_VERSIONS
    """
}
