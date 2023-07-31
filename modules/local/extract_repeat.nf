process EXTRACT_REPEAT {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::perl=5.26.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2' :
        'quay.io/biocontainers/perl:5.26.2' }"

    input:
    tuple val( meta ), path( file )

    output:
    tuple val( meta ), path( "*.bed" )  , emit: bed
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def VERSION     = "1.0" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    extract_repeat.pl $file > ${prefix}_repeats.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/^.*perl //; s/Using.*\$//')
        extract_repeat.pl: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def VERSION     = "1.0" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}_repeats.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/^.*perl //; s/Using.*\$//')
        extract_repeat.pl: $VERSION
    END_VERSIONS
    """
}
