process GET_LARGEST_SCAFF {

    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'ubuntu:20.04' }"

    input:
    tuple val( meta ), path( file )

    output:
    env largest_scaff   , emit: scaff_size
    path "versions.yml" , emit: versions

    shell:
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    $/
    largest_scaff=`cat "${file}" | head -n 1 - | cut -f2`

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    /$

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    largest_scaff=1000000

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    """
}
