process GENERATE_GENOME_FILE {
    tag "${meta.id}"
    label "process_low"

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'ubuntu:20.04' }"

    input:
    tuple val( meta ), path( fai )

    output:
    tuple val( meta ), file( "my.genome" )      , emit: dotgenome
    path "versions.yml"                         , emit: versions

    script:
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    cut -f1,2 $fai | sort -k2,2 -nr > my.genome

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    """

    stub:
    def VERSION = "9.1"
    """
    cut -f1,2 $fai | sort -k2,2 -nr > my.genome

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    """
}
