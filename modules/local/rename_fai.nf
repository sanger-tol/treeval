process RENAME_FAI {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'docker.io/ubuntu:20.04' }"

    input:
    tuple val( meta ), path( fai )
    path( fasta )

    output:
    tuple val( meta ), file( "*fai" )   , emit: newfai
    path "versions.yml"                 , emit: versions

    script:
    def args    = task.ext.args      ?: ''
    def prefix  = task.ext.prefix    ?: "${meta.id}"
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    rename_fai.sh $fasta $fai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.fai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    """
}
