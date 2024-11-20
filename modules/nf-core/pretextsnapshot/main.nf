process PRETEXTSNAPSHOT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :            'docker.io/ubuntu:20.04' }"

    input:
    tuple val(meta), path(pretext_map)

    output:
    tuple val(meta), path('*.{jpeg,png,bmp}'), emit: image
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION = "0.0.4"
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}."
    """
    PretextSnapshot \\
        $args \\
        --memory $task.memory \\
        --map $pretext_map \\
        --prefix $prefix \\
        --folder .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pretextsnapshot: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pretextsnapshot: $VERSION
    END_VERSIONS
    """
}
