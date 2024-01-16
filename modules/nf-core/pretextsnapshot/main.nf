process PRETEXTSNAPSHOT {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
            'docker.io/ubuntu:20.04' }"
    input:
    tuple val(meta), path(pretext_map)

    output:
    tuple val(meta), path('*.{jpeg,png,bmp}'), emit: image
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pretext_path = "${projectDir}/bin/PretextSnapshot/bin/PretextSnapshot"
    """
    ${pretext_path} \\
        $args \\
        --memory $task.memory \\
        --map $pretext_map \\
        --prefix $prefix \\
        --folder .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pretextsnapshot: \$(echo \$(pretext_path --version 2>&1) | sed 's/^.*PretextSnapshot Version //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pretextsnapshot: \$(echo \$(pretext_path --version 2>&1) | sed 's/^.*PretextSnapshot Version //' )
    END_VERSIONS
    """
}
