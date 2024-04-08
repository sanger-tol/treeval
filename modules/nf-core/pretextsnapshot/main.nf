process PRETEXTSNAPSHOT {
    tag "$meta.id"
    label 'process_single'

    container "quay.io/sanger-tol/pretext:0.0.2-yy5-c3"

    input:
    tuple val(meta), path(pretext_map)

    output:
    tuple val(meta), path('*.{jpeg,png,bmp}'), emit: image
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION         = "0.0.4"
    def args            = task.ext.args ?: ''
    def prefix          = task.ext.prefix ?: "${meta.id}"

    """
    PretextSnapshot \\
        $args \\
        --memory $task.memory \\
        --map $pretext_map \\
        --prefix $prefix \\
        --folder .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PretextSnapshot: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PretextSnapshot: $VERSION
    END_VERSIONS
    """
}
