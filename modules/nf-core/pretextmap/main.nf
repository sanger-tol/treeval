
process PRETEXTMAP {
    tag "$meta.id"
    label 'process_single'

    container "quay.io/sanger-tol/pretext:0.0.2-yy5-c3"

    input:
    tuple val(meta), path(input)
    path fasta

    output:
    tuple val(meta), path("*.pretext"), emit: pretext
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION         = "0.1.9"
    def args            = task.ext.args ?: ''
    def prefix          = task.ext.prefix ?: "${meta.id}"
    def reference       = fasta ? "--reference ${fasta}" : ""

    """
    if [[ $input == *.pairs.gz ]]; then
        zcat $input | PretextMap \\
            $args \\
            -o ${prefix}.pretext
    else
        samtools \\
            view \\
            $reference \\
            -h \\
            $input | PretextMap \\
            $args \\
            -o ${prefix}.pretext
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PretextMap: $VERSION
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """

    stub:
    def VERSION         = "0.1.9"
    def prefix          = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pretext

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PretextMap: $VERSION
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
