
process PRETEXTMAP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/sanger-tol/pretext:0.0.8-yy5-c1"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.pretext")  , emit: pretext
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args     ?: ''
    def prefix      = task.ext.prefix   ?: "${meta.id}"
    def reference   = fasta             ? "--reference ${fasta}" : ""

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
        PretextGraph: \$(PretextGraph | grep "Version" | sed 's/Pretext* Version //;')
        PretextMap: \$(PretextMap | grep "Version" | sed 's/PretextMap Version//;')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """

    stub:
    def prefix      = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pretext

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PretextGraph: \$(PretextGraph | grep "Version" | sed 's/Pretext* Version //;')
        PretextMap: \$(PretextMap | grep "Version" | sed 's/PretextMap Version//;')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
