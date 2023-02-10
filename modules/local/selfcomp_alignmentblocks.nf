process SELFCOMP_ALIGNMENTBLOCKS {
    tag "$meta.id"
    label 'process_single'
    publishDir "", enabled: false


    conda "anaconda::pandas=1.4.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'quay.io/biocontainers/pandas:1.4.3' }"

    input:
    tuple val(meta), path(bedfile)

    output:
    tuple val(meta), path("*.block"), emit: blockfile
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    build_alignment_block.py $args -i $bedfile

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/^.*python //; s/Using.*\$//')
        pandas: \$(echo \$(pandas: python -c "import pandas as pd; print(pd.__version__)")
        build_alignment_block.py: \$(build_alignment_block.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}

