process SELFCOMP_ALIGNMENTBLOCKS {
    tag "$meta.id"
    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-548f120dc8914d802c46e110ec27751bc1c5a414:8770fa59aa0ae8b50cbf444255b91c201c883685-0' :
        'quay.io/biocontainers/mulled-v2-548f120dc8914d802c46e110ec27751bc1c5a414:8770fa59aa0ae8b50cbf444255b91c201c883685-0' }"

    input:
    tuple val(meta), path(bedfile)

    output:
    tuple val(meta), path("*.block"), emit: blockfile
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args     ?: ''
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    """
    build_alignment_block.py $args -i $bedfile -o ${prefix}_chained.block

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/^.*python //; s/Using.*\$//')
        pandas: \$(echo \$(pandas: python -c "import pandas as pd; print(pd.__version__)"))
        pybedtools: \$(echo \$(pybedtools: python -c "import pybedtools as pb; print(pb.__version__)"))
        build_alignment_block.py: \$(build_alignment_block.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    """
    touch ${prefix}_chained.block

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/^.*python //; s/Using.*\$//')
        pandas: \$(echo \$(pandas: python -c "import pandas as pd; print(pd.__version__)"))
        pybedtools: \$(echo \$(pybedtools: python -c "import pybedtools as pb; print(pb.__version__)"))
        build_alignment_block.py: \$(build_alignment_block.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
