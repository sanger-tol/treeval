process MAKECMAP_CMAP2BED {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(cmap)
    val enzyme

    output:
    tuple val(meta), path("*.bed"), emit: bedfile
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    grep -v '#' $cmap > ${prefix}_${enzyme}_edited.cmap
    cmap2bed.py -t ${prefix}_${enzyme}_edited.cmap -z $enzyme | sort -k1,1 -k2,2n > ${enzyme}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/^.*python //; s/Using.*\$//')
        cmap2bed.py: \$(cmap2bed.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${enzyme}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/^.*python //; s/Using.*\$//')
        cmap2bed.py: \$(cmap2bed.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
