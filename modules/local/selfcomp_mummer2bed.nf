process SELFCOMP_MUMMER2BED {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'python:3.9' }"

    input:
    tuple val(meta), path(mummerfile)
    val (motiflen)

    output:
    tuple val(meta), path("*.bed"), emit: bedfile
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mummer2bed.py $args -i $mummerfile -l $motiflen > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/^.*python //; s/Using.*\$//')
        mummer2bed.py: \$(mummer2bed.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
