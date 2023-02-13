process SELFCOMP_MAPIDS {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'quay.io/biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(bed)
    path(agp)

    output:
    tuple val(meta), path("*.bed"), emit: bedfile
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mapids.py -i $bed -r $agp > ${prefix}_mapped.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/^.*python //; s/Using.*\$//')
        mapids.py: \$(mapids.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
