process SELFCOMP_MERGEBLOCKS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bedtools"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools%3A2.26.0gx--0' :
        'quay.io/biocontainers/bedtools' }"

    input:
    tuple val(meta), path(blockfile)

    output:
    tuple val(meta), path("*.bed"), emit: mergedblocks
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bedtools merge -i $blockfile -c 4,7,8 -o collapse,min,max -delim "|" -d 100000 > ${meta.id}_merged.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
