process CUSTOM_GETCHROMSIZES {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
    tuple   val(meta), path(fasta)
    val     suffix

    output:
    tuple val(meta), path ("*.${suffix}")   , emit: sizes
    tuple val(meta), path ("*.fa")          , emit: fasta
    tuple val(meta), path ("*.fai")         , emit: fai
    tuple val(meta), path ("*.gzi")         , emit: gzi, optional: true
    path  "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    samtools faidx $fasta -o ${prefix}.fa.fai
    cut -f 1,2 ${prefix}.fa.fai > ${prefix}.${suffix}

    if [[ "${fasta}" != "${prefix}-ref.fa" ]]; then
        mv ${fasta} ${prefix}-ref.fa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getchromsizes: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${prefix}.fai
    touch ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getchromsizes: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
