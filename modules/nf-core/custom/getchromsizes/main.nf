// Forked from the nf-core module to:
//  1. allow selecting a different extension for the `sizes` channel
//  2. force all output files to be named according to the prefix
//  3. rename the input fasta file too and output it so that it can be "published"
process CUSTOM_GETCHROMSIZES {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
    tuple   val(meta), path(fasta, stageAs: 'input/*')
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
    ln -s ${fasta} ${prefix}.fa
    samtools faidx ${prefix}.fa -o ${prefix}.fa.fai
    cut -f 1,2 ${prefix}.fa.fai > ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getchromsizes: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    ln -s ${fasta} ${prefix}.fa
    touch ${prefix}.fa.fai
    touch ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getchromsizes: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
