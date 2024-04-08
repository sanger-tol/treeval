process SUBSAMPLE_BAM {
    tag "${meta.id}"
    label 'process_tiny'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(mergedbam)

    output:
    tuple val(meta), path('*.bam'), emit: subsampled_bam
    path "versions.yml",            emit: versions

    shell:
    def prefix = task.ext.prefix ?: "${meta.id}"
    '''
    percentage=`wc -c !{mergedbam} | cut -d$' ' -f1 | awk '{printf "%.2f\\n", 50000000000 / $0}'`

    samtools view -s $percentage -b !{mergedbam} > !{meta.id}_subsampled.bam

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    '''

    stub:
    """
    touch ${meta.id}_subsampled.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """
}
