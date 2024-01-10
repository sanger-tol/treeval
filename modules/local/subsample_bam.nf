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

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    $/
    bamsize=`wc -c "${mergedbam}" | cut -d$' ' -f1`
    threshold=50000000000
    percentage=`echo "scale=0;$threshold/$bamsize" | bc`

    if [[ $percentage -lt 1 ]]
    then
        samtools view -s $percentage -b ${mergedbam} > ${prefix}_subsampled.bam
    else
        mv ${mergedbam} ${prefix}_subsampled.bam
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    /$

    stub:
    """
    touch ${meta.id}_subsampled.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """
}
