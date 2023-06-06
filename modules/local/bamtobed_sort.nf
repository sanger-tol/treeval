process BAMTOBED_SORT {
    tag "$meta.id"
    label "process_high"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:7f9016e4e52c90f2ea56773311205fd2d61479d9-0' :
        'quay.io/biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:7f9016e4e52c90f2ea56773311205fd2d61479d9-0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bed"), emit: sorted_bed

    script:
    def prefix = args.ext.prefix ?: "${meta.id}"
    def thing = '--parallel=8 -S50G'
    """
    samtools view -@4 -u -F0x400 ${bam} | bamToBed | sort -k4 > ${prefix}_merged_sorted.bed
    """

    stub:
    def prefix = args.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_merged_sorted.bed
    """
}
