process CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT {
    tag "$meta.id"
    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-50d89b457e04ed90fa0cbf8ebc3ae1b9ffbc836b:caf993da1689e8d42f5e4c113ffc9ef81d26df96-0' :
        'quay.io/biocontainers/mulled-v2-50d89b457e04ed90fa0cbf8ebc3ae1b9ffbc836b:caf993da1689e8d42f5e4c113ffc9ef81d26df96-0' }"

    input:
    tuple val(meta), path(cramfile), path(cramindex), path(reffasta), val(from), val(to), val(base), val(chunkid), val(rglines), val(bwaprefix)

    output:
    tuple val(meta), path("*.bam"), emit: mappedbam
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cram_filter -n ${from}-${to} ${cramfile} - | samtools fastq -F0xB00 -nt - | bwa-mem2 mem -p ${bwaprefix} -t16 -5SPCp -H'${rglines}' ${reffasta} - | samtools fixmate -mpu - - | samtools sort --write-index -l1 -@16 -T ${base}_${chunkid}_sort_tmp -o ${prefix}_${base}_${chunkid}_mem.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        staden_io_lib: \$(echo \$(staden_io_lib --version 2>&1) | sed 's/^.*staden_io_lib //; s/Using.*\$//')
        samtools: \$(samtools --version | sed 's/samtools //g')
        bwa-mem2: \$(bwa-mem2 --version | sed 's/bwa-mem2 //g')  
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${base}_${chunkid}_mem.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        staden_io_lib: \$(echo \$(staden_io_lib --version 2>&1) | sed 's/^.*staden_io_lib //; s/Using.*\$//')
        samtools: \$(samtools --version | sed 's/samtools //g')
        bwa-mem2: \$(bwa-mem2 --version | sed 's/bwa-mem2 //g')
    END_VERSIONS
    """
}