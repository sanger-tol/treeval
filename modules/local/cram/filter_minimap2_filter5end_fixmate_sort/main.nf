process CRAM_FILTER_MINIMAP2_FILTER5END_FIXMATE_SORT {
    tag "${meta.id}"
    label "process_high"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1a6fe65bd6674daba65066aa796ed8f5e8b4687b:688e175eb0db54de17822ba7810cc9e20fa06dd5-0' :
        'biocontainers/mulled-v2-1a6fe65bd6674daba65066aa796ed8f5e8b4687b:688e175eb0db54de17822ba7810cc9e20fa06dd5-0' }"

    input:
    tuple val(meta), path(cramfile), path(cramindex), val(from), val(to), val(base), val(chunkid), val(rglines), val(ref), path(reference)
    path(grep_pg_program_file)
    path(filter_five_end_program_file)
    path(awk_filter_reads_program_file)

    output:
    tuple val(meta), path("*.bam"), emit: mappedbam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args1 = task.ext.args1 ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.15" // Staden_io versions break the pipeline
    """
    cram_filter -n ${from}-${to} ${cramfile} - | \\
        samtools fastq ${args} - |  \\
        minimap2 -t${task.cpus} -R '${rglines}' ${args1} ${ref} - |  \\
        ${grep_pg_program_file} |  \\
        perl ${filter_five_end_program_file} |  \\
        ${awk_filter_reads_program_file} |  \\
        samtools fixmate ${args2} - - | \\
        samtools sort ${args3} -@${task.cpus} -T ${base}_${chunkid}_sort_tmp -o ${prefix}_${base}_${chunkid}_mm.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
        minimap2: \$(minimap2 --version | sed 's/minimap2 //g')
        staden_io: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.15" // Staden_io versions break the pipeline
    """
    touch ${prefix}_${base}_${chunkid}_mm.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
        minimap2: \$(echo \$(minimap2 version 2>&1) | sed 's/.* //')
        staden_io: $VERSION
    END_VERSIONS
    """
}
