process CRAM_FILTER_MINIMAP2_FILTER5END_FIXMATE_SORT {
    tag "$meta.id"
    label "process_high"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-50d89b457e04ed90fa0cbf8ebc3ae1b9ffbc836b:caf993da1689e8d42f5e4c113ffc9ef81d26df96-0' :
        'biocontainers/mulled-v2-50d89b457e04ed90fa0cbf8ebc3ae1b9ffbc836b:caf993da1689e8d42f5e4c113ffc9ef81d26df96-0' }"

    input:
    tuple val(meta), path(cramfile), path(cramindex), val(from), val(to), val(base), val(chunkid), val(rglines), val(bwaprefix)
    path ref

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
    def args4 = task.ext.args4 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cram_filter -n ${from}-${to} ${cramfile} - | \\
        samtools fastq - |  \\
        minimap2 -t${task.cpus} '${rglines}' -axsr ${ref} - |  \\
        grep -v "^\@PG" |  \\
        awk '{if(\$1 ~ /^\@/) {print(\$0)} else {if(and(\$2,64)>0) {print(1\$0)} else {print(2\$0)}}}' |  \\
        filter_five_end.pl |  \\
        awk 'BEGIN{OFS="\\t"}{if(\$1 ~ /^\@/) {print(\$0)} else {\$2=and(\$2,compl(2048)); print(substr(\$0,2))}}' |  \\
        samtools fixmate ${args3} - - | \\
        samtools sort ${args4} -@${task.cpus} -T ${base}_${chunkid}_sort_tmp -o ${prefix}_${base}_${chunkid}_mm.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
        minimap2: \$(minimap2 --version | sed 's/minimap2 //g')
    END_VERSIONS
    """
    // temp removal staden_io_lib: \$(echo \$(staden_io_lib --version 2>&1) | sed 's/^.*staden_io_lib //; s/Using.*\$//') CAUSES ERROR

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def base    = "45022_3#2"
    def chunkid = "1"
    """
    touch ${prefix}_${base}_${chunkid}_mm.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
        minimap2: \$(echo \$(minimap2 version 2>&1) | sed 's/.* //')
    END_VERSIONS
    """
}
