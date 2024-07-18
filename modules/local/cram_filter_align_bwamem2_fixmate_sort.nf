process CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT {
    tag "$meta.id"
    label "process_high"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1a6fe65bd6674daba65066aa796ed8f5e8b4687b:688e175eb0db54de17822ba7810cc9e20fa06dd5-0' :
        'biocontainers/mulled-v2-1a6fe65bd6674daba65066aa796ed8f5e8b4687b:688e175eb0db54de17822ba7810cc9e20fa06dd5-0' }"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT module does not support Conda. Please use Docker / Singularity instead."
    }

    input:
    tuple val(meta), path(cramfile), path(cramindex), val(from), val(to), val(base), val(chunkid), val(rglines), val(bwaprefix), path(reference)

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
    // Please be aware one of the tools here required mem = 28 * reference size!!!
    """
    cram_filter -n ${from}-${to} ${cramfile} - | \\
        samtools fastq ${args1} | \\
        bwa-mem2 mem -p ${bwaprefix} -t${task.cpus} -5SPCp -H'${rglines}' - | \\
        samtools fixmate ${args3} - - | \\
        samtools sort ${args4} -@${task.cpus} -T ${base}_${chunkid}_sort_tmp -o ${prefix}_${base}_${chunkid}_mem.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
        bwa-mem2: \$(bwa-mem2 --version | sed 's/bwa-mem2 //g')
    END_VERSIONS
    """
    // temp removal staden_io_lib: \$(echo \$(staden_io_lib --version 2>&1) | sed 's/^.*staden_io_lib //; s/Using.*\$//') CAUSES ERROR

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def base    = "45022_3#2"
    def chunkid = "1"
    """
    touch ${prefix}_${base}_${chunkid}_mem.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
    END_VERSIONS
    """
}
