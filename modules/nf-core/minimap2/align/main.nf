process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    // Note: the versions here need to match the versions used in the mulled container below and minimap2/index
    conda "bioconda::minimap2=2.24 bioconda::samtools=1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' }"

    input:
    tuple val(meta), path(reads)
    path reference
    val bam_format
    val cigar_paf_format
    val cigar_bam
    val bed_format

    output:
    tuple val(meta), path("*.paf"), optional: true, emit: paf
    tuple val(meta), path("*.bam"), optional: true, emit: bam
    tuple val(meta), path("*.bed"), optional: true, emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bam_output = reference.size() > 2.5e9 && bam_format ? "-a | samtools view -b -T ${reference} - > ${prefix}.bam" : reference.size() < 2.5e9 && bam_format ? "-a | samtools view -@ ${task.cpus} -b -h -o ${prefix}.bam" : bed_format ? "| paftools.js splice2bed - > ${prefix}.bed " : "-o ${prefix}.paf"
    def cigar_paf = cigar_paf_format && !bam_format ? "-c" : ''
    def set_cigar_bam = cigar_bam && bam_format ? "-L" : ''

    """
    minimap2 \\
        $args \\
        -t $task.cpus \\
        "${reference ?: reads}" \\
        "$reads" \\
        $cigar_paf \\
        $set_cigar_bam \\
        $bam_output


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bam_output = reference.size() > 2.5e9 && bam_format ? "-a | samtools view -b -T ${reference} - > ${prefix}.bam" : reference.size() < 2.5e9 && bam_format ? "-a | samtools view -@ ${task.cpus} -b -h -o ${prefix}.bam" : "-o ${prefix}.paf"
    def cigar_paf = cigar_paf_format && !bam_format ? "-c" : ''
    def set_cigar_bam = cigar_bam && bam_format ? "-L" : ''
    def extension = bam_format ? "bam" : bed_format ? "bed" : "paf"
    """
    touch ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """

}

