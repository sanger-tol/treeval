process MINIMAP_SAMTOOLS {
    tag "${meta.id}"
    label "process_medium"

    def version = '2.22_1.12'

    container "niemasd/minimap2_samtools:${version}"

    input:
    tuple val( ref_meta ), file( ref )
    tuple val( meta ), file( nuc_file )

    output:
    tuple val( meta ), file( "${nuc_file}.bam" ),     emit: partial_alignment
    path "versions.yml",                              emit: versions

    script:
    def intron_size = task.ext.args ?: '50k'
    """
    minimap2 -ax splice $ref $nuc_file -G $intron_size | samtools view -Sb -T $ref - > ${nuc_file}.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MINIMAP_SAMTOOLS: $version
    END_VERSIONS
    """
}

//USE THE OFFICIAL MINIMAP MODULE