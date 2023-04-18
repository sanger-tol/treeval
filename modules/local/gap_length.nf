process GAP_LENGTH {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'ubuntu:20.04' }"

    input:
    tuple val( meta ), path( file )

    output:
    tuple val( meta ), file( "*bed" ),      emit: bed

    script:
    """
    add_len2gap.sh $file > ${meta.id}_gaplen.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        add_len2gap: \$(add_len2gap.sh -v)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_gaplen.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        add_len2gap: \$(add_len2gap.sh -v)
    END_VERSIONS
    """
}
