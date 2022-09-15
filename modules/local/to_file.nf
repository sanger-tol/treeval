process TO_FILE {
    tag "${sample_id}"
    label "process_small"

    input:
    val ( sample_id )
    val ( r_file )

    output:
    tuple val( sample_id ), file( '*fa' ),    emit: file_path

    script:
    """
    cp "${r_file}" "./${sample_id}-genome.fa"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        TO_FILE: BASH
    END_VERSIONS
    """
}
