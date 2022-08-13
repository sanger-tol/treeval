process CSV_GENERATOR {
    tag "${ch_org}"
    label 'process_small'

    input:
    val ch_dir
    val ch_org

    output:
    path "${ch_org}-data.csv",      emit: csv_path

    script:
    """
    cp "${ch_dir}${ch_org}-data.csv" "${ch_org}-data.csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvgenerator: BASH
    END_VERSIONS
    """
}
