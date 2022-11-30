process CSV_GENERATOR {
    tag "${ch_org}"
    label 'process_low'

    input:
    val ch_org
    val data_dir
    val classT

    output:
    path "${ch_org}-data.csv",      emit: csv_path

    script:
    def csv_loc = "/csv_data/"
    """
    cp "${data_dir}${classT}${csv_loc}${ch_org}-data.csv" "${ch_org}-data.csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvgenerator: BASH
    END_VERSIONS
    """
}
