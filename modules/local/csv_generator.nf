process CSV_GENERATOR {
    tag "${ch_org}"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::coreutils=9.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'ubuntu:20.04' }"

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
