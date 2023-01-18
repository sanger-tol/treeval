process CSV_GENERATOR {
    tag "${ch_org}"
    label 'process_low'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'ubuntu:20.04' }"

    input:
    val ch_org
    val data_dir
    val classT

    output:
    val "${data_dir}${classT}/csv_data/${ch_org}-data.csv",      emit: csv_path

    script:
    def csv_loc = "/csv_data/"
    """
    echo "Path generated = ${data_dir}${classT}${csv_loc}${ch_org}-data.csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(echo \$(bash --version | grep -Eo 'version [[:alnum:].]+' | sed 's/version //'))
    END_VERSIONS
    """
}
