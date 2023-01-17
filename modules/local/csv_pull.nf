process CSV_PULL {
    tag "${ch_org}"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::coreutils=9.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'ubuntu:20.04' }"

    input:
    val ch_org
    path csv_loc

    output:
    path "${ch_org}-data.csv",      emit: csv

    script:
    """
    echo "File downloaded from: ${csv_loc}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(echo \$(bash --version | grep -Eo 'version [[:alnum:].]+' | sed 's/version //'))
    END_VERSIONS
    """
}
