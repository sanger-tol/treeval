process GET_READS_FROM_DIRECTORY {
    tag "${meta.id}"
    label "process_low"

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'ubuntu:20.04' }"

    input:
    tuple val(meta), val(directory_path)

    output:
    tuple val(meta), path("*.fasta.gz"), emit: files

    script:
    """
    if [ ! -d ${directory_path} ] || [ -z "\$(ls -A ${directory_path})" ]
    then
        echo "Directory is empty or doesn't exist"
    else
        cp ${directory_path}*.fasta.gz .
    fi
    """
}
