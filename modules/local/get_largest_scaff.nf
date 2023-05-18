process GET_LARGEST_SCAFF {

    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'ubuntu:20.04' }"

    input:
    tuple val( meta ), path( file )

    output:
    env largest_scaff       , emit: scaff_size

    script:
    """
    largest_scaff=`cut_size.sh $file`
    """
}