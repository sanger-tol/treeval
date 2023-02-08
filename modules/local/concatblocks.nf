process CONCATBLOCKS {
    tag "${meta.id}.bed"
    label "process_single"

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
            'ubuntu:20.04' }"

    input:
    tuple val(meta), path(mergeblocks)

    output:
    tuple val(meta), path("*.bed"), emit: chainfile
    path "versions.yml",              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat $mergeblocks | filter.sh > ${meta.id}_chain.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ubuntu: \$(ubuntu --version | sed 's/Ubuntu //g')
    END_VERSIONS
    """
}
