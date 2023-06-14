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
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    $/
    cat "${mergeblocks}" \
    |awk '{split($4,a,":");print $1"\t"$2"\t"$3"\t"a[1]"\t"$5"\t"$6}'\
    |awk 'sqrt(($3-$2)*($3-$2)) > 5000'\
    |sort -k 1,1 -k2,2n \
    > ${prefix}_chain.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    /$

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}_chain.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    """
}
