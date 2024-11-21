process PAF2BED {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'docker.io/ubuntu:20.04' }"

    input:
    tuple val( meta ), path( file )

    output:
    tuple val( meta ), file( "*_punchlist.bed" ), emit: punchlist
    path "versions.yml"                         , emit: versions

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}_${meta.type}_punchlist"
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    paf_to_bed.sh ${file} ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paf_to_bed:   \$(paf_to_bed.sh -v)
        coreutils:    $VERSION
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}_${meta.type}_punchlist"
    def VERSION = "9.1"  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paf_to_bed:   \$(paf_to_bed.sh -v)
        coreutils:    $VERSION
    END_VERSIONS
    """
}
