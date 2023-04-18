process PAF2BED {
    tag "${meta.id}"
    label "process_low"

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'ubuntu:20.04' }"

    input:
    tuple val( meta ), path( file )

    output:
    tuple val( meta ), file( "*_punchlist.bed" ), emit: punchlist

    script:
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    paf_to_bed12.sh ${file} ${meta.id}_${meta.type}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paf_to_bed12:   \$(paf_to_bed12.sh -v)
        coreutils:      $VERSION
    END_VERSIONS
    """

    stub:
    def VERSION = "9.1"
    """
    touch ${meta.id}_${meta.type}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paf_to_bed12:   \$(paf_to_bed12.sh -v)
        coreutils:      $VERSION
    END_VERSIONS
    """
}