process GFF_TO_BED {
    tag "${meta.id}"
    label "process_low"

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'ubuntu:20.04' }"

    input:
    tuple val( meta ), path( file )

    output:
    tuple val( meta ), file( "*.bed" )  , emit: punchlist
    path "versions.yml"                 , emit: versions

    script:
    """
    gff_to_bed.sh ${file} ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gff_to_bed: \$(gff_to_bed.sh -v)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gff_to_bed: \$(gff_to_bed.sh -v)
    END_VERSIONS
    """
}
