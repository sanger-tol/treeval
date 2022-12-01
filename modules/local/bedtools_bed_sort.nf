process BEDTOOLS_BED_SORT {
    tag "${meta.id}"
    label "process_medium"

    def version = '2.30.0--hc088bd4_0'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the chunkfasta process. Please use docker or singularity containers."
    }
    container "quay.io/repository/biocontainers/bedtools:${version}"

    input:
    tuple val( meta ), path( merged_bam )

    output:
    tuple val( meta ), file( "*.bed" ), emit: sorted_bed
    path "versions.yml",                emit: versions

    script:
    """
    bamToBed \\
        -i $merged_bam \\
        -bed12 | \\
        bedtools sort \\
        -i > $meta.id/.sorted.bed
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        BEDTOOLS_BED_SORT: $version
    END_VERSIONS
    """
}