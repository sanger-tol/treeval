process UCSC_BEDTOBIGBED {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::ucsc-bedtobigbed=377" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
<<<<<<< HEAD
        'https://depot.galaxyproject.org/singularity/ucsc-bedtobigbed:377--ha8a8165_3' :
        'quay.io/biocontainers/ucsc-bedtobigbed:377--ha8a8165_3' }"
=======
        'https://depot.galaxyproject.org/singularity/ucsc-bedtobigbed:377--h446ed27_1' :
        'quay.io/biocontainers/ucsc-bedtobigbed:377--h446ed27_1' }"
>>>>>>> bb842cc (All modules added)

    input:
    tuple val(meta), path(bed)
    path  sizes
<<<<<<< HEAD
    path  autosql
=======
>>>>>>> bb842cc (All modules added)

    output:
    tuple val(meta), path("*.bigBed"), emit: bigbed
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
<<<<<<< HEAD
    def as_option = autosql ? "-as=${autosql}" : ""
=======
>>>>>>> bb842cc (All modules added)
    def VERSION = '377' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    bedToBigBed \\
        $bed \\
        $sizes \\
<<<<<<< HEAD
        $as_option \\
=======
>>>>>>> bb842cc (All modules added)
        $args \\
        ${prefix}.bigBed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """
}
