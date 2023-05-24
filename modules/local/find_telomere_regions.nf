process FIND_TELOMERE_REGIONS {
    tag "${meta.id}"
    label "process_low"

    container 'docker.io/library/gcc:7.1.0'

    input:
    tuple val( meta ), path( file )
    val (telomereseq)

    output:
    tuple val( meta ), file( "*.telomere" ) , emit: telomere
    path "versions.yml"                     , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.0" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def find_telomere = task.ext.find_telomere ?: ''
    """
    find_telomere ${file} $telomereseq > ${prefix}.telomere

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        find_telomere: \$(find_telomere -v)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.0" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def find_telomere = task.ext.find_telomere ?: ''
    """
    touch ${prefix}.telomere

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        find_telomere: \$(find_telomere -v)
    END_VERSIONS
    """

}
