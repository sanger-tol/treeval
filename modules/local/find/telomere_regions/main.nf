process FIND_TELOMERE_REGIONS {
    tag "${meta.id}"
    label 'process_low'

    container 'quay.io/sanger-tol/telomere:0.0.1-c2'

    input:
    tuple val( meta ), path( file )
    val (telomereseq)

    output:
    tuple val( meta ), file( "*.telomere" ) , emit: telomere
    path "versions.yml"                     , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    find_telomere ${file} $telomereseq > ${prefix}.telomere

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        find_telomere: \$(find_telomere -V)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.telomere

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        find_telomere: \$(find_telomere -V)
    END_VERSIONS
    """

}
