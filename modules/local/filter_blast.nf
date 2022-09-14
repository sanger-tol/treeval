process FILTER_BLAST {
    tag "${meta.id} - ${meta.type}"
    label "process_medium"

    def version = '0.002-c1'

    container "quay.io/sanger-tol/genealignment:${version}"

    input:
    tuple val( meta ), file( concat_blast_out )

    output:
    tuple val( meta ), file( "${meta.id}-${meta.type}-*.tsv")   , emit: final_tsv
    path "versions.yml"                                         , emit: versions
    path "*log"

    script:
    def filt_percent = task.ext.args ?: 90.00
    """
    filter_blast.py $meta.id $meta.type $concat_blast_out $filt_percent
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        filter_blast: $version
    END_VERSIONS
    """
}
