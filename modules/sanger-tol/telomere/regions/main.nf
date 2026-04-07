process TELOMERE_REGIONS {
    tag "${meta.id}"
    label 'process_low'

    container 'sanger-tol/telomere:0.0.1-c2'

    input:
    tuple val(meta), path(reference)
    val telomereseq

    output:
    tuple val( meta ), path( "*.telomere" ) , emit: telomere
    tuple val("${task.process}"), val('find_telomere_regions'), val("1.0.0"), topic: versions, emit: versions_telomereregions

    when:
    task.ext.when == null || task.ext.when

    script:

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "TELOMERE_REGIONS module does not support Conda. Please use Docker / Singularity instead."
    }

    def prefix          = task.ext.prefix ?: "${meta.id}"
    """
    find_telomere $reference $telomereseq > ${prefix}.telomere
    """

    stub:

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "TELOMERE_REGIONS module does not support Conda. Please use Docker / Singularity instead."
    }

    def prefix          = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.telomere
    """

}
