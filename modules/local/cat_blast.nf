process CAT_BLAST {
    tag "${meta.id} - ${meta.type}"
    label "process_medium"

    input:
    tuple val(meta), file(input_files)

    output:
    tuple val(meta), file ("${meta.id}-${meta.type}-final_blast.tsv"),     emit: concat_blast

    script:
    """
    cat $input_files > $meta.id-$meta.type-final_blast.tsv
    """
}
