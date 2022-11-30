process CONCAT_GFF {
    tag "${meta.id} - ${meta.type}"
    label "process_medium"

    input:
    tuple val(meta), file(input_files)

    output:
    tuple val(meta), file ("${meta.id}-${meta.type}-all.gff"),     emit: concat_gff

    script:
    """
    cat $input_files > $meta.id-$meta.type-all.gff
    """
}