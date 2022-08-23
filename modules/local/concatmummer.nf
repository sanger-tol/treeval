process CONCATMUMMER {
    tag "${meta} -> ${meta.id}.mummer"
    label "process_small"

    input:
    tuple val(meta), path(coords)

    output:
    tuple val(meta), path("*.mummer"), emit: mummer

    script:
    """
    cat $coords > ${meta.id}.mummer
    """
}
