process PRETEXTANNOTATE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3c/3c63c935c52e8af4e4ba9ac16b4ade1ea5f176a373464faec9ab3f9feb3f7ebc/data'
        : 'community.wave.seqera.io/library/pip_pretextannotate:1fe480fd3b9b3943'}"

    input:
    tuple val(meta),  path(sizes)
    tuple val(meta2), path(snapshot)

    output:
    tuple val(meta), path("*.png"), emit: png
    tuple val(meta), path("*.gif"), emit: gif
    tuple val(meta), path("*.tif"), emit: tif
    tuple val("${task.process}"), val('pretextannotate'), eval("pretextannotate -v | sed 's/pretextannotate. //g'"), topic: versions, emit: versions_pretextannotate

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args     ?: ""
    def prefix  = task.ext.prefix   ?: "${meta.id}"

    """
    pretextannotate \\
        --pretext_file ${snapshot} \\
        --sizes ${sizes} \\
        --output ./ \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.png
    touch ${prefix}.gif
    touch ${prefix}.tif
    """
}
