process MAKECMAP_CMAP2BED {
    tag "$meta.id"
    label 'process_medium'

    def version = '0.001-c1'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the makecmap process. Please use docker or singularity containers."
    }
    container "quay.io/sanger-tol/cmap2bed:${version}"


    input:
    tuple val(meta), path(cmap)
    val enzyme

    output:
    tuple val(meta), path("*.bed"), emit: bedfile
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    grep -v '#' $cmap > ${prefix}_${enzyme}_edited.cmap
    /scripts/wrapper.sh -t ${prefix}_${enzyme}_edited.cmap -z $enzyme | sort -k1,1 -k2,2n > ${prefix}_${enzyme}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cmap2bed: $version
    END_VERSIONS
    """
}
