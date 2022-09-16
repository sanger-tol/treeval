process SELFCOMP_MUMMER2BED {
    tag "$meta.id"
    label 'process_medium'
    def version = '0.001-c1'
    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the mummer2bed process. Please use docker or singularity containers."
    }
    container "quay.io/sanger-tol/selfcomp:${version}"

    input:
    tuple val(meta), path(mummerfile)
    val (motiflen)

    output:
    tuple val(meta), path("*.bed"), emit: bedfile
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mummer2bed.py $args -i $mummerfile -l $motiflen > ${prefix}.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        selfcomp: $version
    END_VERSIONS
    """
}
